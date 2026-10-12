//! Residual equations (design note 0005 §4) with analytic gradients.
//!
//! Unknowns live in a flat slot vector: point `p` has slots `2p` (x) and
//! `2p + 1` (y); circle radii, conic shapes and hidden curve parameters
//! follow. Every
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
    /// `E(θ)_axis − p_axis` for a point on an ellipse (design note 0008
    /// §4.2): `E(θ) = (1 − cos θ − sin θ)·c + cos θ·m + sin θ·n`.
    EllipsePoint {
        /// Centre point.
        c: usize,
        /// Major axis end point.
        m: usize,
        /// Minor axis end point.
        n: usize,
        /// Slot of the hidden parameter `θ`.
        t: usize,
        /// The point.
        p: usize,
        /// 0 for x, 1 for y.
        axis: usize,
    },
    /// A line tangent to an ellipse (design note 0008 §4.3):
    /// `cross(dir, v)/‖dir‖ − sign·√((n̂·a)² + (n̂·b)²)`, `n̂` the line's unit
    /// left normal, `v` the centre minus a point of the line and `a`, `b`
    /// the semi-axes.
    EllipseSupport {
        /// Line direction.
        dir: VecExpr,
        /// Centre minus the line's start.
        v: VecExpr,
        /// First semi-axis.
        a: VecExpr,
        /// Second semi-axis.
        b: VecExpr,
        /// +1 when the centre is on the line's left, −1 otherwise.
        sign: f64,
    },
    /// Directions compared at a shared point when one of them (or both) is
    /// an elliptical arc's tangent there (design note 0008 §4.4); otherwise
    /// as [`EqKind::DirPair`] with [`PairMode::Cross`] or [`PairMode::Dot`].
    TanDirs {
        /// First direction.
        a: TanDir,
        /// Second direction.
        b: TanDir,
        /// Cross (parallel) or dot (perpendicular).
        mode: PairMode,
    },
    /// `C(t)_axis − p_axis` for a point on a conic arc (design note 0009
    /// §4.2): `C(t) = Σ Rₖ(t) Pₖ`, the rational quadratic Bézier with
    /// middle weight `e^u`.
    ConicPoint {
        /// Start point.
        p0: usize,
        /// Control point.
        p1: usize,
        /// End point.
        p2: usize,
        /// Slot of the shape unknown `u = ln w`.
        u: usize,
        /// Slot of the hidden parameter `t ∈ [0, 1]`.
        t: usize,
        /// The point.
        p: usize,
        /// 0 for x, 1 for y.
        axis: usize,
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

/// A direction for [`EqKind::TanDirs`].
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum TanDir {
    /// A linear combination of points.
    Lin(VecExpr),
    /// The tangent of the ellipse on centre `c` and axis ends `m`, `n` at
    /// its point `p`, written without a parameter (design note 0008 §4.4):
    /// with `w = p − c`, `A = m − c`, `B = n − c`,
    /// `T = (−(w·B)/‖B‖²)·A + ((w·A)/‖A‖²)·B`, which is `E′(θ)` at `p`'s
    /// parameter when `p` is on the ellipse and `A ⟂ B`.
    Ellipse {
        /// Centre.
        c: usize,
        /// Major axis end.
        m: usize,
        /// Minor axis end.
        n: usize,
        /// The point.
        p: usize,
    },
}

impl TanDir {
    /// The direction's value, its Jacobian columns `(slot, ∂x/∂slot,
    /// ∂y/∂slot)` if `jac` is given, or `None` where an axis has no length.
    fn eval(&self, v: &[f64], jac: Option<&mut Vec<(usize, f64, f64)>>) -> Option<(f64, f64)> {
        match self {
            TanDir::Lin(e) => {
                if let Some(j) = jac {
                    for &(p, c) in &e.0 {
                        j.push((2 * p, c, 0.0));
                        j.push((2 * p + 1, 0.0, c));
                    }
                }
                Some(e.eval(v))
            }
            &TanDir::Ellipse { c, m, n, p } => {
                let pt = |q: usize| [v[2 * q], v[2 * q + 1]];
                let (pc, pm, pn, pp) = (pt(c), pt(m), pt(n), pt(p));
                let w = [pp[0] - pc[0], pp[1] - pc[1]];
                let a = [pm[0] - pc[0], pm[1] - pc[1]];
                let b = [pn[0] - pc[0], pn[1] - pc[1]];
                let a2 = a[0] * a[0] + a[1] * a[1];
                let b2 = b[0] * b[0] + b[1] * b[1];
                if !(a2 > 0.0 && b2 > 0.0) || !a2.is_finite() || !b2.is_finite() {
                    return None;
                }
                // T = −σ·A + κ·B with σ = (w·B)/b2, κ = (w·A)/a2.
                let sg = (w[0] * b[0] + w[1] * b[1]) / b2;
                let ka = (w[0] * a[0] + w[1] * a[1]) / a2;
                let t = [-sg * a[0] + ka * b[0], -sg * a[1] + ka * b[1]];
                if let Some(j) = jac {
                    // ∂σ/∂w = B/b2, ∂σ/∂B = w/b2 − 2σB/b2;
                    // ∂κ/∂w = A/a2, ∂κ/∂A = w/a2 − 2κA/a2.
                    // dT = −dσ·A − σ·dA + dκ·B + κ·dB.
                    // Columns over (w, A, B) components, row i of T.
                    let mut cols: [[f64; 2]; 6] = [[0.0; 2]; 6];
                    for i in 0..2 {
                        for k in 0..2 {
                            let delta = if i == k { 1.0 } else { 0.0 };
                            // w_k
                            cols[k][i] = -(b[k] / b2) * a[i] + (a[k] / a2) * b[i];
                            // A_k
                            cols[2 + k][i] = -sg * delta + ((w[k] - 2.0 * ka * a[k]) / a2) * b[i];
                            // B_k
                            cols[4 + k][i] = -((w[k] - 2.0 * sg * b[k]) / b2) * a[i] + ka * delta;
                        }
                    }
                    for k in 0..2 {
                        let (cw, ca, cb) = (cols[k], cols[2 + k], cols[4 + k]);
                        j.push((2 * p + k, cw[0], cw[1]));
                        j.push((2 * m + k, ca[0], ca[1]));
                        j.push((2 * n + k, cb[0], cb[1]));
                        // w, A and B all subtract the centre.
                        j.push((
                            2 * c + k,
                            -(cw[0] + ca[0] + cb[0]),
                            -(cw[1] + ca[1] + cb[1]),
                        ));
                    }
                }
                Some((t[0], t[1]))
            }
        }
    }

    /// The points referenced.
    fn points(&self) -> Vec<usize> {
        match self {
            TanDir::Lin(e) => e.points().collect(),
            &TanDir::Ellipse { c, m, n, p } => vec![c, m, n, p],
        }
    }
}

/// The ellipse point `E(θ)_axis` and `∂/∂θ` from slot values.
fn ellipse_axis(v: &[f64], c: usize, m: usize, n: usize, t: usize, axis: usize) -> (f64, f64) {
    let (s, co) = v[t].sin_cos();
    let (xc, xm, xn) = (v[2 * c + axis], v[2 * m + axis], v[2 * n + axis]);
    (
        xc + co * (xm - xc) + s * (xn - xc),
        -s * (xm - xc) + co * (xn - xc),
    )
}

/// A spline as used by equations: its linear map and its points.
#[derive(Clone, Debug)]
pub(crate) struct SplineDef {
    /// The map from points to curve.
    pub map: SplineMap,
    /// Point indices, in the spline's order.
    pub points: Vec<usize>,
}

impl SplineDef {
    /// The point of the curve at `t` (slot values `v`).
    pub fn point_at(&self, v: &[f64], t: f64) -> (f64, f64) {
        let (w, _) = self.map.weights(t);
        let (mut x, mut y) = (0.0, 0.0);
        for (k, &q) in self.points.iter().enumerate() {
            x += w[k] * v[2 * q];
            y += w[k] * v[2 * q + 1];
        }
        (x, y)
    }

    /// The squared distance from point `p` to the curve point at `t`.
    pub fn dist2(&self, v: &[f64], p: usize, t: f64) -> f64 {
        let (x, y) = self.point_at(v, t);
        let (dx, dy) = (x - v[2 * p], y - v[2 * p + 1]);
        dx * dx + dy * dy
    }

    /// The curve parameter of the point closest to point `p` (slot values
    /// `v`). Every knot span (polynomial piece) is sampled; each local
    /// minimum of the sampled distance is refined by golden-section search
    /// between its neighbouring samples, then by safeguarded Newton on
    /// `(C(t) − P)·C′(t) = 0` (P&T eq. 6.3); the nearest wins (the first on
    /// ties). A curve that passes near `P` more than once is thus not
    /// mistaken for its nearest sample's neighbourhood.
    pub fn closest_param(&self, v: &[f64], p: usize) -> f64 {
        const PER_SPAN: usize = 16;
        let mut ts: Vec<f64> = Vec::new();
        for w in self.map.knots.windows(2) {
            if w[1] > w[0] {
                for i in 0..PER_SPAN {
                    ts.push(w[0] + (w[1] - w[0]) * (i as f64 / PER_SPAN as f64));
                }
            }
        }
        ts.push(1.0);
        let d: Vec<f64> = ts.iter().map(|&t| self.dist2(v, p, t)).collect();
        let n = ts.len();
        let mut best = (f64::INFINITY, 0.0);
        for i in 0..n {
            let left = if i > 0 { d[i - 1] } else { f64::INFINITY };
            let right = if i + 1 < n { d[i + 1] } else { f64::INFINITY };
            // Local minima of the samples; of a run of equal samples (a
            // stretch of curve that does not move, such as a zero-length
            // leg), only its two ends.
            if !(d[i] <= left && d[i] <= right) || (d[i] == left && d[i] == right) {
                continue;
            }
            let (a, b) = (ts[i.saturating_sub(1)], ts[(i + 1).min(n - 1)]);
            let (t, dt) = self.refine(v, p, (ts[i], d[i]), a, b);
            if dt < best.0 {
                best = (dt, t);
            }
        }
        best.1
    }

    /// Refines a sampled local minimum `(t, distance²)` of the distance to
    /// point `p` within `[a, b]`; returns the improved `(t, distance²)`.
    fn refine(&self, v: &[f64], p: usize, start: (f64, f64), a: f64, b: f64) -> (f64, f64) {
        let (mut best_t, mut best_d) = start;
        // Golden-section search (Kiefer 1953) on [a, b].
        let g = 0.5 * (5f64.sqrt() - 1.0);
        let (mut lo, mut hi) = (a, b);
        let (mut x1, mut x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
        let (mut f1, mut f2) = (self.dist2(v, p, x1), self.dist2(v, p, x2));
        for _ in 0..40 {
            if f1 <= f2 {
                hi = x2;
                (x2, f2) = (x1, f1);
                x1 = hi - g * (hi - lo);
                f1 = self.dist2(v, p, x1);
            } else {
                lo = x1;
                (x1, f1) = (x2, f2);
                x2 = lo + g * (hi - lo);
                f2 = self.dist2(v, p, x2);
            }
        }
        for (t, d) in [(x1, f1), (x2, f2)] {
            if d < best_d {
                (best_t, best_d) = (t, d);
            }
        }
        // Newton polish, accepting steps that do not increase the distance
        // (near the minimum the distance stalls at rounding level while the
        // parameter still improves).
        for _ in 0..16 {
            let (w, w1, w2) = self.map.weights_ders(best_t, 2);
            let (mut c, mut c1, mut c2) = ([0.0f64; 2], [0.0f64; 2], [0.0f64; 2]);
            for (k, &q) in self.points.iter().enumerate() {
                for axis in 0..2 {
                    let x = v[2 * q + axis];
                    c[axis] += w[k] * x;
                    c1[axis] += w1[k] * x;
                    c2[axis] += w2[k] * x;
                }
            }
            c[0] -= v[2 * p];
            c[1] -= v[2 * p + 1];
            let gr = c[0] * c1[0] + c[1] * c1[1];
            let gp = c1[0] * c1[0] + c1[1] * c1[1] + c[0] * c2[0] + c[1] * c2[1];
            if gp.is_nan() || gp <= 0.0 {
                break;
            }
            let t = (best_t - gr / gp).clamp(0.0, 1.0);
            let d = self.dist2(v, p, t);
            if d.is_nan() || best_d.is_nan() || d > best_d || t == best_t {
                break;
            }
            (best_t, best_d) = (t, d);
        }
        (best_t, best_d)
    }
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
    /// The vectors whose length or direction the residual uses that are
    /// within `degen` of zero (or not finite) at slot values `v`: there
    /// the residual cannot be evaluated (a direction) or has a kink (a
    /// length, whose gradient is then an arbitrary choice).
    pub fn short_vectors(&self, v: &[f64], degen: f64) -> Vec<&VecExpr> {
        let short = |e: &VecExpr| {
            let (x, y) = e.eval(v);
            let n = norm(x, y);
            n.is_nan() || n <= degen
        };
        fn rhs_norms<'a>(r: &'a Rhs, out: &mut Vec<&'a VecExpr>) {
            for (_, t) in &r.terms {
                if let RTerm::Norm(e) = t {
                    out.push(e);
                }
            }
        }
        let mut out: Vec<&VecExpr> = Vec::new();
        match self {
            EqKind::Norm { v: e, rhs, .. } => {
                out.push(e);
                rhs_norms(rhs, &mut out);
            }
            EqKind::SignedDist { dir, rhs, .. } => {
                out.push(dir);
                rhs_norms(rhs, &mut out);
            }
            EqKind::Proj { dir, .. } => out.push(dir),
            EqKind::DirPair { a, b, .. } => {
                out.push(a);
                out.push(b);
            }
            EqKind::EllipseSupport { dir, a, b, .. } => {
                out.push(dir);
                out.push(a);
                out.push(b);
            }
            EqKind::TanDirs { a, b, .. } => {
                for d in [a, b] {
                    if let TanDir::Lin(e) = d {
                        out.push(e);
                    }
                }
            }
            EqKind::Lin { .. }
            | EqKind::AbsLin { .. }
            | EqKind::SplinePoint { .. }
            | EqKind::EllipsePoint { .. }
            | EqKind::ConicPoint { .. } => {}
        }
        out.retain(|e| short(e));
        out
    }

    /// True where an absolute value the residual takes is within `degen`
    /// of zero at slot values `v` (a kink).
    pub fn abs_kinked(&self, v: &[f64], degen: f64) -> bool {
        let near = |x: f64| x.abs() <= degen;
        match self {
            EqKind::AbsLin { terms, .. } => near(terms.iter().map(|&(s, c)| c * v[s]).sum()),
            EqKind::Norm { rhs, .. } | EqKind::SignedDist { rhs, .. } => {
                rhs.abs && near(rhs.raw(v))
            }
            _ => false,
        }
    }

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
            EqKind::SignedDist {
                dir,
                v: w,
                abs,
                rhs,
            } => {
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
                    dir.push_grad(sg * (wy / n - c * dx / n3), sg * (-wx / n - c * dy / n3), g);
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
            &EqKind::EllipsePoint {
                c,
                m,
                n,
                t,
                p,
                axis,
            } => {
                let (e, de) = ellipse_axis(v, c, m, n, t, axis);
                if let Some(g) = grad {
                    let (s, co) = v[t].sin_cos();
                    g.push((t, de));
                    g.push((2 * c + axis, 1.0 - co - s));
                    g.push((2 * m + axis, co));
                    g.push((2 * n + axis, s));
                    g.push((2 * p + axis, -1.0));
                }
                Some(Eval {
                    r: e - v[2 * p + axis],
                    lever: inf,
                })
            }
            &EqKind::ConicPoint {
                p0,
                p1,
                p2,
                u,
                t,
                p,
                axis,
            } => {
                let (w, tv) = (v[u].exp(), v[t]);
                let (x0, x1, x2) = (v[2 * p0 + axis], v[2 * p1 + axis], v[2 * p2 + axis]);
                let s = 1.0 - tv;
                let (b0, b1, b2) = (s * s, 2.0 * tv * s * w, tv * tv);
                let d = b0 + b1 + b2;
                let (r0, r1, r2) = (b0 / d, b1 / d, b2 / d);
                let c = r0 * x0 + r1 * x1 + r2 * x2;
                if let Some(g) = grad {
                    // ∂C/∂t = (N′ − C·D′)/D; ∂C/∂u = w·∂C/∂w = R₁(P₁ − C).
                    let n1 = -2.0 * s * x0 + 2.0 * (1.0 - 2.0 * tv) * w * x1 + 2.0 * tv * x2;
                    let d1 = -2.0 * s + 2.0 * (1.0 - 2.0 * tv) * w + 2.0 * tv;
                    g.push((t, (n1 - c * d1) / d));
                    g.push((u, r1 * (x1 - c)));
                    g.push((2 * p0 + axis, r0));
                    g.push((2 * p1 + axis, r1));
                    g.push((2 * p2 + axis, r2));
                    g.push((2 * p + axis, -1.0));
                }
                Some(Eval {
                    r: c - v[2 * p + axis],
                    lever: inf,
                })
            }
            EqKind::EllipseSupport {
                dir,
                v: w,
                a,
                b,
                sign,
            } => {
                let (dx, dy) = dir.eval(v);
                let nd = norm(dx, dy);
                if nd.is_nan() || nd <= degen {
                    return None;
                }
                // Unit left normal and unit direction.
                let (nx, ny) = (-dy / nd, dx / nd);
                let (ux, uy) = (dx / nd, dy / nd);
                let (wx, wy) = w.eval(v);
                let (ax, ay) = a.eval(v);
                let (bx, by) = b.eval(v);
                let pa = nx * ax + ny * ay;
                let pb = nx * bx + ny * by;
                let h = norm(pa, pb);
                if h.is_nan() || h <= degen {
                    return None;
                }
                let dist = nx * wx + ny * wy;
                if let Some(g) = grad {
                    // ∂(n̂·X)/∂d = ((X.y, −X.x) − (n̂·X)·û)/‖d‖.
                    let gd = |x: f64, y: f64, nxv: f64| ((y - nxv * ux) / nd, (-x - nxv * uy) / nd);
                    let (g1x, g1y) = gd(wx, wy, dist);
                    let (gax, gay) = gd(ax, ay, pa);
                    let (gbx, gby) = gd(bx, by, pb);
                    let (ka, kb) = (sign * pa / h, sign * pb / h);
                    w.push_grad(nx, ny, g);
                    a.push_grad(-ka * nx, -ka * ny, g);
                    b.push_grad(-kb * nx, -kb * ny, g);
                    dir.push_grad(g1x - ka * gax - kb * gbx, g1y - ka * gay - kb * gby, g);
                }
                Some(Eval {
                    r: dist - sign * h,
                    lever: inf,
                })
            }
            EqKind::TanDirs { a, b, mode } => {
                let mut ja = Vec::new();
                let mut jb = Vec::new();
                let want = grad.is_some();
                let (ax, ay) = a.eval(v, want.then_some(&mut ja))?;
                let (bx, by) = b.eval(v, want.then_some(&mut jb))?;
                let na = norm(ax, ay);
                let nb = norm(bx, by);
                if na.is_nan() || nb.is_nan() || na <= degen || nb <= degen {
                    return None;
                }
                let nn = na * nb;
                let cr = ax * by - ay * bx;
                let dt = ax * bx + ay * by;
                let (cc, cd) = match mode {
                    PairMode::Dot => (0.0, 1.0),
                    _ => (1.0, 0.0),
                };
                let r = (cc * cr + cd * dt) / nn;
                if let Some(g) = grad {
                    let (na2, nb2) = (na * na, nb * nb);
                    let gax = (cc * by + cd * bx) / nn - r * ax / na2;
                    let gay = (-cc * bx + cd * by) / nn - r * ay / na2;
                    let gbx = (-cc * ay + cd * ax) / nn - r * bx / nb2;
                    let gby = (cc * ax + cd * ay) / nn - r * by / nb2;
                    for &(s, jx, jy) in &ja {
                        g.push((s, gax * jx + gay * jy));
                    }
                    for &(s, jx, jy) in &jb {
                        g.push((s, gbx * jx + gby * jy));
                    }
                }
                Some(Eval {
                    r,
                    lever: na.min(nb),
                })
            }
            EqKind::SplinePoint { spline, t, p, axis } => {
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
            &EqKind::EllipsePoint {
                c,
                m,
                n,
                t,
                p,
                axis,
            } => {
                slots.push(t);
                for q in [c, m, n, p] {
                    slots.push(2 * q + axis);
                }
            }
            &EqKind::ConicPoint {
                p0,
                p1,
                p2,
                u,
                t,
                p,
                axis,
            } => {
                slots.push(t);
                slots.push(u);
                for q in [p0, p1, p2, p] {
                    slots.push(2 * q + axis);
                }
            }
            EqKind::EllipseSupport { dir, v, a, b, .. } => {
                for e in [dir, v, a, b] {
                    points.extend(e.points());
                }
            }
            EqKind::TanDirs { a, b, .. } => {
                points.extend(a.points());
                points.extend(b.points());
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

    #[test]
    fn closest_param_on_a_segment_and_a_polyline() {
        let def = SplineDef {
            map: crate::spline::SplineMap::control(2, 1).unwrap(),
            points: vec![0, 1],
        };
        // P0 = (0,0), P1 = (1,0), query (0.3, 0.5).
        let v = vec![0.0, 0.0, 1.0, 0.0, 0.3, 0.5];
        let t = def.closest_param(&v, 2);
        assert!((t - 0.3).abs() < 1e-12);
        // Beyond the end: clamped to it.
        let v = vec![0.0, 0.0, 1.0, 0.0, 1.7, 0.2];
        assert_eq!(def.closest_param(&v, 2), 1.0);
        // A polyline (0,0) → (1,0) → (1,1): the query (1.4, 0.6) is closest
        // to the second segment at (1, 0.6), t = 0.5 + 0.6 / 2.
        let def = SplineDef {
            map: crate::spline::SplineMap::control(3, 1).unwrap(),
            points: vec![0, 1, 2],
        };
        let v = vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.4, 0.6];
        let t = def.closest_param(&v, 3);
        assert!((t - 0.8).abs() < 1e-12, "{t}");
        assert!((def.dist2(&v, 3, t) - 0.16).abs() < 1e-12);
        // A zero-length first leg (control points 0, 0, 1, 2): the samples
        // over it tie; the query just past its end lies on the next leg.
        let def = SplineDef {
            map: crate::spline::SplineMap::control(4, 1).unwrap(),
            points: vec![0, 0, 1, 2],
        };
        let v = vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.01, 0.0];
        let t = def.closest_param(&v, 3);
        assert!((t - (1.0 + 0.01) / 3.0).abs() < 1e-12, "{t}");
        assert!(def.dist2(&v, 3, t) < 1e-24);
        // A curve passing the query twice, the nearer pass between samples
        // of a coarse grid: still found.
        let def = SplineDef {
            map: crate::spline::SplineMap::control(5, 1).unwrap(),
            points: vec![0, 1, 2, 3, 4],
        };
        // (0,0) → (10,0) → (10,0.001) → (0,0.001) → (0,0.5); query (5,0.0004).
        let v = vec![
            0.0, 0.0, 10.0, 0.0, 10.0, 0.001, 0.0, 0.001, 0.0, 0.5, 5.0, 0.0004,
        ];
        let t = def.closest_param(&v, 5);
        assert!(
            (def.dist2(&v, 5, t) - 0.0004f64.powi(2)).abs() < 1e-18,
            "{t}"
        );
        assert!((t - 0.125).abs() < 1e-12, "{t}");
    }

    /// Six points (12 slots), a radius slot (12), a radius slot (13), a
    /// spline parameter slot (14). A conic's shape uses slot 12 or 13.
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
            (
                "point on ellipse, x",
                EqKind::EllipsePoint {
                    c: 0,
                    m: 1,
                    n: 2,
                    t: 14,
                    p: 3,
                    axis: 0,
                },
            ),
            (
                "point on ellipse, y",
                EqKind::EllipsePoint {
                    c: 4,
                    m: 1,
                    n: 5,
                    t: 14,
                    p: 0,
                    axis: 1,
                },
            ),
            (
                "point on conic, x",
                EqKind::ConicPoint {
                    p0: 0,
                    p1: 1,
                    p2: 2,
                    u: 12,
                    t: 14,
                    p: 3,
                    axis: 0,
                },
            ),
            (
                "point on conic, y",
                EqKind::ConicPoint {
                    p0: 4,
                    p1: 5,
                    p2: 1,
                    u: 13,
                    t: 14,
                    p: 0,
                    axis: 1,
                },
            ),
            (
                "line tangent to ellipse, centre left",
                EqKind::EllipseSupport {
                    dir: VecExpr::diff(5, 4),
                    v: VecExpr::diff(0, 4),
                    a: VecExpr::diff(1, 0),
                    b: VecExpr::diff(2, 0),
                    sign: 1.0,
                },
            ),
            (
                "line tangent to ellipse, centre right",
                EqKind::EllipseSupport {
                    dir: VecExpr::diff(3, 2),
                    v: VecExpr::diff(4, 2),
                    a: VecExpr::diff(5, 4),
                    b: VecExpr::diff(0, 4),
                    sign: -1.0,
                },
            ),
            (
                "elliptical arc tangent to a line at a point",
                EqKind::TanDirs {
                    a: TanDir::Ellipse {
                        c: 0,
                        m: 1,
                        n: 2,
                        p: 3,
                    },
                    b: TanDir::Lin(VecExpr::diff(5, 4)),
                    mode: PairMode::Cross,
                },
            ),
            (
                "elliptical arc tangent to an arc at a point",
                EqKind::TanDirs {
                    a: TanDir::Lin(VecExpr::diff(3, 5)),
                    b: TanDir::Ellipse {
                        c: 0,
                        m: 1,
                        n: 2,
                        p: 3,
                    },
                    mode: PairMode::Dot,
                },
            ),
            (
                "two elliptical arcs tangent at a point",
                EqKind::TanDirs {
                    a: TanDir::Ellipse {
                        c: 0,
                        m: 1,
                        n: 2,
                        p: 3,
                    },
                    b: TanDir::Ellipse {
                        c: 4,
                        m: 5,
                        n: 1,
                        p: 3,
                    },
                    mode: PairMode::Cross,
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
        let mut dense = [0.0; SLOTS];
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
                EqKind::SignedDist { .. }
                    | EqKind::Proj { .. }
                    | EqKind::DirPair { .. }
                    | EqKind::EllipseSupport { .. }
                    | EqKind::TanDirs { .. }
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
