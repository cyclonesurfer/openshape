//! Numerical solution of one component (design note 0005 §5 step 3, §6, §7).
//!
//! * **Gauss–Newton with minimum-norm steps.** `J Δ = −r` is solved in the
//!   weighted minimum-norm least-squares sense — among the least-squares
//!   solutions, the one minimizing `‖W^{1/2} Δ‖` — through the complete
//!   orthogonal decomposition of `J W^{−1/2}` (Golub & Van Loan §5.5;
//!   Nocedal & Wright, *Numerical Optimization*, 2nd ed., §10.3). Unit
//!   weights make every step the smallest change that fixes the linearized
//!   residuals, so the solver stays on the branch nearest the start; dragged
//!   points carry weight 1e6 and move last.
//! * **Bounds and corners** (`bounded_step`): hidden curve parameters stay
//!   in `[0, 1]` by an active set; a radius may at most halve per step;
//!   at a corner of a degree-1 spline the step uses the one-sided
//!   derivative of the piece it moves into (both pieces are tried), and
//!   the line search stops at the first corner crossed.
//! * **Armijo backtracking** on `½‖r‖²` (c = 1e-4, halving down to 2⁻²⁰;
//!   Nocedal & Wright §3.1). Trial points must be admissible: finite, radii
//!   positive, and geometry within [`REACH`] sketch extents of the centre
//!   and the ±1e6 m coordinate limit.
//! * **Levenberg–Marquardt fallback** when the line search fails: steps from
//!   `(JᵀJ + μW) h = −Jᵀr`, solved as an augmented least-squares problem,
//!   with the gain-ratio update of μ of Madsen, Nielsen & Tingleff, *Methods
//!   for Non-Linear Least Squares Problems*, DTU 2004, §3.2 (eq. 3.14–3.16).
//! * **Convergence** per equation: `|rᵢ| ≤ tolᵢ`, with `tolᵢ = 1e-12`
//!   (scaled) for length residuals and, for angle residuals, the angle whose
//!   lever (the shorter direction involved) moves a point by 1e-12 of the
//!   sketch size, clamped to [1e-12, 1e-10]; never below the rounding noise
//!   of the residual at the current point, never above 1e-10.
//! * **Restarts.** When an iteration stalls, hidden curve parameters are
//!   moved to the foot points of their points (`CompSystem::reproject`)
//!   and the iteration continues. If the direct solve fails, it is retried
//!   with steps limited to turning any direction by [`MAX_TURN`] (the first
//!   attempt may have swung a line through and into another branch's
//!   basin).
//! * **Continuation** (§6): if both fail (or land on the
//!   antipodal branch of an angle), the change from the start is applied in
//!   equal steps of 1/2, halving the step on failure down to 1/64, each
//!   solved from the previous result. Dimensions and fixed values are
//!   interpolated from their current values; other equations follow the
//!   Newton homotopy `r(x) − (1 − s)·r(x₀)`.

use openshape_linalg::{dot, lstsq_min_norm, weighted_min_norm, Matrix};

use crate::build::Problem;
use crate::eq::{norm, EqKind, PairMode, VecExpr};
use crate::reduce::{Component, Reduced};

/// Relative rank threshold for steps and diagnostics (design note §5).
pub(crate) const RANK_TOL: f64 = 1e-10;
/// Upper bound on any equation's convergence tolerance (scaled).
const TOL_CAP: f64 = 1e-10;
/// Base tolerance for angle residuals.
const TOL_ANGLE: f64 = 1e-12;

/// The largest rotation of any direction an equation uses (a line's
/// direction, a radius vector, the vector between two dimensioned points)
/// in one step. The equations are linearized in the coordinates, but they
/// depend on these directions nonlinearly; a step that swings a direction
/// through a large angle (or flips a line by moving one end past the
/// other) leaves the region where the linearization means anything, and
/// can land in the basin of a solution on another branch, or of none. The
/// line search treats such steps like inadmissible ones (a trust region on
/// the directions; cf. Nocedal & Wright §4).
pub(crate) const MAX_TURN: f64 = std::f64::consts::FRAC_PI_4;

/// Foot-point restarts allowed per Gauss–Newton run (see
/// `CompSystem::reproject`).
const MAX_REPROJECTIONS: usize = 4;

/// How far the solver may move geometry: within this many times the
/// sketch's extent (the scaling length) of the sketch's centre.
pub(crate) const REACH: f64 = 1e3;

/// The admissible values of coordinate and radius unknowns (scaled): the
/// solver never leaves this box, so it can never return geometry that input
/// validation would reject (coordinates beyond ±1e6 m) or that has run off
/// along a degenerate direction (for example, control points growing
/// without bound as a point slides towards the clamped end of a spline,
/// where their weights vanish).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Region {
    /// Bounds of x coordinates.
    pub x: (f64, f64),
    /// Bounds of y coordinates.
    pub y: (f64, f64),
    /// Largest radius.
    pub r_max: f64,
}

/// Tolerances for length residuals (scaled), chosen per sketch so that the
/// unscaled residual also meets the model tolerance.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LengthTol {
    /// Convergence target.
    pub base: f64,
    /// Never exceeded, whatever the rounding noise.
    pub cap: f64,
}

/// The equations of one component over its local unknowns.
pub(crate) struct CompSystem<'a> {
    prob: &'a Problem,
    /// Equation kinds (angle values may be overridden during continuation).
    pub kinds: Vec<EqKind>,
    angular: Vec<bool>,
    /// Subtracted from each residual (continuation).
    offsets: Vec<f64>,
    /// Slots of each local unknown.
    members: Vec<Vec<usize>>,
    /// Slot → local unknown, or `usize::MAX`.
    col_of_slot: Vec<usize>,
    /// Metric weights of the local unknowns.
    pub weights: Vec<f64>,
    is_param: Vec<bool>,
    is_radius: Vec<bool>,
    /// Admissible range of each local unknown.
    range: Vec<(f64, f64)>,
    /// For hidden parameters, the spline they run along.
    param_spline: Vec<Option<usize>>,
    /// For hidden parameters, the point they place on the spline.
    param_point: Vec<Option<usize>>,
    /// The vectors whose direction the equations use (line directions,
    /// radius vectors, point-to-point vectors), for the turn limit.
    vectors: Vec<VecExpr>,
    /// Whether steps are limited to turning directions by [`MAX_TURN`].
    turn_limit: bool,
    /// Slot values: constants and the current unknowns.
    pub v: Vec<f64>,
    len_tol: LengthTol,
    grad: Vec<(usize, f64)>,
}

/// Result of one Gauss–Newton run.
#[derive(Clone, Debug)]
pub(crate) struct Outcome {
    /// Final unknowns.
    pub u: Vec<f64>,
    /// Converged (every residual within tolerance, angles on their branch).
    pub converged: bool,
    /// Iterations used.
    pub iterations: usize,
}

impl<'a> CompSystem<'a> {
    /// The system of `comp`, with slot values `v` (constants included) and
    /// per-unknown weights `unknown_weight` (indexed by reduced unknown).
    pub fn new(
        prob: &'a Problem,
        red: &Reduced,
        comp: &Component,
        v: Vec<f64>,
        unknown_weight: &[f64],
        len_tol: LengthTol,
        region: Region,
    ) -> Self {
        let mut col_of_slot = vec![usize::MAX; prob.n_slots()];
        let mut members = Vec::with_capacity(comp.unknowns.len());
        let mut is_param = Vec::with_capacity(comp.unknowns.len());
        let mut is_radius = Vec::with_capacity(comp.unknowns.len());
        let mut weights = Vec::with_capacity(comp.unknowns.len());
        let mut range = Vec::with_capacity(comp.unknowns.len());
        for (k, &u) in comp.unknowns.iter().enumerate() {
            let class = &red.classes[red.unknowns[u]];
            for &s in &class.members {
                col_of_slot[s] = k;
            }
            members.push(class.members.clone());
            let s0 = class.members[0];
            is_param.push(prob.is_param(s0));
            is_radius.push(prob.is_radius(s0));
            range.push(if prob.is_param(s0) {
                (0.0, 1.0)
            } else if prob.is_radius(s0) {
                (0.0, region.r_max)
            } else if s0 % 2 == 0 {
                region.x
            } else {
                region.y
            });
            weights.push(unknown_weight[u]);
        }
        let mut param_spline = vec![None; comp.unknowns.len()];
        let mut param_point = vec![None; comp.unknowns.len()];
        for &i in &comp.eqs {
            if let EqKind::SplinePoint { spline, t, p, .. } = prob.eqs[i].kind {
                if col_of_slot[t] != usize::MAX {
                    param_spline[col_of_slot[t]] = Some(spline);
                    param_point[col_of_slot[t]] = Some(p);
                }
            }
        }
        let mut vectors = Vec::new();
        let mut add = |e: &VecExpr| {
            // Only vectors between two points: a spline's end leg is a
            // combination of many points, and its turning is the curve
            // changing shape, not a branch change.
            if e.0.len() == 2 {
                vectors.push(e.clone());
            }
        };
        for &i in &comp.eqs {
            match &prob.eqs[i].kind {
                EqKind::Norm { v, .. } => add(v),
                EqKind::SignedDist { dir, .. } | EqKind::Proj { dir, .. } => add(dir),
                EqKind::DirPair { a, b, .. } => {
                    add(a);
                    add(b);
                }
                _ => {}
            }
        }
        CompSystem {
            prob,
            kinds: comp.eqs.iter().map(|&i| prob.eqs[i].kind.clone()).collect(),
            angular: comp.eqs.iter().map(|&i| prob.eqs[i].angular).collect(),
            offsets: vec![0.0; comp.eqs.len()],
            members,
            col_of_slot,
            weights,
            is_param,
            is_radius,
            range,
            param_spline,
            param_point,
            vectors,
            turn_limit: false,
            v,
            len_tol,
            grad: Vec::new(),
        }
    }

    /// Number of local unknowns.
    pub fn n(&self) -> usize {
        self.members.len()
    }

    /// Number of equations.
    pub fn m(&self) -> usize {
        self.kinds.len()
    }

    /// Writes `u` into the slot vector.
    pub fn write(&mut self, u: &[f64]) {
        for (k, slots) in self.members.iter().enumerate() {
            for &s in slots {
                self.v[s] = u[k];
            }
        }
    }

    /// Clamps hidden parameters to `[0, 1]`; returns true if any moved.
    pub fn project(&self, u: &mut [f64]) -> bool {
        let mut moved = false;
        for (k, x) in u.iter_mut().enumerate() {
            if self.is_param[k] {
                let c = x.clamp(0.0, 1.0);
                if c != *x {
                    *x = c;
                    moved = true;
                }
            }
        }
        moved
    }

    fn admissible(&self, u: &[f64]) -> bool {
        u.iter()
            .zip(&self.is_radius)
            .zip(&self.range)
            .all(|((x, &rad), &(lo, hi))| {
                x.is_finite() && *x >= lo && *x <= hi && (!rad || *x > self.prob.degen)
            })
    }

    /// Residuals at `u`, or `None` if `u` is not admissible (degenerate
    /// direction, non-positive radius, non-finite value).
    pub fn residuals(&mut self, u: &[f64]) -> Option<Vec<f64>> {
        if !self.admissible(u) {
            return None;
        }
        self.write(u);
        let mut r = Vec::with_capacity(self.m());
        for (i, k) in self.kinds.iter().enumerate() {
            let e = k.eval(&self.v, &self.prob.splines, self.prob.degen, None)?;
            let ri = e.r - self.offsets[i];
            if !ri.is_finite() {
                return None;
            }
            r.push(ri);
        }
        Some(r)
    }

    /// Residuals, Jacobian (`m × n`) and per-equation tolerances at `u`.
    pub fn jacobian(&mut self, u: &[f64]) -> Option<(Vec<f64>, Matrix, Vec<f64>)> {
        if !self.admissible(u) {
            return None;
        }
        self.write(u);
        let (m, n) = (self.m(), self.n());
        let mut r = Vec::with_capacity(m);
        let mut tol = Vec::with_capacity(m);
        let mut j = Matrix::zeros(m, n);
        let mut grad = std::mem::take(&mut self.grad);
        for i in 0..m {
            grad.clear();
            let e = self.kinds[i].eval(
                &self.v,
                &self.prob.splines,
                self.prob.degen,
                Some(&mut grad),
            )?;
            let ri = e.r - self.offsets[i];
            if !ri.is_finite() {
                self.grad = grad;
                return None;
            }
            let mut noise = 0.0;
            let mut gross = 0.0;
            for &(s, g) in &grad {
                noise += g.abs() * self.v[s].abs();
                gross += g.abs();
                let c = self.col_of_slot[s];
                if c != usize::MAX {
                    j[(i, c)] += g;
                }
            }
            // A row that cancels to within RANK_TOL of its own terms (the
            // members of an equality class pulling equally both ways, or a
            // derivative that is zero up to rounding, such as an angle's
            // with respect to a vertical line's length) is zero: alone in
            // a component, its rounding noise would otherwise count as
            // full rank against itself.
            let row: f64 = (0..n).map(|c| j[(i, c)].abs()).sum();
            if row <= RANK_TOL * gross {
                for c in 0..n {
                    j[(i, c)] = 0.0;
                }
            }
            r.push(ri);
            tol.push(equation_tol(self.angular[i], e.lever, noise, self.len_tol));
        }
        self.grad = grad;
        if !j.is_finite() {
            return None;
        }
        Some((r, j, tol))
    }

    /// For each local unknown that is a hidden parameter sitting exactly on
    /// a kink of its spline (a corner of a degree-1 spline), the Jacobian
    /// column of the piece to the left of the kink; the column from
    /// [`Self::jacobian`] is the right-hand piece's. Uses the slot values
    /// last written.
    fn left_columns(&self) -> Vec<Option<Vec<f64>>> {
        let mut out: Vec<Option<Vec<f64>>> = vec![None; self.n()];
        for (k, sp) in self.param_spline.iter().enumerate() {
            let Some(sp) = *sp else { continue };
            let def = &self.prob.splines[sp];
            let t_slot = self.members[k][0];
            let Some(dw) = def.map.left_derivative(self.v[t_slot]) else {
                continue;
            };
            let mut col = vec![0.0; self.m()];
            for (i, kind) in self.kinds.iter().enumerate() {
                if let EqKind::SplinePoint { t, axis, .. } = kind {
                    if self.col_of_slot[*t] == k {
                        let mut dc = 0.0;
                        for (q, &pt) in def.points.iter().enumerate() {
                            dc += dw[q] * self.v[2 * pt + axis];
                        }
                        col[i] = dc;
                    }
                }
            }
            out[k] = Some(col);
        }
        out
    }

    /// Moves each hidden parameter to the closest point of its curve (the
    /// foot point) when that is strictly closer than where it is. Returns
    /// true if any moved.
    ///
    /// The parameter has no meaning of its own: any value that puts the
    /// point on the curve is a solution. When Gauss–Newton stalls, the
    /// foot point is the natural restart — the parameter correction of
    /// point-distance minimization in curve fitting (Wang, Pottmann & Liu,
    /// "Fitting B-spline curves to point clouds by curvature-based squared
    /// distance minimization", ACM TOG 25(2), 2006, §3) — and it can jump
    /// to another piece of the curve, past a corner or a fold, where the
    /// iteration cannot go.
    fn reproject(&mut self, u: &mut [f64]) -> bool {
        self.write(u);
        let mut moved = false;
        for (k, uk) in u.iter_mut().enumerate() {
            let (Some(sp), Some(p)) = (self.param_spline[k], self.param_point[k]) else {
                continue;
            };
            let def = &self.prob.splines[sp];
            let t = def.closest_param(&self.v, p);
            if def.dist2(&self.v, p, t) < def.dist2(&self.v, p, *uk) {
                *uk = t;
                moved = true;
            }
        }
        if moved {
            self.write(u);
        }
        moved
    }

    /// True if some hidden parameter sits exactly on a corner of its spline.
    fn at_corner(&self, u: &[f64]) -> bool {
        self.param_spline
            .iter()
            .zip(u)
            .any(|(sp, &t)| sp.is_some_and(|sp| self.prob.splines[sp].map.kinks.contains(&t)))
    }

    /// The first kink of unknown `k`'s spline strictly beyond `t` in the
    /// direction of motion.
    fn next_kink(&self, k: usize, t: f64, up: bool) -> Option<f64> {
        let kinks = &self.prob.splines[self.param_spline[k]?].map.kinks;
        if up {
            kinks.iter().copied().find(|&x| x > t)
        } else {
            kinks.iter().rev().copied().find(|&x| x < t)
        }
    }

    /// The current values of the direction vectors (see [`MAX_TURN`]).
    fn directions(&mut self, u: &[f64]) -> Vec<(f64, f64)> {
        self.write(u);
        self.vectors.iter().map(|e| e.eval(&self.v)).collect()
    }

    /// True if, at `u`, no direction vector has turned by more than
    /// [`MAX_TURN`] from `before` (vectors shorter than the degeneracy
    /// threshold at either end are not compared).
    fn turn_ok(&mut self, before: &[(f64, f64)], u: &[f64]) -> bool {
        if !self.turn_limit {
            return true;
        }
        self.write(u);
        let degen = self.prob.degen;
        self.vectors.iter().zip(before).all(|(e, &(x0, y0))| {
            let (x1, y1) = e.eval(&self.v);
            if norm(x0, y0) <= degen || norm(x1, y1) <= degen {
                return true;
            }
            let cross = x0 * y1 - y0 * x1;
            let dot = x0 * x1 + y0 * y1;
            cross.abs().atan2(dot) <= MAX_TURN
        })
    }

    /// True if every angle equation is on its intended branch:
    /// `cos(θ − v) > 0`.
    pub fn angles_ok(&mut self, u: &[f64]) -> bool {
        self.write(u);
        self.kinds.iter().all(|k| {
            k.angle_cos(&self.v, self.prob.degen)
                .is_none_or(|c| c > 0.0)
        })
    }

    /// Indices (local) of angle equations on the antipodal branch.
    pub fn wrong_branch(&mut self, u: &[f64]) -> Vec<usize> {
        self.write(u);
        let v = &self.v;
        self.kinds
            .iter()
            .enumerate()
            .filter(|(_, k)| k.angle_cos(v, self.prob.degen).is_some_and(|c| c <= 0.0))
            .map(|(i, _)| i)
            .collect()
    }
}

/// The convergence tolerance of one equation (see the module docs).
pub(crate) fn equation_tol(angular: bool, lever: f64, noise: f64, len: LengthTol) -> f64 {
    let floor = 32.0 * f64::EPSILON * noise;
    if angular {
        let base = if lever > 0.0 {
            (TOL_ANGLE / lever).clamp(TOL_ANGLE, TOL_CAP)
        } else {
            TOL_CAP
        };
        base.max(floor).min(TOL_CAP)
    } else {
        len.base.max(floor).min(len.cap)
    }
}

fn converged(r: &[f64], tol: &[f64]) -> bool {
    r.iter().zip(tol).all(|(ri, ti)| ri.abs() <= *ti)
}

/// The Gauss–Newton step under the unknowns' bounds and at spline corners,
/// by an active-set method in the manner of projected Newton methods
/// (Bertsekas, "Projected Newton methods for optimization problems with
/// simple constraints", SIAM J. Control Optim. 20(2), 1982):
///
/// * a hidden curve parameter at a bound of `[0, 1]` that the step would
///   push out is held there. Without this, a point sliding off the end of a
///   spline would be stuck: every step would put the along-curve correction
///   into the cheap parameter, and the projection would throw it away;
/// * a hidden parameter exactly at a corner of its spline uses the
///   derivative of the piece it moves into: the preferred piece's first;
///   if the step then moves it into the other piece, that piece's; if that
///   step moves it back, it is held at the corner (a generalized-Newton
///   choice of one-sided derivative for a piecewise-smooth residual; Qi &
///   Sun, "A nonsmooth version of Newton's method", Math. Programming 58,
///   1993). [`gauss_newton`] tries both preferences and keeps the better
///   step;
/// * a radius that the step would shrink below half its current value is
///   set to half (a fraction-to-the-boundary rule, Nocedal & Wright §19.2,
///   with the bound `r > 0`) and the rest of the step is recomputed. Without
///   this, the minimum-norm step can aim a small circle's radius below zero
///   on every iteration, and the line search then creeps towards a
///   collapsed circle instead of moving the other unknowns.
///
/// Returns the Jacobian of the directions taken (one-sided columns
/// substituted; for the line search's predicted decrease), the same with
/// the held or set columns zeroed (for the Levenberg–Marquardt fallback),
/// and the step.
fn bounded_step(
    sys: &mut CompSystem<'_>,
    u: &[f64],
    j: Matrix,
    neg: &[f64],
    prefer_left: bool,
) -> (Matrix, Matrix, Vec<f64>) {
    let n = j.cols();
    sys.write(u);
    let left = sys.left_columns();
    let right: Vec<Option<Vec<f64>>> = left
        .iter()
        .enumerate()
        .map(|(k, l)| l.as_ref().map(|_| j.col(k).to_vec()))
        .collect();
    let mut jd = j;
    // Side of each corner parameter: the preferred piece first, the other
    // after one switch.
    let mut on_left = vec![false; n];
    let mut switched = vec![false; n];
    if prefer_left {
        for (k, l) in left.iter().enumerate() {
            if let Some(col) = l {
                on_left[k] = true;
                jd.col_mut(k).copy_from_slice(col);
            }
        }
    }
    // Fixed step components (held parameters: 0; set radii: −r/2).
    let mut fixed: Vec<Option<f64>> = vec![None; n];
    loop {
        let mut jz = jd.clone();
        let mut rhs = neg.to_vec();
        for (k, f) in fixed.iter().enumerate() {
            if let Some(value) = *f {
                let col = jz.col_mut(k);
                for (b, c) in rhs.iter_mut().zip(col.iter()) {
                    *b -= c * value;
                }
                col.fill(0.0);
            }
        }
        let mut step = weighted_min_norm(&jz, &rhs, &sys.weights, RANK_TOL).x;
        for (x, f) in step.iter_mut().zip(&fixed) {
            if let Some(v) = f {
                *x = *v;
            }
        }
        let mut changed = false;
        for k in 0..n {
            if fixed[k].is_some() {
                continue;
            }
            if let (Some(lcol), Some(rcol)) = (&left[k], &right[k]) {
                let wrong_way = if on_left[k] {
                    step[k] > 0.0
                } else {
                    step[k] < 0.0
                };
                if wrong_way {
                    if switched[k] {
                        fixed[k] = Some(0.0);
                    } else {
                        switched[k] = true;
                        on_left[k] = !on_left[k];
                        jd.col_mut(k)
                            .copy_from_slice(if on_left[k] { lcol } else { rcol });
                    }
                    changed = true;
                    continue;
                }
            }
            let value = if sys.is_param[k]
                && ((u[k] <= 0.0 && step[k] < 0.0) || (u[k] >= 1.0 && step[k] > 0.0))
            {
                0.0
            } else if sys.is_radius[k] && u[k] + step[k] < 0.5 * u[k] {
                -0.5 * u[k]
            } else {
                continue;
            };
            fixed[k] = Some(value);
            changed = true;
        }
        if !changed {
            return (jd, jz, step);
        }
    }
}

/// The largest step fraction in `(0, 1]` that moves no hidden parameter
/// past a corner of its spline, and the corners reached at that fraction.
fn corner_limit(sys: &CompSystem<'_>, u: &[f64], step: &[f64]) -> (f64, Vec<(usize, f64)>) {
    let mut alpha = 1.0f64;
    let mut hits: Vec<(usize, f64, f64)> = Vec::new();
    for k in 0..u.len() {
        if !sys.is_param[k] || step[k] == 0.0 {
            continue;
        }
        let Some(kink) = sys.next_kink(k, u[k], step[k] > 0.0) else {
            continue;
        };
        let a = (kink - u[k]) / step[k];
        if a < 1.0 {
            hits.push((k, kink, a));
            alpha = alpha.min(a);
        }
    }
    let at: Vec<(usize, f64)> = hits
        .into_iter()
        .filter(|h| h.2 == alpha)
        .map(|h| (h.0, h.1))
        .collect();
    (alpha, at)
}

/// The Levenberg–Marquardt step `h = W^{−1/2} y` with
/// `[J W^{−1/2}; √μ I] y ≈ [−r; 0]`.
fn lm_step(j: &Matrix, r: &[f64], w: &[f64], mu: f64) -> Vec<f64> {
    let (m, n) = (j.rows(), j.cols());
    let s: Vec<f64> = w.iter().map(|x| 1.0 / x.sqrt()).collect();
    let sm = mu.sqrt();
    let a = Matrix::from_fn(m + n, n, |i, k| {
        if i < m {
            j[(i, k)] * s[k]
        } else if i - m == k {
            sm
        } else {
            0.0
        }
    });
    let mut b = vec![0.0; m + n];
    for i in 0..m {
        b[i] = -r[i];
    }
    let y = lstsq_min_norm(&a, &b, 0.0);
    y.iter().zip(&s).map(|(yi, si)| yi * si).collect()
}

/// Armijo backtracking along `step` from `u` (residuals `r`, `f = ½‖r‖²`,
/// Jacobian of the step's directions `j`): the first step fraction, from 1
/// (or the first spline corner) halving down to 2⁻²⁰, that is admissible,
/// within the turn limit, and decreases `f` sufficiently. Returns the
/// accepted point and its `f`.
fn armijo(
    sys: &mut CompSystem<'_>,
    u: &[f64],
    r: &[f64],
    f: f64,
    j: &Matrix,
    step: &[f64],
    dirs: &[(f64, f64)],
) -> Option<(Vec<f64>, f64)> {
    if !step.iter().all(|s| s.is_finite()) {
        return None;
    }
    let slope = dot(r, &j.mul_vec(step));
    if slope.is_nan() || slope >= 0.0 {
        return None;
    }
    // Stop at the first spline corner crossed: the linear model holds only
    // up to it.
    let (alpha0, corners) = corner_limit(sys, u, step);
    let mut alpha = alpha0;
    for _ in 0..=20 {
        let mut trial: Vec<f64> = u.iter().zip(step).map(|(a, d)| a + alpha * d).collect();
        let snapped = alpha == alpha0 && !corners.is_empty();
        if snapped {
            for &(k, kink) in &corners {
                trial[k] = kink;
            }
        }
        let clamped = sys.project(&mut trial) || snapped;
        if trial == u {
            return None;
        }
        if sys.turn_ok(dirs, &trial) {
            if let Some(rt) = sys.residuals(&trial) {
                let ft = 0.5 * dot(&rt, &rt);
                let pred = if clamped {
                    let disp: Vec<f64> = trial.iter().zip(u).map(|(a, b)| a - b).collect();
                    dot(r, &j.mul_vec(&disp))
                } else {
                    alpha * slope
                };
                if pred < 0.0 && ft <= f + 1e-4 * pred {
                    return Some((trial, ft));
                }
            }
        }
        alpha *= 0.5;
    }
    None
}

/// Damped Gauss–Newton from `u` (see the module docs).
pub(crate) fn gauss_newton(sys: &mut CompSystem<'_>, mut u: Vec<f64>, max_it: usize) -> Outcome {
    let mut its = 0usize;
    let mut mu: Option<f64> = None;
    let mut nu = 2.0f64;
    let mut stall = 0usize;
    let mut reprojections = 0usize;
    loop {
        let Some((r, j, tol)) = sys.jacobian(&u) else {
            return Outcome {
                u,
                converged: false,
                iterations: its,
            };
        };
        if converged(&r, &tol) {
            let ok = sys.angles_ok(&u);
            return Outcome {
                u,
                converged: ok,
                iterations: its,
            };
        }
        if its >= max_it {
            return Outcome {
                u,
                converged: false,
                iterations: its,
            };
        }
        its += 1;
        let f = 0.5 * dot(&r, &r);
        let neg: Vec<f64> = r.iter().map(|x| -x).collect();
        let dirs = sys.directions(&u);
        let at_corner = sys.at_corner(&u);
        let j_copy = at_corner.then(|| j.clone());
        let (j, mut jz, step) = bounded_step(sys, &u, j, &neg, false);
        let mut next = armijo(sys, &u, &r, f, &j, &step, &dirs);
        if let Some(j_copy) = j_copy {
            // At a spline corner: also try preferring the left-hand piece,
            // and keep whichever step reduces the residual more.
            let (j2, jz2, step2) = bounded_step(sys, &u, j_copy, &neg, true);
            if step2 != step {
                if let Some(alt) = armijo(sys, &u, &r, f, &j2, &step2, &dirs) {
                    if next.as_ref().is_none_or(|n| alt.1 < n.1) {
                        next = Some(alt);
                        jz = jz2;
                    }
                }
            }
        }
        if next.is_none() {
            // Levenberg–Marquardt for this iteration.
            let mut m = mu.unwrap_or_else(|| {
                let mut dmax = 0.0f64;
                for k in 0..j.cols() {
                    let c = j.col(k);
                    dmax = dmax.max(dot(c, c) / sys.weights[k]);
                }
                if dmax > 0.0 {
                    1e-3 * dmax
                } else {
                    1e-3
                }
            });
            for _ in 0..64 {
                let h = lm_step(&jz, &r, &sys.weights, m);
                let mut trial: Vec<f64> = u.iter().zip(&h).map(|(a, d)| a + d).collect();
                sys.project(&mut trial);
                if trial == u || !h.iter().all(|x| x.is_finite()) {
                    break;
                }
                let disp: Vec<f64> = trial.iter().zip(&u).map(|(a, b)| a - b).collect();
                let jd = jz.mul_vec(&disp);
                let pred = -dot(&r, &jd) - 0.5 * dot(&jd, &jd);
                let turned = !sys.turn_ok(&dirs, &trial);
                if let Some(rt) = (!turned).then(|| sys.residuals(&trial)).flatten() {
                    let ft = 0.5 * dot(&rt, &rt);
                    if pred > 0.0 {
                        let rho = (f - ft) / pred;
                        if rho > 0.0 {
                            let t = 2.0 * rho - 1.0;
                            m *= (1.0 / 3.0f64).max(1.0 - t * t * t);
                            nu = 2.0;
                            next = Some((trial, ft));
                            break;
                        }
                    }
                }
                m *= nu;
                nu *= 2.0;
                if m.is_nan() || m >= 1e300 {
                    break;
                }
            }
            mu = Some(m.min(1e300));
        }
        let Some((un, fnew)) = next else {
            // No admissible decrease: a stationary point of ½‖r‖². Restart
            // from the foot points if they differ.
            if reprojections < MAX_REPROJECTIONS && sys.reproject(&mut u) {
                reprojections += 1;
                stall = 0;
                mu = None;
                continue;
            }
            return Outcome {
                u,
                converged: false,
                iterations: its,
            };
        };
        if f - fnew <= 1e-12 * f {
            stall += 1;
        } else {
            stall = 0;
        }
        u = un;
        if stall >= 3 {
            // Progress has become negligible: stop at this (near-)stationary
            // point, re-checking convergence there first.
            if let Some((r, _, tol)) = sys.jacobian(&u) {
                if converged(&r, &tol) && sys.angles_ok(&u) {
                    return Outcome {
                        u,
                        converged: true,
                        iterations: its,
                    };
                }
            }
            if reprojections < MAX_REPROJECTIONS && sys.reproject(&mut u) {
                reprojections += 1;
                stall = 0;
                mu = None;
                continue;
            }
            return Outcome {
                u,
                converged: false,
                iterations: its,
            };
        }
    }
}

/// What continuation interpolates (design note 0005 §6).
struct Stage {
    /// Constant slots: (slot, start value, target value).
    consts: Vec<(usize, f64, f64)>,
    /// Residuals at the start (with constants at their start values).
    r0: Vec<f64>,
    /// For angle equations: (start angle, target value).
    angles: Vec<Option<(f64, f64)>>,
}

impl CompSystem<'_> {
    fn set_stage(&mut self, st: &Stage, s: f64) {
        for &(slot, c0, c1) in &st.consts {
            self.v[slot] = if s >= 1.0 { c1 } else { c0 + s * (c1 - c0) };
        }
        for (i, k) in self.kinds.iter_mut().enumerate() {
            match (&st.angles[i], k) {
                (
                    Some((theta0, target)),
                    EqKind::DirPair {
                        mode: PairMode::AngleSin { value },
                        ..
                    },
                ) => {
                    *value = if s >= 1.0 {
                        *target
                    } else {
                        // Shortest way round, in (−π, π].
                        let mut d = target - theta0;
                        if d > std::f64::consts::PI {
                            d -= std::f64::consts::TAU;
                        } else if d <= -std::f64::consts::PI {
                            d += std::f64::consts::TAU;
                        }
                        theta0 + s * d
                    };
                    self.offsets[i] = 0.0;
                }
                _ => {
                    self.offsets[i] = if s >= 1.0 { 0.0 } else { (1.0 - s) * st.r0[i] };
                }
            }
        }
    }
}

/// Result of solving a component.
#[derive(Clone, Debug)]
pub(crate) struct CompResult {
    /// The solution (or the direct attempt's final point on failure).
    pub u: Vec<f64>,
    /// Solved.
    pub converged: bool,
    /// Iterations over all attempts.
    pub iterations: usize,
}

/// Solves a component from `u0`: directly, then by continuation.
///
/// `const_seed` gives, for constant slots referenced by the component, the
/// value they had in the starting geometry (continuation starts there).
pub(crate) fn solve_component(
    sys: &mut CompSystem<'_>,
    u0: Vec<f64>,
    const_seed: &[(usize, f64)],
    max_it: usize,
) -> CompResult {
    let direct = gauss_newton(sys, u0.clone(), max_it);
    let mut iterations = direct.iterations;
    if direct.converged {
        return CompResult {
            u: direct.u,
            converged: true,
            iterations,
        };
    }
    // Again with the turn limit: the direct attempt may have overshot into
    // the basin of another branch.
    sys.turn_limit = true;
    let limited = gauss_newton(sys, u0.clone(), max_it);
    sys.turn_limit = false;
    iterations += limited.iterations;
    if limited.converged {
        return CompResult {
            u: limited.u,
            converged: true,
            iterations,
        };
    }
    // Continuation from the starting geometry.
    let consts: Vec<(usize, f64, f64)> = const_seed
        .iter()
        .map(|&(slot, c0)| (slot, c0, sys.v[slot]))
        .collect();
    let targets: Vec<f64> = consts.iter().map(|c| c.2).collect();
    for &(slot, c0, _) in &consts {
        sys.v[slot] = c0;
    }
    let r0 = sys.residuals(&u0);
    let angles: Vec<Option<(f64, f64)>> = {
        sys.write(&u0);
        sys.kinds
            .iter()
            .map(|k| match k {
                EqKind::DirPair {
                    mode: PairMode::AngleSin { value },
                    ..
                } => k.angle_theta(&sys.v).map(|t| (t, *value)),
                _ => None,
            })
            .collect()
    };
    let restore = |sys: &mut CompSystem<'_>| {
        for (&(slot, _, _), &t) in consts.iter().zip(&targets) {
            sys.v[slot] = t;
        }
    };
    let Some(r0) = r0 else {
        restore(sys);
        sys.set_stage(
            &Stage {
                consts: Vec::new(),
                r0: vec![0.0; sys.m()],
                angles: vec![None; sys.m()],
            },
            1.0,
        );
        return CompResult {
            u: direct.u,
            converged: false,
            iterations,
        };
    };
    let stage = Stage {
        consts: consts.clone(),
        r0,
        angles,
    };
    let mut s = 0.0f64;
    let mut h = 0.5f64;
    let mut u = u0;
    let mut ok = false;
    while h >= 1.0 / 64.0 {
        let s_next = if s + h >= 1.0 { 1.0 } else { s + h };
        sys.set_stage(&stage, s_next);
        let out = gauss_newton(sys, u.clone(), max_it);
        iterations += out.iterations;
        if out.converged {
            s = s_next;
            u = out.u;
            if s >= 1.0 {
                ok = true;
                break;
            }
        } else {
            h *= 0.5;
        }
    }
    // Leave the system at the target (stage 1).
    sys.set_stage(&stage, 1.0);
    if ok {
        CompResult {
            u,
            converged: true,
            iterations,
        }
    } else {
        CompResult {
            u: direct.u,
            converged: false,
            iterations,
        }
    }
}
