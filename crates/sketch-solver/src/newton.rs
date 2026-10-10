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
//! * **Armijo backtracking** on `½‖r‖²` (c = 1e-4, halving down to 2⁻²⁰;
//!   Nocedal & Wright §3.1).
//! * **Levenberg–Marquardt fallback** when the line search fails: steps from
//!   `(JᵀJ + μW) h = −Jᵀr`, solved as an augmented least-squares problem,
//!   with the gain-ratio update of μ of Madsen, Nielsen & Tingleff, *Methods
//!   for Non-Linear Least Squares Problems*, DTU 2004, §3.2 (eq. 3.14–3.16).
//! * **Convergence** per equation: `|rᵢ| ≤ tolᵢ`, with `tolᵢ = 1e-12`
//!   (scaled) for length residuals and, for angle residuals, the angle whose
//!   lever (the shorter direction involved) moves a point by 1e-12 of the
//!   sketch size, clamped to [1e-12, 1e-10]; never below the rounding noise
//!   of the residual at the current point, never above 1e-10.
//! * **Continuation** (§6): if the direct solve fails (or lands on the
//!   antipodal branch of an angle), the change from the start is applied in
//!   equal steps of 1/2, halving the step on failure down to 1/64, each
//!   solved from the previous result. Dimensions and fixed values are
//!   interpolated from their current values; other equations follow the
//!   Newton homotopy `r(x) − (1 − s)·r(x₀)`.

use openshape_linalg::{dot, lstsq_min_norm, weighted_min_norm, Matrix};

use crate::build::Problem;
use crate::eq::{EqKind, PairMode};
use crate::reduce::{Component, Reduced};

/// Relative rank threshold for steps and diagnostics (design note §5).
pub(crate) const RANK_TOL: f64 = 1e-10;
/// Upper bound on any equation's convergence tolerance (scaled).
const TOL_CAP: f64 = 1e-10;
/// Base tolerance for angle residuals.
const TOL_ANGLE: f64 = 1e-12;

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
    ) -> Self {
        let mut col_of_slot = vec![usize::MAX; prob.n_slots()];
        let mut members = Vec::with_capacity(comp.unknowns.len());
        let mut is_param = Vec::with_capacity(comp.unknowns.len());
        let mut is_radius = Vec::with_capacity(comp.unknowns.len());
        let mut weights = Vec::with_capacity(comp.unknowns.len());
        for (k, &u) in comp.unknowns.iter().enumerate() {
            let class = &red.classes[red.unknowns[u]];
            for &s in &class.members {
                col_of_slot[s] = k;
            }
            members.push(class.members.clone());
            is_param.push(prob.is_param(class.members[0]));
            is_radius.push(prob.is_radius(class.members[0]));
            weights.push(unknown_weight[u]);
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
            .all(|(x, &rad)| x.is_finite() && (!rad || *x > self.prob.degen))
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
            let e = self.kinds[i].eval(&self.v, &self.prob.splines, self.prob.degen, Some(&mut grad))?;
            let ri = e.r - self.offsets[i];
            if !ri.is_finite() {
                self.grad = grad;
                return None;
            }
            let mut noise = 0.0;
            for &(s, g) in &grad {
                noise += g.abs() * self.v[s].abs();
                let c = self.col_of_slot[s];
                if c != usize::MAX {
                    j[(i, c)] += g;
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

    /// True if every angle equation is on its intended branch:
    /// `cos(θ − v) > 0`.
    pub fn angles_ok(&mut self, u: &[f64]) -> bool {
        self.write(u);
        self.kinds
            .iter()
            .all(|k| k.angle_cos(&self.v, self.prob.degen).is_none_or(|c| c > 0.0))
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

/// The Gauss–Newton step with hidden curve parameters held at a bound of
/// `[0, 1]` when the step would push them out (the active set of projected
/// Newton methods: Bertsekas, "Projected Newton methods for optimization
/// problems with simple constraints", SIAM J. Control Optim. 20(2), 1982).
/// Without this, a point sliding off the end of a spline would be stuck:
/// every step would put the along-curve correction into the cheap
/// parameter, and the projection would throw it away.
///
/// Returns the Jacobian with the held columns zeroed (for the line search
/// and the Levenberg–Marquardt fallback) and the step (zero in those
/// columns).
fn bounded_step(sys: &CompSystem<'_>, u: &[f64], mut j: Matrix, neg: &[f64]) -> (Matrix, Vec<f64>) {
    let n = j.cols();
    let mut held = vec![false; n];
    loop {
        let mut step = weighted_min_norm(&j, neg, &sys.weights, RANK_TOL).x;
        for (x, &h) in step.iter_mut().zip(&held) {
            if h {
                *x = 0.0;
            }
        }
        let mut changed = false;
        for k in 0..n {
            if sys.is_param[k]
                && !held[k]
                && ((u[k] <= 0.0 && step[k] < 0.0) || (u[k] >= 1.0 && step[k] > 0.0))
            {
                held[k] = true;
                changed = true;
                j.col_mut(k).fill(0.0);
            }
        }
        if !changed {
            return (j, step);
        }
    }
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

/// Damped Gauss–Newton from `u` (see the module docs).
pub(crate) fn gauss_newton(sys: &mut CompSystem<'_>, mut u: Vec<f64>, max_it: usize) -> Outcome {
    let mut its = 0usize;
    let mut mu: Option<f64> = None;
    let mut nu = 2.0f64;
    let mut stall = 0usize;
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
        let (j, step) = bounded_step(sys, &u, j, &neg);
        let mut next: Option<(Vec<f64>, f64)> = None;
        if step.iter().all(|s| s.is_finite()) {
            let slope = dot(&r, &j.mul_vec(&step));
            if slope < 0.0 {
                let mut alpha = 1.0f64;
                for _ in 0..=20 {
                    let mut trial: Vec<f64> =
                        u.iter().zip(&step).map(|(a, d)| a + alpha * d).collect();
                    let clamped = sys.project(&mut trial);
                    if trial == u {
                        break;
                    }
                    if let Some(rt) = sys.residuals(&trial) {
                        let ft = 0.5 * dot(&rt, &rt);
                        let pred = if clamped {
                            let disp: Vec<f64> = trial.iter().zip(&u).map(|(a, b)| a - b).collect();
                            dot(&r, &j.mul_vec(&disp))
                        } else {
                            alpha * slope
                        };
                        if pred < 0.0 && ft <= f + 1e-4 * pred {
                            next = Some((trial, ft));
                            break;
                        }
                    }
                    alpha *= 0.5;
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
                let h = lm_step(&j, &r, &sys.weights, m);
                let mut trial: Vec<f64> = u.iter().zip(&h).map(|(a, d)| a + d).collect();
                sys.project(&mut trial);
                if trial == u || !h.iter().all(|x| x.is_finite()) {
                    break;
                }
                let disp: Vec<f64> = trial.iter().zip(&u).map(|(a, b)| a - b).collect();
                let jd = j.mul_vec(&disp);
                let pred = -dot(&r, &jd) - 0.5 * dot(&jd, &jd);
                if let Some(rt) = sys.residuals(&trial) {
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
                if !(m < 1e300) {
                    break;
                }
            }
            mu = Some(m.min(1e300));
        }
        let Some((un, fnew)) = next else {
            // No admissible decrease: a stationary point of ½‖r‖².
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
