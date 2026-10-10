//! The solve pipeline (design note 0005 §5): validate, scale and build,
//! eliminate and decompose, solve each component, write back, check against
//! the model tolerance, diagnose.

use std::collections::BTreeMap;

use glam::DVec2;
use openshape_linalg::{select_columns_in_order, Matrix};

use crate::build::{build, scale_of, Problem, Scale, DEGENERATE};
use crate::diagnose::{conflicts, diagnose, CompDiag};
use crate::eq::SplineDef;
use crate::index::{Circ, EntKind, Index};
use crate::model::{Geometry, Sketch, SketchKey};
use crate::newton::{equation_tol, solve_component, CompSystem, LengthTol, RANK_TOL};
use crate::reduce::{reduce, Reduced, Status};
use crate::spline::combine;
use crate::{EntityState, Solution, SolveError, SolveOptions, SolveStatus};

/// Model linear tolerance (metres), ADR 0006.
const LINEAR_TOL: f64 = 1e-7;
/// Model angular tolerance (radians), ADR 0006.
const ANGULAR_TOL: f64 = 1e-9;
/// Metric weight of dragged points (design note §7).
const DRAG_WEIGHT: f64 = 1e6;
/// Metric weight of hidden curve parameters: they move freely, the points
/// they place carry the cost of moving.
const PARAM_WEIGHT: f64 = 1e-4;

/// A failure before any solving: every entity is reported free and no
/// constraint is counted.
pub(crate) fn early_failure<K: SketchKey>(sketch: &Sketch<K>, err: SolveError) -> Solution<K> {
    let mut entities = BTreeMap::new();
    let mut dof = 0usize;
    for e in &sketch.entities {
        entities.insert(e.id, EntityState::Free);
        dof += match e.geometry {
            Geometry::Point { .. } => 2,
            Geometry::Circle { .. } => 1,
            _ => 0,
        };
    }
    Solution {
        status: SolveStatus::Failed(err),
        sketch: sketch.clone(),
        dof,
        entities,
        redundant: Vec::new(),
        conflicting: Vec::new(),
        max_residual: 0.0,
        iterations: 0,
    }
}

/// Every slot of each entity (points: x, y; circles add the radius).
fn entity_slots<K: SketchKey>(index: &Index<K>, prob: &Problem, ent: usize) -> Vec<usize> {
    let pt = |p: usize| [2 * p, 2 * p + 1];
    match index.kinds[ent] {
        EntKind::Point(p) => pt(p).to_vec(),
        EntKind::Line(l) => [pt(l.s), pt(l.e)].concat(),
        EntKind::Circular(Circ::Circle { c, r }) => {
            let mut v = pt(c).to_vec();
            v.push(prob.radius_slot(r));
            v
        }
        EntKind::Circular(Circ::Arc { c, s, e }) => [pt(c), pt(s), pt(e)].concat(),
        EntKind::Spline(s) => index.splines[s].points.iter().flat_map(|&p| pt(p)).collect(),
    }
}

/// The curve parameter of the point of `def` closest to point `p` (slot
/// values `v`): dense sampling, then safeguarded Newton on
/// `(C(t) − P)·C′(t) = 0` (P&T eq. 6.3).
fn closest_param(def: &SplineDef, v: &[f64], p: usize) -> f64 {
    let pts: Vec<DVec2> = def
        .points
        .iter()
        .map(|&q| DVec2::new(v[2 * q], v[2 * q + 1]))
        .collect();
    let target = DVec2::new(v[2 * p], v[2 * p + 1]);
    let dist2 = |t: f64| {
        let (w, _) = def.map.weights(t);
        (combine(&w, &pts) - target).length_squared()
    };
    let samples = 16 * def.map.count + 16;
    let mut best_t = 0.0;
    let mut best_d = f64::INFINITY;
    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let d = dist2(t);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    for _ in 0..32 {
        let (w, w1, w2) = def.map.weights_ders(best_t, 2);
        let c = combine(&w, &pts) - target;
        let c1 = combine(&w1, &pts);
        let c2 = combine(&w2, &pts);
        let g = c.dot(c1);
        let gp = c1.dot(c1) + c.dot(c2);
        if !(gp > 0.0) {
            break;
        }
        let t = (best_t - g / gp).clamp(0.0, 1.0);
        let d = dist2(t);
        if !(d < best_d) {
            break;
        }
        best_d = d;
        best_t = t;
    }
    best_t
}

/// The message for a degenerate equation.
fn degenerate_message<K: SketchKey>(sketch: &Sketch<K>, prob: &Problem, eq: usize) -> String {
    let e = &prob.eqs[eq];
    match (e.owner, e.arc) {
        (Some(c), _) => format!(
            "constraint {:?} needs a direction that has zero length (a zero-length line, coincident points or a zero-radius arc)",
            sketch.constraints[c].id
        ),
        (None, Some(a)) => format!("arc {:?} has zero radius", sketch.entities[a].id),
        _ => "degenerate geometry".to_string(),
    }
}

/// Evaluates every equation of `prob` at slot values `v`, returning
/// `(residual, tolerance)` per equation (`None` where degenerate).
fn check_all(prob: &Problem, v: &[f64], lin: f64, ang: f64) -> Vec<Option<(f64, f64, bool)>> {
    let mut grad = Vec::new();
    prob.eqs
        .iter()
        .map(|e| {
            grad.clear();
            let ev = e.kind.eval(v, &prob.splines, prob.degen, Some(&mut grad))?;
            let noise: f64 = grad.iter().map(|&(s, g)| g.abs() * v[s].abs()).sum();
            let model = if e.angular { ang } else { lin };
            let tol = model.max(32.0 * f64::EPSILON * noise);
            let branch_ok = e.kind.angle_cos(v, prob.degen).is_none_or(|c| c > 0.0);
            Some((ev.r, tol, branch_ok))
        })
        .collect()
}

/// Absolute slot values of a sketch (positions and radii) with the given
/// hidden parameters.
fn absolute_slots<K: SketchKey>(index: &Index<K>, positions: &[DVec2], radii: &[f64], params: &[f64]) -> Vec<f64> {
    let mut v = Vec::with_capacity(2 * positions.len() + radii.len() + params.len());
    for p in positions {
        v.push(p.x);
        v.push(p.y);
    }
    v.extend_from_slice(radii);
    v.extend_from_slice(params);
    debug_assert_eq!(index.points.len(), positions.len());
    v
}

/// Largest absolute residual (metres or radians) over the evaluable
/// equations.
fn max_abs_residual(checks: &[Option<(f64, f64, bool)>]) -> f64 {
    checks
        .iter()
        .flatten()
        .map(|(r, _, _)| r.abs())
        .filter(|r| r.is_finite())
        .fold(0.0, f64::max)
}

pub(crate) fn solve<K: SketchKey>(sketch: &Sketch<K>, options: &SolveOptions<K>) -> Solution<K> {
    let index = match Index::new(sketch) {
        Ok(i) => i,
        Err(e) => return early_failure(sketch, e),
    };
    let cons = match index.constraints(sketch) {
        Ok(c) => c,
        Err(e) => return early_failure(sketch, e),
    };
    let drag = match index.drag(&options.drag) {
        Ok(d) => d,
        Err(e) => return early_failure(sketch, e),
    };
    let scale = scale_of(&index, &cons, &drag);
    let prob = build(&index, &cons, scale);
    let lt = (0.5 * LINEAR_TOL / scale.l).min(1e-12);
    let len_tol = LengthTol {
        base: lt,
        cap: (0.5 * LINEAR_TOL / scale.l).min(1e-10),
    };
    let red = reduce(&prob, len_tol.base);
    let n_slots = prob.n_slots();
    let np = prob.n_points;

    // Starting values: absolute, scaled, and metric weights per slot.
    let mut abs_seed = vec![0.0; n_slots];
    let mut slot_w = vec![1.0; n_slots];
    let mut dragged = vec![false; n_slots];
    for (p, at) in index.positions.iter().enumerate() {
        abs_seed[2 * p] = at.x;
        abs_seed[2 * p + 1] = at.y;
    }
    for &(p, at) in &drag {
        abs_seed[2 * p] = at.x;
        abs_seed[2 * p + 1] = at.y;
        for s in [2 * p, 2 * p + 1] {
            dragged[s] = true;
            slot_w[s] = DRAG_WEIGHT;
        }
    }
    for (r, rad) in index.radii.iter().enumerate() {
        abs_seed[prob.radius_slot(r)] = *rad;
    }
    for k in 0..prob.params.len() {
        slot_w[prob.param_slot(k)] = PARAM_WEIGHT;
    }
    let to_scaled = |s: usize, a: f64| -> f64 {
        if s < 2 * np {
            if s % 2 == 0 {
                (a - scale.cx) / scale.l
            } else {
                (a - scale.cy) / scale.l
            }
        } else {
            a / scale.l
        }
    };
    let seed: Vec<f64> = (0..n_slots).map(|s| to_scaled(s, abs_seed[s])).collect();

    // Class starting values: dragged members win; otherwise the mean (exact
    // when all members agree).
    let class_seed: Vec<f64> = red
        .classes
        .iter()
        .map(|cl| {
            let dm: Vec<usize> = cl.members.iter().copied().filter(|&s| dragged[s]).collect();
            let set = if dm.is_empty() { cl.members.clone() } else { dm };
            if set.iter().all(|&s| seed[s].to_bits() == seed[set[0]].to_bits()) {
                seed[set[0]]
            } else {
                set.iter().map(|&s| seed[s]).sum::<f64>() / set.len() as f64
            }
        })
        .collect();
    let mut v = vec![0.0; n_slots];
    for (c, cl) in red.classes.iter().enumerate() {
        let val = cl.constant.map_or(class_seed[c], |k| k.0);
        for &s in &cl.members {
            v[s] = val;
        }
    }
    // Hidden curve parameters start at the closest point.
    for (k, par) in prob.params.iter().enumerate() {
        let t = closest_param(&prob.splines[par.spline], &v, par.point);
        v[prob.param_slot(k)] = t;
    }
    let mut u_global: Vec<f64> = red
        .unknowns
        .iter()
        .map(|&c| v[red.classes[c].members[0]])
        .collect();
    let u_seed = u_global.clone();
    let unknown_w: Vec<f64> = red
        .unknowns
        .iter()
        .map(|&c| red.classes[c].members.iter().map(|&s| slot_w[s]).sum())
        .collect();

    // Degenerate geometry at the start.
    for (i, e) in prob.eqs.iter().enumerate() {
        if matches!(red.status[i], Status::General | Status::Demoted)
            && e.kind.eval(&v, &prob.splines, prob.degen, None).is_none()
        {
            let msg = degenerate_message(sketch, &prob, i);
            return early_failure(sketch, SolveError::Degenerate(msg));
        }
    }

    // Solve every component.
    let mut iterations = 0usize;
    let mut failed = false;
    let mut conflicting_eqs: Vec<usize> = Vec::new();
    let mut redundant_eqs: Vec<usize> = Vec::new();
    let mut slot_free = vec![false; n_slots];
    let mut dof = 0usize;
    let mut comp_rank: Vec<usize> = Vec::with_capacity(red.components.len());
    for comp in &red.components {
        let mut sys = CompSystem::new(&prob, &red, comp, v.clone(), &unknown_w, len_tol);
        let u0: Vec<f64> = comp.unknowns.iter().map(|&u| u_global[u]).collect();
        // Constant slots referenced, with their values in the start geometry.
        let mut const_seed: Vec<(usize, f64)> = Vec::new();
        for &e in &comp.eqs {
            for s in prob.eqs[e].kind.slots(&prob.splines) {
                let c = red.class_of[s];
                if red.classes[c].constant.is_some() {
                    const_seed.push((s, class_seed[c]));
                }
            }
        }
        const_seed.sort_by_key(|x| x.0);
        const_seed.dedup_by_key(|x| x.0);
        let res = solve_component(&mut sys, u0, &const_seed, options.max_iterations);
        iterations += res.iterations;
        for (k, &u) in comp.unknowns.iter().enumerate() {
            u_global[u] = res.u[k];
        }
        let hidden: Vec<bool> = comp
            .unknowns
            .iter()
            .map(|&u| prob.is_param(red.classes[red.unknowns[u]].members[0]))
            .collect();
        let d: CompDiag = diagnose(&mut sys, &res.u, &hidden);
        if !res.converged {
            failed = true;
            for i in conflicts(&mut sys, &res.u) {
                conflicting_eqs.push(comp.eqs[i]);
            }
        }
        dof += d.dof;
        comp_rank.push(d.rank);
        for (k, &u) in comp.unknowns.iter().enumerate() {
            if d.free[k] {
                for &s in &red.classes[red.unknowns[u]].members {
                    slot_free[s] = true;
                }
            }
        }
        for i in d.redundant {
            redundant_eqs.push(comp.eqs[i]);
        }
    }
    for &u in &red.free_unknowns {
        dof += 1;
        for &s in &red.classes[red.unknowns[u]].members {
            slot_free[s] = true;
        }
    }
    // Equations whose unknowns were all eliminated: satisfied (hence implied
    // by the eliminations) or in conflict with them.
    {
        let mut grad = Vec::new();
        for &i in &red.constant_eqs {
            let e = &prob.eqs[i];
            grad.clear();
            match e.kind.eval(&v, &prob.splines, prob.degen, Some(&mut grad)) {
                Some(ev) => {
                    let noise: f64 = grad.iter().map(|&(s, g)| g.abs() * v[s].abs()).sum();
                    let tol = equation_tol(e.angular, ev.lever, noise, len_tol);
                    let branch_ok = e.kind.angle_cos(&v, prob.degen).is_none_or(|c| c > 0.0);
                    if ev.r.abs() <= tol && branch_ok {
                        redundant_eqs.push(i);
                    } else {
                        failed = true;
                        conflicting_eqs.push(i);
                    }
                }
                None => {
                    failed = true;
                    conflicting_eqs.push(i);
                }
            }
        }
    }
    for (i, st) in red.status.iter().enumerate() {
        if *st == Status::Redundant {
            redundant_eqs.push(i);
        }
    }

    // Final slot values (scaled), for diagnostics of the full system.
    let mut v_final = v.clone();
    for (u, &c) in red.unknowns.iter().enumerate() {
        for &s in &red.classes[c].members {
            v_final[s] = u_global[u];
        }
    }
    let redundant_eqs = priority_redundancy(&prob, &red, &v_final, redundant_eqs, &comp_rank);

    // Entity states.
    let mut entities = BTreeMap::new();
    for (i, e) in sketch.entities.iter().enumerate() {
        let free = entity_slots(&index, &prob, i).iter().any(|&s| slot_free[s]);
        entities.insert(
            e.id,
            if free {
                EntityState::Free
            } else {
                EntityState::Constrained
            },
        );
    }
    let owners = |eqs: &[usize]| -> Vec<K> {
        let mut pos: Vec<usize> = eqs.iter().filter_map(|&i| prob.eqs[i].owner).collect();
        pos.sort_unstable();
        pos.dedup();
        pos.into_iter().map(|c| sketch.constraints[c].id).collect()
    };
    let redundant = owners(&redundant_eqs);
    let conflicting = owners(&conflicting_eqs);

    // The absolute system, for the model-tolerance check and max_residual.
    let mut abs_prob = build(&index, &cons, Scale::identity());
    abs_prob.degen = DEGENERATE * scale.l;
    let params: Vec<f64> = (0..prob.params.len())
        .map(|k| v_final[prob.param_slot(k)])
        .collect();
    let input_params: Vec<f64> = (0..prob.params.len())
        .map(|k| v[prob.param_slot(k)])
        .collect();
    let input_residual = || {
        let va = absolute_slots(&index, &index.positions, &index.radii, &input_params);
        max_abs_residual(&check_all(&abs_prob, &va, LINEAR_TOL, ANGULAR_TOL))
    };

    if failed {
        let err = if conflicting.is_empty() {
            SolveError::NotConverged
        } else {
            SolveError::Conflicting
        };
        return Solution {
            status: SolveStatus::Failed(err),
            sketch: sketch.clone(),
            dof,
            entities,
            redundant,
            conflicting,
            max_residual: input_residual(),
            iterations,
        };
    }

    // Write back.
    let abs_value = |s: usize| -> f64 {
        let c = red.class_of[s];
        let cl = &red.classes[c];
        if let Some((_, a)) = cl.constant {
            return a;
        }
        let u = red.unknown_of_class[c].expect("non-constant class has an unknown");
        let unchanged = u_global[u].to_bits() == u_seed[u].to_bits()
            && cl
                .members
                .iter()
                .all(|&m| abs_seed[m].to_bits() == abs_seed[cl.members[0]].to_bits());
        if unchanged {
            abs_seed[s]
        } else if s < 2 * np {
            if s % 2 == 0 {
                scale.cx + scale.l * u_global[u]
            } else {
                scale.cy + scale.l * u_global[u]
            }
        } else {
            scale.l * u_global[u]
        }
    };
    let positions: Vec<DVec2> = (0..np)
        .map(|p| DVec2::new(abs_value(2 * p), abs_value(2 * p + 1)))
        .collect();
    let radii: Vec<f64> = (0..prob.n_circles)
        .map(|r| abs_value(prob.radius_slot(r)))
        .collect();
    let va = absolute_slots(&index, &positions, &radii, &params);
    let checks = check_all(&abs_prob, &va, LINEAR_TOL, ANGULAR_TOL);
    let all_ok = checks
        .iter()
        .all(|c| matches!(c, Some((r, tol, true)) if r.abs() <= *tol))
        && positions.iter().all(|p| p.is_finite())
        && radii.iter().all(|r| r.is_finite() && *r > 0.0);
    if !all_ok {
        return Solution {
            status: SolveStatus::Failed(SolveError::NotConverged),
            sketch: sketch.clone(),
            dof,
            entities,
            redundant,
            conflicting,
            max_residual: input_residual(),
            iterations,
        };
    }
    let mut out = sketch.clone();
    for (p, &ent) in index.points.iter().enumerate() {
        if let Geometry::Point { at } = &mut out.entities[ent].geometry {
            *at = positions[p];
        }
    }
    for (r, &ent) in index.circles.iter().enumerate() {
        if let Geometry::Circle { radius, .. } = &mut out.entities[ent].geometry {
            *radius = radii[r];
        }
    }
    let status = if redundant.is_empty() {
        SolveStatus::Solved
    } else {
        SolveStatus::Redundant
    };
    Solution {
        status,
        sketch: out,
        dof,
        entities,
        redundant,
        conflicting: Vec::new(),
        max_residual: max_abs_residual(&checks),
        iterations,
    }
}

/// Makes the redundancy report follow priority order across eliminations.
///
/// The reduced analysis applies every elimination first. That matches the
/// priority order unless an equation it found redundant precedes, in the
/// same connected set of slots, an elimination that it depends on; then the
/// later elimination is the one implied. For such sets the independent
/// equations are chosen afresh, in priority order, from the full Jacobian
/// (rank-capped so the count matches the reduced analysis).
fn priority_redundancy(
    prob: &Problem,
    red: &Reduced,
    v: &[f64],
    redundant: Vec<usize>,
    comp_rank: &[usize],
) -> Vec<usize> {
    let n = prob.n_slots();
    let neq = prob.eqs.len();
    let slots: Vec<Vec<usize>> = prob.eqs.iter().map(|e| e.kind.slots(&prob.splines)).collect();
    let mut uf: Vec<usize> = (0..n).collect();
    fn find(uf: &mut [usize], mut x: usize) -> usize {
        while uf[x] != x {
            uf[x] = uf[uf[x]];
            x = uf[x];
        }
        x
    }
    for ss in &slots {
        for w in ss.windows(2) {
            let (a, b) = (find(&mut uf, w[0]), find(&mut uf, w[1]));
            if a != b {
                let (lo, hi) = (a.min(b), a.max(b));
                uf[hi] = lo;
            }
        }
    }
    let root_of_eq: Vec<Option<usize>> = slots
        .iter()
        .map(|ss| ss.first().map(|&s| find(&mut uf, s)))
        .collect();
    let mut is_red = vec![false; neq];
    for &i in &redundant {
        is_red[i] = true;
    }
    // Per full set: earliest redundant numerical equation, latest elimination.
    let mut first_red: BTreeMap<usize, usize> = BTreeMap::new();
    let mut last_sub: BTreeMap<usize, usize> = BTreeMap::new();
    for i in 0..neq {
        let Some(root) = root_of_eq[i] else { continue };
        match red.status[i] {
            Status::Accepted | Status::Redundant => {
                last_sub.insert(root, i);
            }
            Status::General | Status::Demoted => {
                if is_red[i] {
                    first_red.entry(root).or_insert(i);
                }
            }
        }
    }
    let redo: Vec<usize> = first_red
        .iter()
        .filter(|(root, fr)| last_sub.get(root).is_some_and(|ls| ls > fr))
        .map(|(root, _)| *root)
        .collect();
    if redo.is_empty() {
        let mut r = redundant;
        r.sort_unstable();
        r.dedup();
        return r;
    }
    let mut out: Vec<usize> = redundant
        .into_iter()
        .filter(|&i| root_of_eq[i].is_none_or(|r| !redo.contains(&r)))
        .collect();
    for root in redo {
        let eqs: Vec<usize> = (0..neq).filter(|&i| root_of_eq[i] == Some(root)).collect();
        let set_slots: Vec<usize> = (0..n).filter(|&s| find(&mut uf, s) == root).collect();
        let mut col = vec![usize::MAX; n];
        for (k, &s) in set_slots.iter().enumerate() {
            col[s] = k;
        }
        // Rank of the full set: independent eliminations plus the ranks of
        // the reduced components inside it.
        let mut rank = eqs
            .iter()
            .filter(|&&i| red.status[i] == Status::Accepted)
            .count();
        for (ci, comp) in red.components.iter().enumerate() {
            if comp.eqs.first().is_some_and(|&e| root_of_eq[e] == Some(root)) {
                rank += comp_rank[ci];
            }
        }
        // Jᵀ: one column per equation, in priority order.
        let mut jt = Matrix::zeros(set_slots.len(), eqs.len());
        let mut grad = Vec::new();
        for (c, &i) in eqs.iter().enumerate() {
            grad.clear();
            if prob.eqs[i]
                .kind
                .eval(v, &prob.splines, prob.degen, Some(&mut grad))
                .is_some()
            {
                for &(s, g) in &grad {
                    jt[(col[s], c)] += g;
                }
            }
        }
        let sel = select_columns_in_order(&jt, RANK_TOL, Some(rank));
        out.extend(sel.rejected.iter().map(|&c| eqs[c]));
    }
    out.sort_unstable();
    out.dedup();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closest_param_on_a_segment() {
        let def = SplineDef {
            map: crate::spline::SplineMap::control(2, 1).unwrap(),
            points: vec![0, 1],
        };
        // P0 = (0,0), P1 = (1,0), query (0.3, 0.5).
        let v = vec![0.0, 0.0, 1.0, 0.0, 0.3, 0.5];
        let t = closest_param(&def, &v, 2);
        assert!((t - 0.3).abs() < 1e-12);
    }
}
