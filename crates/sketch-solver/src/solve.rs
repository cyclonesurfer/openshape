//! The solve pipeline (design note 0005 §5): validate, scale and build,
//! eliminate and decompose, solve each component, diagnose, write back and
//! check the result against the model tolerance ([`accept`]: also
//! non-degenerate at the result's own scale, so that it can be solved
//! again). If solving fails but the undragged input already passes that
//! check, the input is the answer. A result with a collapsed curve
//! ([`collapsed`]) is no answer: the constraint to blame is searched for
//! ([`Part::attribute`], design note 0005 §9) and reported as
//! conflicting, or the curve is reported degenerate. Nor is one that leaves
//! an offset invalid ([`offset_problem`]): the offset is reported.

use std::collections::BTreeMap;

use glam::DVec2;
use openshape_linalg::{dot, norm2, residual, select_columns_in_order, CompleteOrthogonal, Matrix};

use crate::build::{build, scale_of, ParamCurve, Problem, Scale, DEGENERATE};
use crate::diagnose::{conflicts, diagnose, CompDiag, FREE_TOL};
use crate::ellipse;
use crate::eq::EqKind;
use crate::index::{Circ, Con, Ell, EntKind, Index, OCurve, OEq, MAX_MAGNITUDE};
use crate::model::{Geometry, Sketch, SketchKey};
use crate::newton::{
    equation_tol, solve_component, CompSystem, LengthTol, Region, RANK_TOL, REACH,
};
use crate::offset::{self, Invalid, OffsetPiece, Piece};
use crate::reduce::{peel, reduce, Reduced, Status};
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
        EntKind::Ellipse(_) => entity_points(index, ent).into_iter().flat_map(pt).collect(),
        EntKind::Spline(s) => index.splines[s]
            .points
            .iter()
            .flat_map(|&p| pt(p))
            .collect(),
    }
}

/// The message for a degenerate equation.
fn degenerate_message<K: SketchKey>(sketch: &Sketch<K>, prob: &Problem, eq: usize) -> String {
    let e = &prob.eqs[eq];
    match (e.owner, e.arc) {
        (Some(c), _) => format!(
            "constraint {:?} needs a direction that has zero length (a zero-length line, coincident points or a zero-radius arc)",
            sketch.constraints[c].id
        ),
        (None, Some(a)) => match sketch.entities[a].geometry {
            Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } => {
                format!("ellipse {:?} has an axis of zero length", sketch.entities[a].id)
            }
            _ => format!("arc {:?} has zero radius", sketch.entities[a].id),
        },
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
fn absolute_slots<K: SketchKey>(
    index: &Index<K>,
    positions: &[DVec2],
    radii: &[f64],
    params: &[f64],
) -> Vec<f64> {
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

/// A block of columns of the whole system's null space: a component's
/// orthonormal basis, or a free unknown's own direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum NullBlock {
    Component(usize),
    Free(usize),
}

pub(crate) fn solve<K: SketchKey>(sketch: &Sketch<K>, options: &SolveOptions<K>) -> Solution<K> {
    solve_inner(sketch, options, true)
}

/// The pipeline; `attribute` enables the search for the constraint that
/// collapses a curve ([`Part::attribute`], which solves without it).
fn solve_inner<K: SketchKey>(
    sketch: &Sketch<K>,
    options: &SolveOptions<K>,
    attribute: bool,
) -> Solution<K> {
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
    let mut red = reduce(&prob, len_tol.base);
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
            let set = if dm.is_empty() {
                cl.members.clone()
            } else {
                dm
            };
            if set
                .iter()
                .all(|&s| seed[s].to_bits() == seed[set[0]].to_bits())
            {
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
        let t = match par.curve {
            ParamCurve::Spline(s) => prob.splines[s].closest_param(&v, par.point),
            ParamCurve::Ellipse(el) => {
                let at = |q: usize| DVec2::new(v[2 * q], v[2 * q + 1]);
                let c = at(el.c);
                ellipse::closest_param(c, at(el.m) - c, at(el.n) - c, at(par.point))
            }
        };
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
    // Images of derived constraints move first (`newton::images_first`):
    // an unknown is an image when every slot it stands for is, and none is
    // dragged (a dragged image moves last, like any dragged point).
    let mut image_slot = vec![false; n_slots];
    for con in &cons {
        match con {
            Con::Mirror(_, d) => {
                for &(_, q) in &d.points {
                    image_slot[2 * q] = true;
                    image_slot[2 * q + 1] = true;
                }
            }
            Con::LinearPattern(_, pt) | Con::CircularPattern(_, _, pt) => {
                for &(_, q, _) in &pt.points {
                    image_slot[2 * q] = true;
                    image_slot[2 * q + 1] = true;
                }
            }
            Con::Offset(o) => {
                for &(x, ..) in &o.points {
                    image_slot[2 * x] = true;
                    image_slot[2 * x + 1] = true;
                }
                for &(_, img, _) in &o.pieces {
                    if let OCurve::Circle { r, .. } = img {
                        image_slot[prob.radius_slot(r)] = true;
                    }
                }
            }
            _ => {}
        }
    }
    let unknown_image: Vec<bool> = red
        .unknowns
        .iter()
        .map(|&c| {
            red.classes[c]
                .members
                .iter()
                .all(|&s| image_slot[s] && !dragged[s])
        })
        .collect();
    // Images that only their defining linear equation references are
    // computed from it after the solve instead (`reduce::peel`).
    peel(&prob, &mut red, &unknown_image);

    // Degenerate geometry at the start: arcs without a radius, and
    // directions that constraints need.
    let span = |a: usize, b: usize| (v[2 * b] - v[2 * a]).hypot(v[2 * b + 1] - v[2 * a + 1]);
    let gone = |x: f64| x.is_nan() || x <= prob.degen;
    for (i, k) in index.kinds.iter().enumerate() {
        let msg = match *k {
            EntKind::Circular(Circ::Arc { c, s, .. }) if gone(span(c, s)) => {
                format!("arc {:?} has zero radius", sketch.entities[i].id)
            }
            EntKind::Ellipse(el) if gone(span(el.c, el.m)) || gone(span(el.c, el.n)) => {
                format!(
                    "ellipse {:?} has an axis of zero length",
                    sketch.entities[i].id
                )
            }
            _ => continue,
        };
        return early_failure(sketch, SolveError::Degenerate(msg));
    }
    for (i, e) in prob.eqs.iter().enumerate() {
        if matches!(red.status[i], Status::General | Status::Demoted)
            && e.kind.eval(&v, &prob.splines, prob.degen, None).is_none()
        {
            let msg = degenerate_message(sketch, &prob, i);
            return early_failure(sketch, SolveError::Degenerate(msg));
        }
    }

    // Where geometry may go (scaled): within REACH sketch extents of the
    // centre, and within the coordinate limit of input validation.
    let bounds = |c: f64| {
        (
            ((-MAX_MAGNITUDE - c) / scale.l).max(-REACH),
            ((MAX_MAGNITUDE - c) / scale.l).min(REACH),
        )
    };
    let region = Region {
        x: bounds(scale.cx),
        y: bounds(scale.cy),
        r_max: (MAX_MAGNITUDE / scale.l).min(REACH),
    };

    // Solve every component.
    let mut iterations = 0usize;
    let mut failed = false;
    let mut conflicting_eqs: Vec<usize> = Vec::new();
    let mut redundant_eqs: Vec<usize> = Vec::new();
    let mut slot_free = vec![false; n_slots];
    let mut dof = 0usize;
    let mut comp_rank: Vec<usize> = Vec::with_capacity(red.components.len());
    // Each unknown's row of the null space of the whole system, as blocks:
    // (component, its row there) or a free unknown's own column.
    let mut null_rows: Vec<Vec<(NullBlock, Vec<f64>)>> = vec![Vec::new(); red.unknowns.len()];
    for (ci, comp) in red.components.iter().enumerate() {
        let mut sys = CompSystem::new(
            &prob,
            &red,
            comp,
            v.clone(),
            &unknown_w,
            &unknown_image,
            len_tol,
            region,
        );
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
        if !red.peeled.is_empty() {
            for (k, &u) in comp.unknowns.iter().enumerate() {
                null_rows[u] = vec![(NullBlock::Component(ci), d.null.row(k))];
            }
        }
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
        null_rows[u] = vec![(NullBlock::Free(u), vec![1.0])];
    }
    // Peeled images, from the unknowns they depend on: their values, and
    // their null-space rows (the same combination of the others' rows).
    for pe in &red.peeled {
        let EqKind::Lin { terms, rhs } = &prob.eqs[pe.eq].kind else {
            unreachable!("only linear equations are peeled")
        };
        let mut value = *rhs;
        let mut row: BTreeMap<NullBlock, Vec<f64>> = BTreeMap::new();
        for &(s, c) in terms {
            let class = red.class_of[s];
            match red.unknown_of_class[class] {
                Some(w) if w == pe.unknown => {}
                Some(w) => {
                    value -= c * u_global[w];
                    for (block, r) in &null_rows[w] {
                        let acc = row.entry(*block).or_insert_with(|| vec![0.0; r.len()]);
                        for (a, x) in acc.iter_mut().zip(r) {
                            *a -= c / pe.coef * x;
                        }
                    }
                }
                None => value -= c * red.classes[class].constant.map_or(0.0, |k| k.0),
            }
        }
        u_global[pe.unknown] = value / pe.coef;
        let norm: f64 = row
            .values()
            .flat_map(|r| r.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        if norm > FREE_TOL {
            for &s in &red.classes[red.unknowns[pe.unknown]].members {
                slot_free[s] = true;
            }
        }
        null_rows[pe.unknown] = row.into_iter().collect();
    }
    // Equations whose unknowns were all eliminated: satisfied (hence implied
    // by the eliminations) or in conflict with them.
    let mut const_conflicts: Vec<usize> = Vec::new();
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
                        const_conflicts.push(i);
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
    let conflicting = owners(&conflicting_eqs);
    // A constraint with several equations can be partly implied and partly
    // in conflict (a second fix of a point at another x but the same y):
    // it is reported as conflicting only.
    let redundant: Vec<K> = owners(&redundant_eqs)
        .into_iter()
        .filter(|k| !conflicting.contains(k))
        .collect();

    // The absolute system, for the model-tolerance check and max_residual.
    let mut abs_prob = build(&index, &cons, Scale::identity());
    let params: Vec<f64> = (0..prob.params.len())
        .map(|k| v_final[prob.param_slot(k)])
        .collect();
    let input_params: Vec<f64> = (0..prob.params.len())
        .map(|k| v[prob.param_slot(k)])
        .collect();
    let degen_in = DEGENERATE * scale.l * ACCEPT_MARGIN;
    let (input_ok, input_residual) = accept(
        &mut abs_prob,
        &index,
        &index.positions,
        &index.radii,
        &input_params,
        degen_in,
    );
    let input_ok =
        input_ok && offset_problem(sketch, &index, &cons, &index.positions, &index.radii).is_none();

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
    // The result must also be non-degenerate at its own scale, which is
    // the scale it will be solved at next.
    let out_scale = {
        let mut out_index = index.clone();
        out_index.positions.clone_from(&positions);
        out_index.radii.clone_from(&radii);
        scale_of(&out_index, &cons, &[])
    };
    let degen_out = DEGENERATE * scale.l.max(out_scale.l) * ACCEPT_MARGIN;
    // Constraints that collapse a curve (a line to a point, a circle or arc
    // to its centre, an arc onto itself) have no valid solution, however
    // well the equations hold: blame the constraint whose removal avoids
    // the collapse, like a conflict, or report the curve degenerate. The
    // input fallback does not apply to a converged result: an input this
    // close to the collapse is degenerate under its constraints too. A
    // solve that failed with no conflict found while heading into a
    // collapse (a radius that cannot reach zero) is blamed the same way
    // when the input is no solution, and otherwise left as it is.
    let fallback = input_ok && drag.is_empty();
    let collapse = if (!failed || conflicting.is_empty()) && !(failed && fallback) {
        collapsed(
            &index,
            &positions,
            &radii,
            Some((&index.positions, degen_in)),
            degen_out,
        )
    } else {
        Vec::new()
    };
    if let Some(&first) = collapse.first() {
        let blamed = if attribute {
            let reached = Reached {
                prob: &prob,
                v: &v_final,
                degen: degen_out,
                degen_in,
            };
            let (blamed, its) = blame_collapse(sketch, options, &index, &reached, &collapse);
            iterations += its;
            blamed
        } else {
            Err(first)
        };
        match blamed {
            Ok(culprits) => {
                let conflicting: Vec<K> =
                    culprits.iter().map(|&j| sketch.constraints[j].id).collect();
                return Solution {
                    status: SolveStatus::Failed(SolveError::Conflicting),
                    sketch: sketch.clone(),
                    dof,
                    entities,
                    redundant: redundant
                        .into_iter()
                        .filter(|k| !conflicting.contains(k))
                        .collect(),
                    conflicting,
                    max_residual: input_residual,
                    iterations,
                };
            }
            Err((ent, how)) if !failed => {
                let kind = match index.kinds[ent] {
                    EntKind::Line(_) => "line",
                    EntKind::Circular(Circ::Circle { .. }) => "circle",
                    EntKind::Ellipse(Ell { arc: None, .. }) => "ellipse",
                    EntKind::Ellipse(_) => "elliptical arc",
                    _ => "arc",
                };
                let msg = format!("{kind} {:?} {}", sketch.entities[ent].id, how.says());
                return early_failure(sketch, SolveError::Degenerate(msg));
            }
            Err(_) => {}
        }
    }
    let (out_ok, max_residual) = if failed {
        (false, f64::NAN)
    } else {
        accept(
            &mut abs_prob,
            &index,
            &positions,
            &radii,
            &params,
            degen_out,
        )
    };
    let offset_error = if out_ok {
        offset_problem(sketch, &index, &cons, &positions, &radii)
    } else {
        None
    };
    let out_ok = out_ok && offset_error.is_none();
    if !out_ok {
        // The input itself may already be a solution (to the model
        // tolerance) that iterating could not polish or that a converged
        // but degenerate result would replace, e.g. a sketch solved at
        // another scale; then it is the answer. Not for drags: a drag that
        // cannot be followed fails.
        if fallback {
            // Equations in conflict at the solver's tolerance hold at the
            // input to the model tolerance: those whose unknowns were all
            // eliminated are implied by the eliminations, and the others
            // were diagnosed redundant or not by rank as usual.
            let mut eqs = redundant_eqs.clone();
            eqs.extend(&const_conflicts);
            let redundant = owners(&eqs);
            return Solution {
                status: if redundant.is_empty() {
                    SolveStatus::Solved
                } else {
                    SolveStatus::Redundant
                },
                sketch: sketch.clone(),
                dof,
                entities,
                redundant,
                conflicting: Vec::new(),
                max_residual: input_residual,
                iterations,
            };
        }
        let err = if let Some(why) = offset_error {
            SolveError::Degenerate(why)
        } else if conflicting.is_empty() {
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
            max_residual: input_residual,
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
        max_residual,
        iterations,
    }
}

/// Margin on the degeneracy threshold for accepting a result, so that it
/// is not degenerate when solved again (at a scale that may differ).
const ACCEPT_MARGIN: f64 = 8.0;

/// The final check of a result (design note 0005 §5): every equation holds
/// to the model tolerance (on its branch), no direction a constraint needs
/// is shorter than `degen` (metres), no curve has collapsed at that
/// threshold ([`collapsed`], an arc may end where it starts), and
/// coordinates and radii are within the input limits. Returns whether it
/// passes and the largest absolute residual.
fn accept<K: SketchKey>(
    abs_prob: &mut Problem,
    index: &Index<K>,
    positions: &[DVec2],
    radii: &[f64],
    params: &[f64],
    degen: f64,
) -> (bool, f64) {
    abs_prob.degen = degen;
    let va = absolute_slots(index, positions, radii, params);
    let checks = check_all(abs_prob, &va, LINEAR_TOL, ANGULAR_TOL);
    let ok = collapsed(index, positions, radii, None, degen).is_empty()
        && checks
            .iter()
            .all(|c| matches!(c, Some((r, tol, true)) if r.abs() <= *tol))
        && positions
            .iter()
            .all(|p| p.is_finite() && p.x.abs() <= MAX_MAGNITUDE && p.y.abs() <= MAX_MAGNITUDE)
        && radii
            .iter()
            .all(|r| r.is_finite() && *r > degen && *r <= MAX_MAGNITUDE);
    (ok, max_abs_residual(&checks))
}

/// The first offset that the geometry leaves invalid (design note 0007
/// §10.5), as a message naming it: an image end off its offset curve (a
/// smooth joint whose sources no longer meet tangentially), a piece turned
/// back by too large a distance, or pieces that cross.
fn offset_problem<K: SketchKey>(
    sketch: &Sketch<K>,
    index: &Index<K>,
    cons: &[Con],
    positions: &[DVec2],
    radii: &[f64],
) -> Option<String> {
    for (j, con) in cons.iter().enumerate() {
        let Con::Offset(o) = con else { continue };
        let piece = |c: OCurve| match c {
            OCurve::Line(l) => Piece::Line {
                a: positions[l.s],
                b: positions[l.e],
            },
            OCurve::Arc { c, s, e } => Piece::Arc {
                c: positions[c],
                a: positions[s],
                b: positions[e],
            },
            OCurve::Circle { c, r } => Piece::Circle {
                c: positions[c],
                r: radii[r],
            },
        };
        // Image ends held on a normal: free ends and smooth joints.
        let held: Vec<usize> = o
            .points
            .iter()
            .filter(|(_, _, eq)| matches!(eq, OEq::Normal(..)))
            .map(|(x, ..)| *x)
            .collect();
        let pieces: Vec<OffsetPiece> = o
            .pieces
            .iter()
            .map(|&(src, img, left)| OffsetPiece {
                src: piece(src),
                img: piece(img),
                left,
                ends: img.ends(),
                held: img.ends().map_or((false, false), |(s, e)| {
                    (held.contains(&s), held.contains(&e))
                }),
            })
            .collect();
        let Err(why) = offset::check(&pieces, o.distance, LINEAR_TOL) else {
            continue;
        };
        let id = sketch.constraints[j].id;
        let image = |k: usize| index.keys[o.entities[k].1];
        let source = |k: usize| index.keys[o.entities[k].0];
        return Some(match why {
            Invalid::OffCurve(k, end) => format!(
                "the curves at a smooth joint of an offset no longer meet tangentially, so \
                 their offsets do not meet (offset {id:?}, at the {} of {:?})",
                if end { "end" } else { "start" },
                image(k)
            ),
            Invalid::TurnedBack(k) => format!(
                "the offset distance is too large: the offset of a curve would turn back \
                 or vanish (offset {id:?}, curve {:?})",
                source(k)
            ),
            Invalid::Cross(a, b) => format!(
                "an offset would cross itself (offset {id:?}, curves {:?} and {:?})",
                image(a),
                image(b)
            ),
        });
    }
    None
}

/// The point indices of entity `ent`.
fn entity_points<K: SketchKey>(index: &Index<K>, ent: usize) -> Vec<usize> {
    match index.kinds[ent] {
        EntKind::Point(p) => vec![p],
        EntKind::Line(l) => vec![l.s, l.e],
        EntKind::Circular(Circ::Circle { c, .. }) => vec![c],
        EntKind::Circular(Circ::Arc { c, s, e }) => vec![c, s, e],
        EntKind::Ellipse(el) => {
            let mut v = vec![el.c, el.m, el.n];
            if let Some((s, e)) = el.arc {
                v.extend([s, e]);
            }
            v
        }
        EntKind::Spline(s) => index.splines[s].points.clone(),
    }
}

/// Union–find root with path halving.
fn find(uf: &mut [usize], mut x: usize) -> usize {
    while uf[x] != x {
        uf[x] = uf[uf[x]];
        x = uf[x];
    }
    x
}

/// Joins the sets of `a` and `b` under the smaller root.
fn union(uf: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(uf, a), find(uf, b));
    if ra != rb {
        uf[ra.max(rb)] = ra.min(rb);
    }
}

/// The connected parts of a sketch, as a root per point index: points are
/// joined through the curves they define and the constraints that
/// reference them.
fn point_roots<K: SketchKey>(sketch: &Sketch<K>, index: &Index<K>) -> Vec<usize> {
    fn join(uf: &mut [usize], pts: &[usize]) {
        for w in pts.windows(2) {
            union(uf, w[0], w[1]);
        }
    }
    let mut uf: Vec<usize> = (0..index.positions.len()).collect();
    for ent in 0..index.kinds.len() {
        join(&mut uf, &entity_points(index, ent));
    }
    for c in &sketch.constraints {
        let pts: Vec<usize> = c
            .constraint
            .referenced()
            .iter()
            .filter_map(|k| index.ent.get(k))
            .flat_map(|&ent| entity_points(index, ent))
            .collect();
        join(&mut uf, &pts);
    }
    (0..uf.len()).map(|p| find(&mut uf, p)).collect()
}

/// How a curve has collapsed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Collapse {
    /// A line whose end points meet.
    Length,
    /// A circle or arc without a radius.
    Radius,
    /// An ellipse or elliptical arc with an axis of zero length (the
    /// major axis when `true`).
    Axis(bool),
    /// An arc whose end has come onto its start.
    Closed,
}

impl Collapse {
    /// What happened, for messages.
    fn says(self) -> &'static str {
        match self {
            Collapse::Length => "has zero length",
            Collapse::Radius => "has zero radius",
            Collapse::Axis(_) => "has an axis of zero length",
            Collapse::Closed => "ends where it starts",
        }
    }
}

/// The curves that have collapsed at `positions` and `radii`, in entity
/// order: a line no longer than `degen` (metres), a circle or arc whose
/// radius is no longer, or an arc whose end lies within `degen` of its
/// start where it did not in the input (`input`: the input positions and
/// their threshold; `None` when `positions` are the input, whose arcs may
/// end where they start). Splines are not checked: their control points
/// may coincide.
fn collapsed<K: SketchKey>(
    index: &Index<K>,
    positions: &[DVec2],
    radii: &[f64],
    input: Option<(&[DVec2], f64)>,
    degen: f64,
) -> Vec<(usize, Collapse)> {
    let short = |a: usize, b: usize| (positions[b] - positions[a]).length() <= degen;
    let mut out = Vec::new();
    for (ent, k) in index.kinds.iter().enumerate() {
        let how = match *k {
            EntKind::Line(l) if short(l.s, l.e) => Collapse::Length,
            EntKind::Circular(Circ::Circle { r, .. }) if radii[r].abs() <= degen => {
                Collapse::Radius
            }
            EntKind::Circular(Circ::Arc { c, s, e }) if short(c, s) || short(c, e) => {
                Collapse::Radius
            }
            EntKind::Circular(Circ::Arc { s, e, .. })
                if short(s, e) && input.is_some_and(|(at, d)| (at[e] - at[s]).length() > d) =>
            {
                Collapse::Closed
            }
            EntKind::Ellipse(el) if short(el.c, el.m) => Collapse::Axis(true),
            EntKind::Ellipse(el) if short(el.c, el.n) => Collapse::Axis(false),
            EntKind::Ellipse(Ell {
                arc: Some((s, e)), ..
            }) if short(s, e) && input.is_some_and(|(at, d)| (at[e] - at[s]).length() > d) => {
                Collapse::Closed
            }
            _ => continue,
        };
        out.push((ent, how));
    }
    out
}

/// The quantity that has collapsed, as gradient rows over slots: the
/// vector between the two points that met (two rows), or the radius.
fn collapse_rows<K: SketchKey>(
    index: &Index<K>,
    prob: &Problem,
    ent: usize,
    how: Collapse,
) -> Vec<Vec<(usize, f64)>> {
    let between = |a: usize, b: usize| {
        vec![
            vec![(2 * b, 1.0), (2 * a, -1.0)],
            vec![(2 * b + 1, 1.0), (2 * a + 1, -1.0)],
        ]
    };
    match (index.kinds[ent], how) {
        (EntKind::Line(l), _) => between(l.s, l.e),
        (EntKind::Circular(Circ::Circle { r, .. }), _) => vec![vec![(prob.radius_slot(r), 1.0)]],
        (EntKind::Circular(Circ::Arc { s, e, .. }), Collapse::Closed) => between(s, e),
        (EntKind::Circular(Circ::Arc { c, s, .. }), _) => between(c, s),
        (EntKind::Ellipse(el), Collapse::Axis(major)) => {
            between(el.c, if major { el.m } else { el.n })
        }
        (
            EntKind::Ellipse(Ell {
                arc: Some((s, e)), ..
            }),
            Collapse::Closed,
        ) => between(s, e),
        _ => Vec::new(),
    }
}

/// Relative size below which a certificate coefficient or residual is zero
/// in [`holding`].
const HOLD_TOL: f64 = 1e-8;

/// Which constraints hold a collapsed curve collapsed, to first order at
/// the state reached (scaled slot values `v`), indexed by constraint
/// position.
///
/// The collapsed quantity (`rows`: the gradients of its components) is
/// pinned by the equations `eqs` when its gradient lies in the span of
/// theirs, and the constraints holding it are those with a nonzero
/// coefficient in the minimum-norm combination: removing any other leaves
/// that combination intact. Where only some direction of a collapsed
/// vector is pinned, that direction is used; where none is, nothing holds
/// it. An equation that uses a direction or length that has collapsed as
/// well (it cannot be evaluated there, or has a kink) gives no gradient:
/// its constraint is counted as holding, and so are the constraints that
/// hold that vector collapsed, found the same way.
fn holding(
    prob: &Problem,
    v: &[f64],
    eqs: &[usize],
    rows: &[Vec<(usize, f64)>],
    n_constraints: usize,
) -> Vec<bool> {
    let mut out = vec![false; n_constraints];
    let mut quantities: Vec<Vec<Vec<(usize, f64)>>> = vec![rows.to_vec()];
    let mut used: Vec<usize> = Vec::new();
    let mut grads: Vec<Vec<(usize, f64)>> = Vec::new();
    for &i in eqs {
        let kind = &prob.eqs[i].kind;
        let short = kind.short_vectors(v, prob.degen);
        let mut grad = Vec::new();
        let ok = short.is_empty()
            && !kind.abs_kinked(v, prob.degen)
            && kind
                .eval(v, &prob.splines, prob.degen, Some(&mut grad))
                .is_some();
        if ok {
            used.push(i);
            grads.push(grad);
            continue;
        }
        if let Some(c) = prob.eqs[i].owner {
            out[c] = true;
        }
        for w in short {
            quantities.push(vec![
                w.0.iter().map(|&(p, c)| (2 * p, c)).collect(),
                w.0.iter().map(|&(p, c)| (2 * p + 1, c)).collect(),
            ]);
        }
    }
    let mut col = BTreeMap::new();
    for (s, _) in grads.iter().chain(quantities.iter().flatten()).flatten() {
        let n = col.len();
        col.entry(*s).or_insert(n);
    }
    // Equation gradients as columns: A λ = g asks for a combination.
    let mut a = Matrix::zeros(col.len(), used.len());
    for (k, grad) in grads.iter().enumerate() {
        for &(s, g) in grad {
            a[(col[&s], k)] += g;
        }
    }
    let cod = CompleteOrthogonal::new(&a, RANK_TOL);
    for rows in &quantities {
        for lambda in pinned_combinations(&a, &cod, &col, rows) {
            let big = lambda.iter().fold(0.0f64, |m, l| m.max(l.abs()));
            for (k, l) in lambda.iter().enumerate() {
                if l.abs() > HOLD_TOL * big {
                    if let Some(c) = prob.eqs[used[k]].owner {
                        out[c] = true;
                    }
                }
            }
        }
    }
    out
}

/// For [`holding`]: the minimum-norm combinations `λ` of the columns of `a`
/// (factored as `cod`; `col` maps slots to rows of `a`) that equal the
/// gradient rows `rows` of a collapsed quantity: one per row when all are
/// pinned, else (for a vector) one for its direction closest to pinned if
/// that one is, else none.
fn pinned_combinations(
    a: &Matrix,
    cod: &CompleteOrthogonal,
    col: &BTreeMap<usize, usize>,
    rows: &[Vec<(usize, f64)>],
) -> Vec<Vec<f64>> {
    let mut lambdas = Vec::with_capacity(rows.len());
    let mut resids = Vec::with_capacity(rows.len());
    let mut size = 0.0f64;
    for row in rows {
        let mut g = vec![0.0; col.len()];
        for &(s, c) in row {
            g[col[&s]] += c;
        }
        size = size.max(norm2(&g));
        let lambda = cod.solve(&g);
        resids.push(residual(a, &lambda, &g));
        lambdas.push(lambda);
    }
    let pinned = |r: &[f64]| norm2(r) <= HOLD_TOL * size;
    if resids.iter().all(|r| pinned(r)) {
        return lambdas;
    }
    if resids.len() != 2 {
        return Vec::new();
    }
    // The direction (x, y) of the collapsed vector closest to pinned: the
    // eigenvector of the smaller eigenvalue of the residuals' Gram matrix.
    let (p, q, r) = (
        dot(&resids[0], &resids[0]),
        dot(&resids[0], &resids[1]),
        dot(&resids[1], &resids[1]),
    );
    let mu = 0.5 * (p + r) - (0.25 * (p - r) * (p - r) + q * q).sqrt();
    let (x, y) = if q != 0.0 {
        (mu - r, q)
    } else if p <= r {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };
    let h = x.hypot(y);
    let mix = |u: &[f64], w: &[f64]| -> Vec<f64> {
        u.iter().zip(w).map(|(a, b)| (x * a + y * b) / h).collect()
    };
    if pinned(&mix(&resids[0], &resids[1])) {
        vec![mix(&lambdas[0], &lambdas[1])]
    } else {
        Vec::new()
    }
}

/// Most trial solves spent on blaming the collapse of one part.
const MAX_TRIALS: usize = 32;

/// The state a solve reached, for blaming a collapse.
struct Reached<'a> {
    /// The scaled problem.
    prob: &'a Problem,
    /// Its slot values reached.
    v: &'a [f64],
    /// Collapse threshold of the result (metres).
    degen: f64,
    /// Collapse threshold of the input (metres).
    degen_in: f64,
}

/// Blames a constraint for each connected part of the sketch with a
/// collapsed curve ([`Part::attribute`]). Returns the constraints'
/// positions in priority order, or the first collapsed curve of a part
/// where no single removal helps; and the iterations spent.
fn blame_collapse<K: SketchKey>(
    sketch: &Sketch<K>,
    options: &SolveOptions<K>,
    index: &Index<K>,
    reached: &Reached<'_>,
    collapse: &[(usize, Collapse)],
) -> (Result<Vec<usize>, (usize, Collapse)>, usize) {
    let roots = point_roots(sketch, index);
    let mut culprits: Vec<usize> = Vec::new();
    let mut done: Vec<usize> = Vec::new();
    let mut iterations = 0;
    for &(ent, how) in collapse {
        let root = roots[entity_points(index, ent)[0]];
        if done.contains(&root) {
            continue;
        }
        done.push(root);
        let part = Part {
            sketch,
            options,
            index,
            roots: &roots,
            root,
        };
        let (found, its) = part.attribute(reached, ent, how);
        iterations += its;
        match found {
            Some(j) => culprits.push(j),
            None => return (Err((ent, how)), iterations),
        }
    }
    culprits.sort_unstable();
    (Ok(culprits), iterations)
}

/// A connected part of a sketch: the entities and constraints whose points
/// have root `root` in `roots` ([`point_roots`]).
struct Part<'a, K: SketchKey> {
    sketch: &'a Sketch<K>,
    options: &'a SolveOptions<K>,
    index: &'a Index<K>,
    roots: &'a [usize],
    root: usize,
}

impl<K: SketchKey> Part<'_, K> {
    fn has_entity(&self, ent: usize) -> bool {
        self.roots[entity_points(self.index, ent)[0]] == self.root
    }

    fn has_constraint(&self, j: usize) -> bool {
        let refs = self.sketch.constraints[j].constraint.referenced();
        refs.first()
            .and_then(|k| self.index.ent.get(k))
            .is_some_and(|&ent| self.has_entity(ent))
    }

    /// Finds the constraint to blame for the collapse of curve `ent` in
    /// this part (design note 0005 §9). The part's constraints that hold
    /// the curve collapsed ([`holding`], at the state reached) are tried in
    /// reverse priority order (latest first), at most [`MAX_TRIALS`] of
    /// them: the first whose removal lets the part solve with no curve
    /// collapsed (also at the threshold of the whole result) is returned,
    /// as a position in `sketch.constraints`. Each trial solves only the
    /// part, without this search. Also returns the iterations spent.
    fn attribute(
        &self,
        reached: &Reached<'_>,
        ent: usize,
        how: Collapse,
    ) -> (Option<usize>, usize) {
        let Reached {
            prob,
            v,
            degen,
            degen_in,
        } = *reached;
        let sketch = self.sketch;
        let cons: Vec<usize> = (0..sketch.constraints.len())
            .filter(|&j| self.has_constraint(j))
            .collect();
        let eqs: Vec<usize> = (0..prob.eqs.len())
            .filter(|&i| match (prob.eqs[i].owner, prob.eqs[i].arc) {
                (Some(j), _) => self.has_constraint(j),
                (None, Some(a)) => self.has_entity(a),
                _ => false,
            })
            .collect();
        let rows = collapse_rows(self.index, prob, ent, how);
        let held = holding(prob, v, &eqs, &rows, sketch.constraints.len());
        let candidates: Vec<usize> = cons
            .iter()
            .rev()
            .copied()
            .filter(|&j| held[j])
            .take(MAX_TRIALS)
            .collect();
        if candidates.is_empty() {
            return (None, 0);
        }
        let mut part = Sketch {
            entities: sketch
                .entities
                .iter()
                .enumerate()
                .filter(|&(i, _)| self.has_entity(i))
                .map(|(_, e)| e.clone())
                .collect(),
            constraints: Vec::new(),
        };
        let part_options = SolveOptions {
            drag: self
                .options
                .drag
                .iter()
                .filter(|(k, _)| self.index.ent.get(k).is_some_and(|&e| self.has_entity(e)))
                .copied()
                .collect(),
            max_iterations: self.options.max_iterations,
        };
        let Ok(part_index) = Index::new(&part) else {
            return (None, 0);
        };
        let mut iterations = 0;
        for j in candidates {
            part.constraints = cons
                .iter()
                .filter(|&&c| c != j)
                .map(|&c| sketch.constraints[c].clone())
                .collect();
            let out = solve_inner(&part, &part_options, false);
            iterations += out.iterations;
            if !matches!(out.status, SolveStatus::Solved | SolveStatus::Redundant) {
                continue;
            }
            let Ok(out_index) = Index::new(&out.sketch) else {
                continue;
            };
            let input = Some((part_index.positions.as_slice(), degen_in));
            if collapsed(
                &out_index,
                &out_index.positions,
                &out_index.radii,
                input,
                degen,
            )
            .is_empty()
            {
                return (Some(j), iterations);
            }
        }
        (None, iterations)
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
    let slots: Vec<Vec<usize>> = prob
        .eqs
        .iter()
        .map(|e| e.kind.slots(&prob.splines))
        .collect();
    let mut uf: Vec<usize> = (0..n).collect();
    for ss in &slots {
        for w in ss.windows(2) {
            union(&mut uf, w[0], w[1]);
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
            Status::Explicit => {}
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
        // Rank of the full set: independent eliminations and peeled
        // equations plus the ranks of the reduced components inside it.
        let mut rank = eqs
            .iter()
            .filter(|&&i| matches!(red.status[i], Status::Accepted | Status::Explicit))
            .count();
        for (ci, comp) in red.components.iter().enumerate() {
            if comp
                .eqs
                .first()
                .is_some_and(|&e| root_of_eq[e] == Some(root))
            {
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
