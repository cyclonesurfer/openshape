//! Profiles of selected regions and point location (design note 0006, §3).

use std::f64::consts::TAU;

use glam::{DVec2, DVec3};
use openshape_nurbs::NurbsCurve;

use super::curve::Crv;
use super::measure::{RegionShape, Side};
use super::{RegionLoop, SketchArrangement};
use crate::{KernelError, ProfileCurve, Tolerance};

fn bad(m: impl std::fmt::Display) -> KernelError {
    KernelError::InvalidInput(m.to_string())
}

/// Half-edge id: `2·edge` forwards, `2·edge + 1` backwards.
fn he(e: usize, fwd: bool) -> usize {
    2 * e + usize::from(!fwd)
}

/// True if `[t0, t1]` is a usable parameter interval of `c` (see the module
/// docs; a circle's interval may wrap once past 2π).
fn in_domain(c: &ProfileCurve, t0: f64, t1: f64) -> bool {
    const SLACK: f64 = 1e-9;
    if !(t0.is_finite() && t1.is_finite() && t0 < t1) {
        return false;
    }
    match c {
        ProfileCurve::Line { .. } => t0 >= 0.0 && t1 <= 1.0,
        ProfileCurve::Arc { sweep_angle, .. } => {
            t0 >= 0.0 && t1 <= sweep_angle.abs().min(TAU) + SLACK
        }
        ProfileCurve::Circle { .. } => (0.0..TAU).contains(&t0) && t1 - t0 <= TAU + SLACK,
        ProfileCurve::Spline {
            degree,
            knots,
            points,
            ..
        } => {
            let n = points.len();
            n > *degree && knots.len() == n + degree + 1 && t0 >= knots[*degree] && t1 <= knots[n]
        }
    }
}

/// Checks that every index in the arrangement is in range and every edge's
/// parameters are usable with `curves`.
fn check(a: &SketchArrangement, curves: &[ProfileCurve]) -> Result<(), KernelError> {
    let nv = a.vertices.len();
    let nr = a.regions.len();
    for (k, e) in a.edges.iter().enumerate() {
        if e.curve >= curves.len()
            || e.start >= nv
            || e.end >= nv
            || !in_domain(&curves[e.curve], e.t0, e.t1)
            || e.left.is_some_and(|r| r >= nr)
            || e.right.is_some_and(|r| r >= nr)
        {
            return Err(bad(format!(
                "arrangement edge {k} does not match the curves"
            )));
        }
    }
    for (r, reg) in a.regions.iter().enumerate() {
        for l in std::iter::once(&reg.outer).chain(&reg.holes) {
            if l.edges.is_empty() || l.edges.iter().any(|&(e, _)| e >= a.edges.len()) {
                return Err(bad(format!("arrangement region {r} has an invalid loop")));
            }
        }
    }
    Ok(())
}

/// The piece `[t0, t1]` of `c` as a profile curve, in the direction of travel.
fn trimmed(c: &ProfileCurve, t0: f64, t1: f64, fwd: bool) -> Result<ProfileCurve, KernelError> {
    let arc = |center: DVec2, radius: f64, a0: f64, sgn: f64| {
        let (s, e) = (a0 + sgn * t0, a0 + sgn * t1);
        let (s, e) = if fwd { (s, e) } else { (e, s) };
        ProfileCurve::Arc {
            center,
            radius,
            start_angle: s,
            sweep_angle: e - s,
        }
    };
    Ok(match c {
        ProfileCurve::Line { start, end } => {
            let p = |t: f64| {
                if t == 0.0 {
                    *start
                } else if t == 1.0 {
                    *end
                } else {
                    *start + (*end - *start) * t
                }
            };
            let (a, b) = if fwd { (p(t0), p(t1)) } else { (p(t1), p(t0)) };
            ProfileCurve::Line { start: a, end: b }
        }
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => arc(
            *center,
            *radius,
            start_angle.rem_euclid(TAU),
            sweep_angle.signum(),
        ),
        ProfileCurve::Circle { center, radius } => arc(*center, *radius, 0.0, 1.0),
        ProfileCurve::Spline {
            degree,
            knots,
            points,
            weights,
        } => {
            let p3: Vec<DVec3> = points.iter().map(|p| p.extend(0.0)).collect();
            let mut n = NurbsCurve::new(*degree, knots.clone(), p3, weights.clone())
                .map_err(|e| bad(format!("spline: {e}")))?;
            let (a, b) = n.domain();
            if t1 < b {
                n = n.try_split(t1).map_err(|e| bad(format!("spline: {e}")))?.0;
            }
            if t0 > a {
                n = n.try_split(t0).map_err(|e| bad(format!("spline: {e}")))?.1;
            }
            if !fwd {
                n = n.reversed();
            }
            ProfileCurve::Spline {
                degree: n.degree(),
                knots: n.knots().as_slice().to_vec(),
                points: n.control_points().iter().map(|p| p.truncate()).collect(),
                weights: if n.is_rational() {
                    Some(n.weights())
                } else {
                    None
                },
            }
        }
    })
}

/// The profile builder (`profile.rs`) merges end points closer than this
/// many tolerances.
const PROFILE_SNAP: f64 = 10.0;

/// The boundary of the union of `regions` as trimmed profile curves (outer
/// loops counter-clockwise, holes clockwise, each loop in order).
pub(super) fn region_profile(
    a: &SketchArrangement,
    curves: &[ProfileCurve],
    regions: &[usize],
    tol: &Tolerance,
) -> Result<Vec<ProfileCurve>, KernelError> {
    if regions.is_empty() {
        return Err(bad("no region selected"));
    }
    let nr = a.regions.len();
    if let Some(&r) = regions.iter().find(|&&r| r >= nr) {
        return Err(bad(format!(
            "region {r} does not exist (the sketch has {nr})"
        )));
    }
    if !(tol.linear.is_finite() && tol.linear > 0.0) {
        return Err(bad("the linear tolerance must be positive"));
    }
    check(a, curves)?;
    let crvs: Vec<Crv> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| Crv::new(i, c, tol.linear))
        .collect::<Result<_, _>>()?;
    let snap = PROFILE_SNAP * tol.linear;
    let mut selected = vec![false; nr];
    for &r in regions {
        selected[r] = true;
    }
    let ne = a.edges.len();
    let edge_start = |h: usize| {
        let e = &a.edges[h / 2];
        if h % 2 == 0 {
            e.start
        } else {
            e.end
        }
    };
    let edge_end = |h: usize| edge_start(h ^ 1);
    // Successor of each half-edge in its (selected) region's loop.
    let mut succ = vec![usize::MAX; 2 * ne];
    for (r, reg) in a.regions.iter().enumerate() {
        if !selected[r] {
            continue;
        }
        for l in std::iter::once(&reg.outer).chain(&reg.holes) {
            let hs: Vec<usize> = l.edges.iter().map(|&(e, f)| he(e, f)).collect();
            for (i, &h) in hs.iter().enumerate() {
                let e = &a.edges[h / 2];
                let face = if h % 2 == 0 { e.left } else { e.right };
                let n = hs[(i + 1) % hs.len()];
                if face != Some(r) || edge_end(h) != edge_start(n) || succ[h] != usize::MAX {
                    return Err(bad(format!(
                        "arrangement region {r} has an inconsistent loop"
                    )));
                }
                succ[h] = n;
            }
        }
    }
    let other_selected = |h: usize| {
        let e = &a.edges[h / 2];
        let other = if h % 2 == 0 { e.right } else { e.left };
        other.is_some_and(|r| selected[r])
    };
    let boundary = |h: usize| succ[h] != usize::MAX && !other_selected(h);
    // Chain boundary half-edges into loops; at a vertex, rotate clockwise
    // through interior half-edges (twin, then its successor) to the next
    // boundary half-edge.
    let mut used = vec![false; 2 * ne];
    let mut out = Vec::new();
    for h0 in 0..2 * ne {
        if used[h0] || !boundary(h0) {
            continue;
        }
        let mut lp: Vec<usize> = Vec::new();
        let mut h = h0;
        loop {
            if used[h] {
                return Err(bad("arrangement loops are inconsistent"));
            }
            used[h] = true;
            lp.push(h);
            let mut n = succ[h];
            let mut guard = 0;
            while !boundary(n) {
                let tw = n ^ 1;
                if succ[tw] == usize::MAX || guard > 2 * ne {
                    return Err(bad("arrangement loops are inconsistent"));
                }
                n = succ[tw];
                guard += 1;
            }
            if n == h0 {
                break;
            }
            h = n;
        }
        let sides: Vec<Side> = lp
            .iter()
            .map(|&h| {
                let e = &a.edges[h / 2];
                Side {
                    curve: e.curve,
                    t0: e.t0,
                    t1: e.t1,
                    fwd: h % 2 == 0,
                }
            })
            .collect();
        for sd in resolve_loop(&crvs, sides, snap)? {
            out.push(trimmed(&curves[sd.curve], sd.t0, sd.t1, sd.fwd)?);
        }
    }
    if out.is_empty() {
        return Err(bad("the selected regions have no boundary"));
    }
    Ok(out)
}

/// Largest distance from `p` to a point of `sd`, bounded through the
/// control points of its Bézier pieces (so never below the true distance).
fn reach(crvs: &[Crv], sd: &Side, p: DVec2) -> f64 {
    crvs[sd.curve]
        .pieces(sd.t0, sd.t1)
        .iter()
        .flat_map(|q| q.bez.hull())
        .map(|q| (q - p).length())
        .fold(0.0, f64::max)
}

/// Where to cut a piece whose ends nearly meet, and the distance from the
/// cut to the nearer end: the parameter midpoint unless a sampled point
/// (piece boundaries included) is more than twice as far from both ends.
fn split_point(crvs: &[Crv], sd: &Side) -> (f64, f64) {
    const SAMPLES: usize = 32;
    let crv = &crvs[sd.curve];
    let (s, e) = (sd.start(crvs), sd.end(crvs));
    let dist = |t: f64| {
        let q = crv.point(t);
        (q - s).length().min((q - e).length())
    };
    let mid = 0.5 * (sd.t0 + sd.t1);
    let mut best = (mid, dist(mid));
    let at_mid = best.1;
    for piece in crv.pieces(sd.t0, sd.t1) {
        for k in 0..SAMPLES {
            let t = piece.param_at(crv, k as f64 / SAMPLES as f64);
            if !(t > sd.t0 && t < sd.t1) {
                continue;
            }
            let d = dist(t);
            if d > best.1 {
                best = (t, d);
            }
        }
    }
    if at_mid >= 0.5 * best.1 {
        (mid, at_mid)
    } else {
        best
    }
}

/// Makes a boundary loop usable by the profile builder, which merges end
/// points closer than `snap`. A piece that stays within `snap` of its start
/// is degenerate to the builder; a run of such pieces that stays within
/// `2·snap` of where it starts and ends within `snap` of it is left out
/// (the neighbouring pieces then meet within `snap`, and the boundary moves
/// by at most `2·snap`). A piece whose own ends are within `snap` (a closed
/// edge, or one that leaves and comes back) is cut in two (see
/// [`split_point`]). Anything else below the resolution (a whole
/// loop, or one that stays within `2·snap` of a point once its tiny pieces
/// are left out, or a longer chain of tiny pieces) is
/// [`KernelError::InvalidInput`]: leaving it out could change the region
/// by more than the resolution.
fn resolve_loop(crvs: &[Crv], lp: Vec<Side>, snap: f64) -> Result<Vec<Side>, KernelError> {
    let n = lp.len();
    let tiny: Vec<bool> = lp
        .iter()
        .map(|sd| reach(crvs, sd, sd.start(crvs)) <= snap)
        .collect();
    let near = |p: DVec2| format!("({:.9}, {:.9})", p.x, p.y);
    if tiny.iter().all(|&t| t) {
        let p = lp[0].start(crvs);
        return Err(bad(format!(
            "a boundary loop of the selected regions near {} is smaller than the profile \
             resolution ({snap:e})",
            near(p)
        )));
    }
    let mut keep = vec![true; n];
    // Runs of tiny pieces, starting after a piece that is not tiny.
    let first = (0..n).find(|&i| !tiny[i]).unwrap_or(0);
    let mut i = 1;
    while i <= n {
        let k = (first + i) % n;
        if !tiny[k] {
            i += 1;
            continue;
        }
        let mut len = 0;
        while len < n && tiny[(k + len) % n] {
            len += 1;
        }
        let p = lp[k].start(crvs);
        let last = &lp[(k + len - 1) % n];
        let within = (0..len).all(|j| reach(crvs, &lp[(k + j) % n], p) <= 2.0 * snap);
        if !within || (last.end(crvs) - p).length() > snap {
            return Err(bad(format!(
                "boundary edges of the selected regions near {} are shorter than the profile \
                 resolution ({snap:e})",
                near(p)
            )));
        }
        for j in 0..len {
            keep[(k + j) % n] = false;
        }
        i += len;
    }
    // What is left can still stay within twice the resolution of where it
    // starts (two arcs a little longer than the resolution crossing at a
    // narrow angle, once a short third side is left out): too small for the
    // builder as a whole.
    let p = lp[first].start(crvs);
    if (0..n)
        .filter(|&i| keep[i])
        .all(|i| reach(crvs, &lp[i], p) <= 2.0 * snap)
    {
        return Err(bad(format!(
            "a boundary loop of the selected regions near {} is smaller than the profile \
             resolution ({snap:e})",
            near(p)
        )));
    }
    let mut out = Vec::with_capacity(n + 1);
    for (sd, kept) in lp.into_iter().zip(keep) {
        if !kept {
            continue;
        }
        if (sd.end(crvs) - sd.start(crvs)).length() <= snap {
            let (m, d) = split_point(crvs, &sd);
            if d <= snap {
                return Err(bad(format!(
                    "a boundary edge of curve {} near {} is smaller than the profile \
                     resolution ({snap:e})",
                    sd.curve,
                    near(sd.start(crvs))
                )));
            }
            let (a, b) = (Side { t1: m, ..sd }, Side { t0: m, ..sd });
            if sd.fwd {
                out.extend([a, b]);
            } else {
                out.extend([b, a]);
            }
        } else {
            out.push(sd);
        }
    }
    Ok(out)
}

fn sides(a: &SketchArrangement, l: &RegionLoop) -> Vec<Side> {
    l.edges
        .iter()
        .map(|&(e, fwd)| {
            let ed = &a.edges[e];
            Side {
                curve: ed.curve,
                t0: ed.t0,
                t1: ed.t1,
                fwd,
            }
        })
        .collect()
}

/// The region whose interior contains `point` farther than the tolerance
/// from its boundary.
pub(super) fn region_containing(
    a: &SketchArrangement,
    curves: &[ProfileCurve],
    point: DVec2,
    tol: &Tolerance,
) -> Option<usize> {
    if !point.is_finite() || !(tol.linear.is_finite() && tol.linear > 0.0) {
        return None;
    }
    check(a, curves).ok()?;
    let crvs: Vec<Crv> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| Crv::new(i, c, tol.linear))
        .collect::<Result<_, _>>()
        .ok()?;
    for (r, reg) in a.regions.iter().enumerate() {
        let mut loops = vec![sides(a, &reg.outer)];
        loops.extend(reg.holes.iter().map(|h| sides(a, h)));
        if RegionShape::new(&crvs, loops).contains(point, tol.linear) {
            return Some(r);
        }
    }
    None
}
