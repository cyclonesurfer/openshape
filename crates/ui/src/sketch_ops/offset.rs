//! Offset (design note 0007 §10): a chain of lines and arcs, or a circle,
//! moved a distance to one side and held there by one derived offset
//! constraint.
//!
//! The chain is the clicked curve and the curves joined to it end to end
//! through points where exactly two curves meet. Each piece's offset is its
//! parallel line or concentric arc. Consecutive offsets share one point: at
//! a smooth joint it lies on the joint's normal, at a corner where the
//! extended offsets meet (the meeting nearest the joint, after Tiller &
//! Hanson), and a free end lies on its end's normal. Image arcs and circles
//! share their source's centre. The images get no constraints of their own:
//! the offset constraint fixes every image point, so they add no degrees of
//! freedom, and its distance is a dimension. Whether the result is valid
//! (no piece turned back, no crossing) is the solver's call, made on the
//! committed sketch.

use std::collections::{BTreeMap, BTreeSet};

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{SketchData, SketchId};
use openshape_sketch_solver::ResolvedCurve;

use super::trim::{Plan, Preview};
use super::{coincidence_classes, leaving_direction, meet, perp, Batch, Carrier, SAME_POINT};

/// Joints whose tangents differ by less than this (sine of the angle) are
/// smooth: the sources meet tangentially, held so by a constraint.
pub const SMOOTH_SIN: f64 = 1e-7;

/// A chain of curves joined end to end, in order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Chain {
    /// The curves, each with `true` when the chain runs through it from its
    /// end to its start.
    pub pieces: Vec<(SketchId, bool)>,
    /// The chain closes: its last curve ends where its first starts.
    pub closed: bool,
}

/// A point the chain may pass through: a coincidence class, or a point in
/// none.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Node {
    Class(usize),
    Point(SketchId),
}

/// The ends of a line, arc, elliptical arc, conic or spline.
fn ends(g: &Geometry<SketchId>) -> Option<(SketchId, SketchId)> {
    match g {
        Geometry::Line { start, end }
        | Geometry::Arc { start, end, .. }
        | Geometry::EllipticalArc { start, end, .. }
        | Geometry::Conic { start, end, .. } => Some((*start, *end)),
        Geometry::ControlSpline { control: p, .. }
        | Geometry::InterpolatedSpline { through: p, .. } => Some((*p.first()?, *p.last()?)),
        _ => None,
    }
}

/// The chain through `start` (design note 0007 §10.1): curves of the same
/// kind (construction or not) joined end to end through points where
/// exactly two of them end. Splines, elliptical arcs and conics end a chain
/// (they cannot be offset); a circle is a chain of its own. The chain runs in `start`'s direction.
pub fn chain(sketch: &SketchData, start: SketchId) -> Result<Chain, String> {
    let e = sketch.entity(start).ok_or("The curve no longer exists")?;
    match e.geometry {
        Geometry::Circle { .. } => {
            return Ok(Chain {
                pieces: vec![(start, false)],
                closed: true,
            })
        }
        Geometry::Line { .. } | Geometry::Arc { .. } => {}
        Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. } => {
            return Err("A spline cannot be offset: its offset is not a spline".into())
        }
        Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } => {
            return Err("An ellipse cannot be offset: its offset is not an ellipse".into())
        }
        Geometry::Conic { .. } => {
            return Err("A conic cannot be offset: its offset is not a conic".into())
        }
        Geometry::Point { .. } => return Err("Click a line, arc or circle to offset".into()),
        Geometry::Text { .. } => return Err("Text cannot be offset".into()),
    }
    let construction = e.construction;
    let classes = coincidence_classes(sketch);
    let node = |p: SketchId| classes.get(&p).map_or(Node::Point(p), |c| Node::Class(*c));
    // The curve ends at each point: (curve, true for its end).
    let mut at: BTreeMap<Node, Vec<(SketchId, bool)>> = BTreeMap::new();
    let mut end_points: BTreeMap<SketchId, (SketchId, SketchId)> = BTreeMap::new();
    for e in &sketch.entities {
        if e.construction != construction {
            continue;
        }
        let Some((a, b)) = ends(&e.geometry) else {
            continue;
        };
        end_points.insert(e.id, (a, b));
        at.entry(node(a)).or_default().push((e.id, false));
        at.entry(node(b)).or_default().push((e.id, true));
    }
    let offsettable = |id: SketchId| {
        matches!(
            sketch.entity(id).map(|e| &e.geometry),
            Some(Geometry::Line { .. }) | Some(Geometry::Arc { .. })
        )
    };
    // The one other curve ending where `curve` does (at its end or start).
    let next = |curve: SketchId, at_end: bool| -> Option<(SketchId, bool)> {
        let (a, b) = end_points[&curve];
        let p = if at_end { b } else { a };
        let others: Vec<(SketchId, bool)> = at[&node(p)]
            .iter()
            .copied()
            .filter(|&(k, e)| !(k == curve && e == at_end))
            .collect();
        match others.as_slice() {
            [(k, e)] if offsettable(*k) => Some((*k, *e)),
            _ => None,
        }
    };
    let mut seen = BTreeSet::from([start]);
    let mut forward = vec![(start, false)];
    let mut closed = false;
    // Forwards from the start curve's end: a curve entered at its end is
    // run backwards.
    let mut cur = (start, true);
    while let Some((k, entered_at_end)) = next(cur.0, cur.1) {
        if k == start {
            closed = !entered_at_end;
            break;
        }
        if !seen.insert(k) {
            break;
        }
        forward.push((k, entered_at_end));
        cur = (k, !entered_at_end);
    }
    let mut backward = Vec::new();
    if !closed {
        // Backwards from its start: a curve left at its start is run
        // backwards.
        let mut cur = (start, false);
        while let Some((k, left_at_end)) = next(cur.0, cur.1) {
            if !seen.insert(k) {
                break;
            }
            backward.push((k, !left_at_end));
            cur = (k, !left_at_end);
        }
    }
    backward.reverse();
    backward.extend(forward);
    Ok(Chain {
        pieces: backward,
        closed,
    })
}

/// The distance from `at` to a curve and whether `at` is on its left (a
/// line's left in its direction, an arc's or circle's inside).
pub fn side_and_distance(curve: &ResolvedCurve, at: DVec2) -> Option<(f64, bool)> {
    match *curve {
        ResolvedCurve::Line { start, end } => {
            let u = (end - start).try_normalize()?;
            let s = perp(u).dot(at - start);
            Some((s.abs(), s > 0.0))
        }
        ResolvedCurve::Arc { center, radius, .. } | ResolvedCurve::Circle { center, radius } => {
            let r = (at - center).length();
            Some(((r - radius).abs(), r < radius))
        }
        ResolvedCurve::Spline { .. }
        | ResolvedCurve::Ellipse { .. }
        | ResolvedCurve::EllipticalArc { .. }
        | ResolvedCurve::Conic { .. } => None,
    }
}

/// One piece of the chain, solved.
struct Leg {
    id: SketchId,
    /// Run from its end to its start.
    rev: bool,
    /// The chain's start and end of this piece.
    v0: DVec2,
    v1: DVec2,
    /// The direction of travel there (unit).
    t0: DVec2,
    t1: DVec2,
    /// Its offset.
    carrier: Carrier,
    /// The offset is on the piece's own left.
    own_left: bool,
}

/// An offset ready to commit.
#[derive(Clone, Debug, PartialEq)]
pub struct Offsetted {
    /// The batch and its preview.
    pub plan: Plan,
    /// The offset constraint's id within the batch.
    pub constraint: SketchId,
    /// Curves offset.
    pub pieces: usize,
}

impl Offsetted {
    /// The status-bar message after committing.
    pub fn message(&self) -> String {
        match self.pieces {
            1 => "Offset 1 curve".to_string(),
            n => format!("Offset {n} curves"),
        }
    }
}

/// Offsets `chain` by `distance` to its left (`left`, in the direction of
/// its first curve as clicked) or right. `curves` are the sketch's solved
/// curves.
pub fn offset(
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    chain: &Chain,
    distance: f64,
    left: bool,
) -> Result<Offsetted, String> {
    let d = distance;
    if !(d.is_finite() && d > SAME_POINT) {
        return Err("The offset distance must be positive".into());
    }
    let resolved = |id: SketchId| {
        curves
            .iter()
            .find(|(k, _)| *k == id)
            .map(|(_, c)| c.clone())
            .ok_or_else(|| "A curve of the chain has no solved shape".to_string())
    };
    let sgn = if left { 1.0 } else { -1.0 };
    let mut batch = Batch::default();
    let mut preview = Vec::new();
    let mut pairs = Vec::new();

    // A circle: a concentric circle sharing its centre.
    if let [(id, _)] = chain.pieces.as_slice() {
        if let Some(Geometry::Circle { center, .. }) = sketch.entity(*id).map(|e| &e.geometry) {
            let ResolvedCurve::Circle { center: c, radius } = resolved(*id)? else {
                return Err("The circle has no solved shape".into());
            };
            let r = radius - sgn * d;
            if r <= SAME_POINT {
                return Err("The distance is larger than the circle's radius".into());
            }
            let construction = sketch.entity(*id).is_some_and(|e| e.construction);
            let img = batch.entity(
                Geometry::Circle {
                    center: *center,
                    radius: r,
                },
                construction,
            );
            pairs.push((*id, img, left));
            preview.push(ResolvedCurve::Circle {
                center: c,
                radius: r,
            });
            let constraint = batch.constraint(Constraint::Offset {
                distance: d,
                pieces: pairs,
                smooth: Vec::new(),
            });
            return Ok(Offsetted {
                plan: Plan {
                    edits: batch.edits,
                    removed: 0,
                    preview: Preview::Shapes {
                        curves: preview,
                        points: Vec::new(),
                    },
                },
                constraint,
                pieces: 1,
            });
        }
    }

    let mut legs = Vec::with_capacity(chain.pieces.len());
    for &(id, rev) in &chain.pieces {
        let curve = resolved(id)?;
        let (own_s, own_e) = match curve {
            ResolvedCurve::Line { start, end } => (start, end),
            ResolvedCurve::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => (
                center + DVec2::from_angle(start_angle) * radius,
                center + DVec2::from_angle(start_angle + sweep_angle) * radius,
            ),
            _ => return Err("Only lines and arcs can be offset in a chain".into()),
        };
        let leave = |at_end: bool| {
            leaving_direction(&curve, at_end)
                .ok_or_else(|| "A curve of the chain has no length".to_string())
        };
        let (v0, v1) = if rev { (own_e, own_s) } else { (own_s, own_e) };
        let t0 = -leave(rev)?;
        let t1 = leave(!rev)?;
        let own_left = left != rev;
        let carrier = match curve {
            ResolvedCurve::Arc { center, radius, .. } => {
                let r = if own_left { radius - d } else { radius + d };
                if r <= SAME_POINT {
                    return Err(
                        "The distance is larger than the radius of an arc in the chain".into(),
                    );
                }
                Carrier::Circle { c: center, r }
            }
            _ => Carrier::Line {
                p: v0 + perp(t0) * (d * sgn),
                d: t0,
            },
        };
        legs.push(Leg {
            id,
            rev,
            v0,
            v1,
            t0,
            t1,
            carrier,
            own_left,
        });
    }
    let n = legs.len();
    // Where consecutive offsets meet: joint k is between legs k − 1 and k
    // (for a closed chain, joint 0 is between the last leg and the first).
    let join = |a: &Leg, b: &Leg| -> Result<(DVec2, bool), String> {
        let v = a.v1;
        if a.t1.dot(b.t0) > 0.0 && a.t1.perp_dot(b.t0).abs() < SMOOTH_SIN {
            return Ok((v + perp(a.t1) * (d * sgn), true));
        }
        meet(a.carrier, b.carrier)
            .into_iter()
            .min_by(|x, y| (*x - v).length().total_cmp(&(*y - v).length()))
            .map(|x| (x, false))
            .ok_or_else(|| {
                "The offsets of two curves meeting at a corner do not meet: try a smaller distance"
                    .to_string()
            })
    };
    let mut joints: Vec<(DVec2, bool)> = Vec::new();
    if chain.closed {
        for k in 0..n {
            joints.push(join(&legs[(k + n - 1) % n], &legs[k])?);
        }
    } else {
        joints.push((legs[0].v0 + perp(legs[0].t0) * (d * sgn), false));
        for k in 1..n {
            joints.push(join(&legs[k - 1], &legs[k])?);
        }
        joints.push((legs[n - 1].v1 + perp(legs[n - 1].t1) * (d * sgn), false));
    }
    let points: Vec<SketchId> = joints.iter().map(|(at, _)| batch.point(*at)).collect();
    let smooth: Vec<SketchId> = joints
        .iter()
        .zip(&points)
        .filter(|((_, s), _)| *s)
        .map(|(_, p)| *p)
        .collect();
    for (k, leg) in legs.iter().enumerate() {
        let (a, b) = (k, if chain.closed { (k + 1) % n } else { k + 1 });
        let (s, e) = if leg.rev { (b, a) } else { (a, b) };
        let entity = sketch
            .entity(leg.id)
            .ok_or("A curve of the chain no longer exists")?;
        let (ps, pe) = (joints[s].0, joints[e].0);
        let (geometry, shape) = match entity.geometry {
            Geometry::Line { .. } => (
                Geometry::Line {
                    start: points[s],
                    end: points[e],
                },
                ResolvedCurve::Line { start: ps, end: pe },
            ),
            Geometry::Arc { center, .. } => {
                let Carrier::Circle { c, r } = leg.carrier else {
                    unreachable!("an arc's offset is a circle")
                };
                let (a0, a1) = ((ps - c).to_angle(), (pe - c).to_angle());
                let mut sweep = (a1 - a0).rem_euclid(std::f64::consts::TAU);
                if sweep == 0.0 {
                    sweep = std::f64::consts::TAU;
                }
                (
                    Geometry::Arc {
                        center,
                        start: points[s],
                        end: points[e],
                    },
                    ResolvedCurve::Arc {
                        center: c,
                        radius: r,
                        start_angle: a0,
                        sweep_angle: sweep,
                    },
                )
            }
            _ => return Err("Only lines and arcs can be offset in a chain".into()),
        };
        let img = batch.entity(geometry, entity.construction);
        pairs.push((leg.id, img, leg.own_left));
        preview.push(shape);
    }
    let constraint = batch.constraint(Constraint::Offset {
        distance: d,
        pieces: pairs,
        smooth,
    });
    Ok(Offsetted {
        plan: Plan {
            edits: batch.edits,
            removed: 0,
            preview: Preview::Shapes {
                curves: preview,
                points: Vec::new(),
            },
        },
        constraint,
        pieces: n,
    })
}
