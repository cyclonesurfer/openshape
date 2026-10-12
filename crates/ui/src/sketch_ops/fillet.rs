//! Sketch fillet (design note 0007 §8): rounds the corner where two lines
//! or arcs end with an arc of a given radius, tangent to both.
//!
//! The fillet's centre is where the two curves' offsets by the radius
//! meet (lines offset to parallel lines, circles to concentric circles);
//! of the candidates, the one whose tangent points lie on both curves,
//! continue them smoothly, and are nearest the corner is used.

use std::collections::BTreeSet;
use std::f64::consts::TAU;

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, DistanceKind, Geometry, TangentKind};
use openshape_core_model::{SketchData, SketchEdit, SketchId};
use openshape_sketch_solver::ResolvedCurve;

use super::trim::{Plan, Preview, CUT_GAP, TANGENT_SIN};
use super::{coincidence_classes, meet, perp, Batch, Carrier};

/// One of the two curves at the corner.
#[derive(Clone, Debug)]
struct Side {
    curve: SketchId,
    /// Its end point at the corner.
    near: SketchId,
    /// Its other end point.
    far: SketchId,
    /// True when the corner is the curve's end (not its start).
    corner_is_end: bool,
    /// Where the curve lies.
    carrier: Carrier,
    /// The corner and far end positions.
    v: DVec2,
    f: DVec2,
    /// An arc's sweep (counter-clockwise, radians); unused for lines.
    sweep: f64,
}

impl Side {
    /// How far along the side from the corner `t` lies (metres), when it
    /// lies strictly inside the side.
    fn along(&self, t: DVec2) -> Option<f64> {
        let s = match self.carrier {
            Carrier::Line { .. } => (t - self.v).dot((self.f - self.v).normalize()),
            Carrier::Circle { c, r } => {
                let (a, b) = ((self.v - c).to_angle(), (t - c).to_angle());
                let turn = if self.corner_is_end { a - b } else { b - a };
                turn.rem_euclid(TAU) * r
            }
        };
        let length = match self.carrier {
            Carrier::Line { .. } => (self.f - self.v).length(),
            Carrier::Circle { r, .. } => self.sweep * r,
        };
        (s > CUT_GAP && s < length - CUT_GAP).then_some(s)
    }

    /// The unit direction at `t` pointing towards the corner.
    fn towards_corner(&self, t: DVec2) -> DVec2 {
        match self.carrier {
            Carrier::Line { .. } => (self.v - self.f).normalize(),
            Carrier::Circle { c, .. } => {
                let ccw = perp(t - c).normalize();
                if self.corner_is_end {
                    ccw
                } else {
                    -ccw
                }
            }
        }
    }

    /// The closest point of the carrier to `x`.
    fn foot(&self, x: DVec2) -> DVec2 {
        match self.carrier {
            Carrier::Line { p, d } => p + d * (x - p).dot(d),
            Carrier::Circle { c, r } => c + (x - c).normalize() * r,
        }
    }

    /// The carrier offset by `r` to either side.
    fn offsets(&self, r: f64) -> Vec<Carrier> {
        match self.carrier {
            Carrier::Line { p, d } => vec![
                Carrier::Line {
                    p: p + perp(d) * r,
                    d,
                },
                Carrier::Line {
                    p: p - perp(d) * r,
                    d,
                },
            ],
            Carrier::Circle { c, r: big } => {
                let mut v = vec![Carrier::Circle { c, r: big + r }];
                if big - r > CUT_GAP {
                    v.push(Carrier::Circle { c, r: big - r });
                }
                v
            }
        }
    }
}

/// The fillet geometry: centre, the tangent points on each side, and the
/// counter-clockwise arc from `start` to `end`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilletGeometry {
    /// Centre.
    pub center: DVec2,
    /// Tangent point on the first side.
    pub t1: DVec2,
    /// Tangent point on the second side.
    pub t2: DVec2,
    /// True when the counter-clockwise arc runs from `t1` to `t2`.
    pub t1_first: bool,
}

fn geometry(a: &Side, b: &Side, r: f64) -> Option<FilletGeometry> {
    let mut best: Option<(f64, FilletGeometry)> = None;
    for oa in a.offsets(r) {
        for ob in b.offsets(r) {
            for o in meet(oa, ob) {
                let (t1, t2) = (a.foot(o), b.foot(o));
                let (Some(s1), Some(s2)) = (a.along(t1), b.along(t2)) else {
                    continue;
                };
                // Travelling along the first side towards the corner, the
                // fillet must continue smoothly and arrive on the second
                // side heading away from the corner.
                let into = a.towards_corner(t1);
                let ccw = into.dot(perp(t1 - o)) > 0.0;
                let turn = if ccw { 1.0 } else { -1.0 };
                let arrive = perp(t2 - o).normalize() * turn;
                let leave = -b.towards_corner(t2);
                if arrive.dot(leave) < 0.5 {
                    continue;
                }
                let g = FilletGeometry {
                    center: o,
                    t1,
                    t2,
                    t1_first: ccw,
                };
                if best.is_none_or(|(d, _)| s1 + s2 < d) {
                    best = Some((s1 + s2, g));
                }
            }
        }
    }
    best.map(|(_, g)| g)
}

/// Plans a fillet of radius `radius` (metres) at the corner containing
/// point `corner` (design note 0007 §8). `curves` are the sketch's solved
/// curves.
pub fn fillet(
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    corner: SketchId,
    radius: f64,
) -> Result<Plan, String> {
    if !(radius.is_finite() && radius > CUT_GAP) {
        return Err("The fillet radius must be positive".into());
    }
    let classes = coincidence_classes(sketch);
    let at_corner = |p: SketchId| {
        p == corner
            || matches!((classes.get(&p), classes.get(&corner)), (Some(x), Some(y)) if x == y)
    };
    // The curves ending at the corner.
    let mut sides = Vec::new();
    for e in &sketch.entities {
        let (start, end) = match &e.geometry {
            Geometry::Line { start, end } | Geometry::Arc { start, end, .. } => (*start, *end),
            Geometry::ControlSpline { control: p, .. }
            | Geometry::InterpolatedSpline { through: p, .. } => {
                if p.first().is_some_and(|x| at_corner(*x))
                    || p.last().is_some_and(|x| at_corner(*x))
                {
                    return Err("Fillets join lines and arcs; a spline ends here".into());
                }
                continue;
            }
            _ => continue,
        };
        for (near, far, corner_is_end) in [(start, end, false), (end, start, true)] {
            if at_corner(near) {
                sides.push((e.id, near, far, corner_is_end));
            }
        }
    }
    // Construction curves ending at a corner of real ones (a centre-point
    // rectangle's diagonals) stay joined to its virtual sharp.
    if sides.len() > 2 {
        let real: Vec<_> = sides
            .iter()
            .copied()
            .filter(|(id, ..)| sketch.entity(*id).is_some_and(|e| !e.construction))
            .collect();
        if real.len() == 2 {
            sides = real;
        }
    }
    let [(c1, n1, f1, e1), (c2, n2, f2, e2)] = sides[..] else {
        return Err(format!(
            "A fillet needs exactly two lines or arcs ending at the point; {} end here",
            sides.len()
        ));
    };
    if c1 == c2 {
        return Err("Both ends of one curve meet here; there is no corner to round".into());
    }
    if n1 == n2 {
        return Err("The two curves share one end point and cannot be pulled apart".into());
    }
    let side = |curve: SketchId, near: SketchId, far: SketchId, corner_is_end: bool| {
        let shape = curves
            .iter()
            .find(|(id, _)| *id == curve)
            .map(|(_, c)| c)
            .ok_or("the curve no longer exists")?;
        let v = sketch.point(near).ok_or("missing end point")?;
        let f = sketch.point(far).ok_or("missing end point")?;
        let (carrier, sweep) = match shape {
            ResolvedCurve::Line { start, end } => {
                let d = (*end - *start).normalize_or_zero();
                if d == DVec2::ZERO {
                    return Err("A line at the corner has no length".to_string());
                }
                (Carrier::Line { p: *start, d }, 0.0)
            }
            ResolvedCurve::Arc {
                center,
                radius,
                sweep_angle,
                ..
            } => (
                Carrier::Circle {
                    c: *center,
                    r: *radius,
                },
                *sweep_angle,
            ),
            _ => return Err("Fillets join lines and arcs".to_string()),
        };
        Ok(Side {
            curve,
            near,
            far,
            corner_is_end,
            carrier,
            v,
            f,
            sweep,
        })
    };
    let a = side(c1, n1, f1, e1)?;
    let b = side(c2, n2, f2, e2)?;
    // A smooth joint has no corner.
    let (into, out) = (a.towards_corner(a.v), -b.towards_corner(b.v));
    if into.perp_dot(out).abs() < TANGENT_SIN {
        return Err(if into.dot(out) > 0.0 {
            "The curves already meet smoothly here; there is no corner to round".into()
        } else {
            "The curves fold back on each other here".into()
        });
    }
    let g = geometry(&a, &b, radius)
        .ok_or("The radius is too large for these curves (the fillet would not fit on them)")?;
    plan(sketch, &a, &b, &g, radius, &at_corner)
}

fn plan(
    sketch: &SketchData,
    a: &Side,
    b: &Side,
    g: &FilletGeometry,
    radius: f64,
    at_corner: &dyn Fn(SketchId) -> bool,
) -> Result<Plan, String> {
    let mut delete: BTreeSet<SketchId> = BTreeSet::new();
    let mut removed = 0;
    let mut batch = Batch::default();
    let mut pending: Vec<Constraint<SketchId>> = Vec::new();
    let ends = [a.near, b.near];
    let is_line = |s: &Side| matches!(s.carrier, Carrier::Line { .. });
    // What has to move to the virtual sharp, if anything.
    let mut to_sharp: Vec<(SketchId, Constraint<SketchId>)> = Vec::new();
    let mut lengths: Vec<(SketchId, SketchId, f64)> = Vec::new();
    let mut joined: Vec<SketchId> = Vec::new();
    for c in &sketch.constraints {
        let refs = c.constraint.referenced();
        match &c.constraint {
            // The corner's own coincidences dissolve; a point that is not
            // a curve end stays joined to the virtual sharp.
            Constraint::Coincident { a: p, b: q } if at_corner(*p) && at_corner(*q) => {
                delete.insert(c.id);
                for x in [*p, *q] {
                    if !ends.contains(&x) && !joined.contains(&x) {
                        joined.push(x);
                    }
                }
            }
            Constraint::Length { line, value } if *line == a.curve || *line == b.curve => {
                delete.insert(c.id);
                let far = if *line == a.curve { a.far } else { b.far };
                lengths.push((c.id, far, *value));
            }
            Constraint::Equal { a: x, b: y }
                if (is_line(a) && (*x == a.curve || *y == a.curve))
                    || (is_line(b) && (*x == b.curve || *y == b.curve)) =>
            {
                if delete.insert(c.id) {
                    removed += 1;
                }
            }
            Constraint::Midpoint { line, .. } if *line == a.curve || *line == b.curve => {
                if delete.insert(c.id) {
                    removed += 1;
                }
            }
            _ if refs.iter().any(|k| ends.contains(k)) => {
                delete.insert(c.id);
                to_sharp.push((c.id, c.constraint.clone()));
            }
            _ => {}
        }
    }
    let sharp = if to_sharp.is_empty() && lengths.is_empty() && joined.is_empty() {
        None
    } else {
        let s = batch.point(a.v);
        pending.push(Constraint::PointOnCurve {
            point: s,
            curve: a.curve,
        });
        pending.push(Constraint::PointOnCurve {
            point: s,
            curve: b.curve,
        });
        Some(s)
    };
    if let Some(s) = sharp {
        for (_, c) in to_sharp {
            let m = c.map_keys(|k| if ends.contains(&k) { s } else { k });
            let r = m.referenced();
            // A constraint between the two ends says nothing once both are
            // the sharp; one putting an end on the other curve repeats what
            // keeps the sharp there.
            let collapsed = r.iter().filter(|k| **k == s).count() > 1;
            let repeats = matches!(m, Constraint::PointOnCurve { point, curve }
                if point == s && (curve == a.curve || curve == b.curve));
            if collapsed {
                removed += 1;
            } else if !repeats {
                pending.push(m);
            }
        }
        for (_, far, value) in lengths {
            pending.push(Constraint::Distance {
                a: s,
                b: far,
                kind: DistanceKind::Aligned,
                value,
            });
        }
        for p in joined {
            pending.push(Constraint::Coincident { a: s, b: p });
        }
    }
    // The fillet arc, joined and tangent to both sides.
    let construction = [a.curve, b.curve]
        .iter()
        .all(|c| sketch.entity(*c).is_some_and(|e| e.construction));
    let center = batch.point(g.center);
    let (first, second) = if g.t1_first { (a, b) } else { (b, a) };
    let (p_first, p_second) = if g.t1_first {
        (g.t1, g.t2)
    } else {
        (g.t2, g.t1)
    };
    let start = batch.point(p_first);
    let end = batch.point(p_second);
    let arc = batch.entity(Geometry::Arc { center, start, end }, construction);
    pending.push(Constraint::Coincident {
        a: start,
        b: first.near,
    });
    pending.push(Constraint::Coincident {
        a: end,
        b: second.near,
    });
    for s in [first, second] {
        pending.push(Constraint::Tangent {
            a: s.curve,
            b: arc,
            kind: TangentKind::AtPoint { point: s.near },
        });
    }
    pending.push(Constraint::Radius {
        curve: arc,
        value: radius,
    });
    for c in pending {
        batch.constraint(c);
    }
    let mut edits = Vec::new();
    if !delete.is_empty() {
        edits.push(SketchEdit::Delete {
            items: delete.into_iter().collect(),
        });
    }
    edits.push(SketchEdit::MovePoints {
        points: vec![(a.near, g.t1), (b.near, g.t2)],
    });
    edits.extend(batch.edits);
    let (s, e) = (p_first - g.center, p_second - g.center);
    let start_angle = s.to_angle();
    Ok(Plan {
        edits,
        removed,
        preview: Preview::Piece(ResolvedCurve::Arc {
            center: g.center,
            radius,
            start_angle,
            sweep_angle: (e.to_angle() - start_angle).rem_euclid(TAU),
        }),
    })
}
