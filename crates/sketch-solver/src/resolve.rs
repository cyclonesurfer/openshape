//! Explicit curve geometry for the kernel, and measurement of dimensions.

use std::f64::consts::{PI, TAU};
use std::fmt;

use glam::DVec2;

use crate::eq::signed_angle;
use crate::index::{Circ, Con, EntKind, Index, LineRef};
use crate::model::{Constraint, DistanceKind, Sketch, SketchKey};
use crate::{ResolvedCurve, SolveError};

/// The explicit geometry of every curve entity, in entity order.
pub(crate) fn resolve_curves<K: SketchKey>(
    sketch: &Sketch<K>,
) -> Result<Vec<(K, bool, ResolvedCurve)>, SolveError> {
    let index = Index::new(sketch)?;
    let pos = |p: usize| index.positions[p];
    let mut out = Vec::new();
    for (i, e) in sketch.entities.iter().enumerate() {
        let id = e.id;
        let curve = match index.kinds[i] {
            EntKind::Point(_) => continue,
            EntKind::Line(l) => {
                let (s, t) = (pos(l.s), pos(l.e));
                if s == t {
                    return Err(SolveError::Degenerate(format!(
                        "line {id:?} has zero length"
                    )));
                }
                ResolvedCurve::Line { start: s, end: t }
            }
            EntKind::Circular(Circ::Circle { c, r }) => ResolvedCurve::Circle {
                center: pos(c),
                radius: index.radius(r),
            },
            EntKind::Circular(Circ::Arc { c, s, e }) => {
                let center = pos(c);
                let ds = pos(s) - center;
                let de = pos(e) - center;
                let radius = ds.length();
                if !(radius > 0.0) || de == DVec2::ZERO {
                    return Err(SolveError::Degenerate(format!(
                        "arc {id:?} has zero radius"
                    )));
                }
                let a0 = signed_angle(ds.y, ds.x);
                let a1 = signed_angle(de.y, de.x);
                let mut sweep = a1 - a0;
                if sweep < 0.0 {
                    sweep += TAU;
                }
                if !(sweep > 0.0 && sweep < TAU) {
                    return Err(SolveError::Degenerate(format!(
                        "arc {id:?} ends where it starts"
                    )));
                }
                ResolvedCurve::Arc {
                    center,
                    radius,
                    start_angle: a0,
                    sweep_angle: sweep,
                }
            }
            EntKind::Spline(s) => {
                let sp = &index.splines[s];
                let pts: Vec<DVec2> = sp.points.iter().map(|&p| pos(p)).collect();
                ResolvedCurve::Spline {
                    degree: sp.map.degree,
                    knots: sp.map.knots.clone(),
                    points: sp.map.control_points(&pts),
                }
            }
        };
        out.push((id, e.construction, curve));
    }
    Ok(out)
}

/// Debug label used in messages about the constraint being measured.
struct Measured;

impl fmt::Debug for Measured {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("being measured")
    }
}

fn dir(index_pos: &[DVec2], l: LineRef) -> DVec2 {
    index_pos[l.e] - index_pos[l.s]
}

fn nonzero(d: DVec2, what: &str) -> Result<DVec2, SolveError> {
    if d.length() > 0.0 {
        Ok(d)
    } else {
        Err(SolveError::Degenerate(format!("{what} has zero length")))
    }
}

fn radius_of<K: SketchKey>(index: &Index<K>, c: Circ) -> f64 {
    match c {
        Circ::Circle { r, .. } => index.radius(r),
        Circ::Arc { c, s, .. } => (index.positions[s] - index.positions[c]).length(),
    }
}

/// The current value of a dimension, with the conventions of its equation.
pub(crate) fn measure<K: SketchKey>(
    sketch: &Sketch<K>,
    constraint: &Constraint<K>,
) -> Result<f64, SolveError> {
    if constraint.value().is_none() {
        return Err(SolveError::Invalid(
            "measure: the constraint is not a dimension".to_string(),
        ));
    }
    let index = Index::new(sketch)?;
    // The stored value is irrelevant here; validate with a placeholder.
    let mut probe = constraint.clone();
    probe.set_value(1.0);
    let con = index.constraint_dbg(&Measured, &probe, &[])?;
    let p = &index.positions;
    let line_dist = |q: DVec2, l: LineRef| -> Result<f64, SolveError> {
        let d = nonzero(dir(p, l), "the line")?;
        Ok(d.perp_dot(q - p[l.s]).abs() / d.length())
    };
    Ok(match con {
        Con::DistPP(a, b, kind, _) => {
            let d = p[b] - p[a];
            match kind {
                DistanceKind::Aligned => d.length(),
                DistanceKind::Horizontal => d.x.abs(),
                DistanceKind::Vertical => d.y.abs(),
            }
        }
        Con::DistPL(q, l, _) => line_dist(p[q], l)?,
        Con::DistLL(a, b, _) => line_dist(p[a.s], b)?,
        Con::Length(l, _) => dir(p, l).length(),
        Con::Angle(a, b, _) => {
            let da = nonzero(dir(p, a), "the first line")?;
            let db = nonzero(dir(p, b), "the second line")?;
            let th = signed_angle(da.perp_dot(db), da.dot(db));
            debug_assert!(th > -PI && th <= PI);
            th
        }
        Con::Radius(c, _) => radius_of(&index, c),
        Con::Diameter(c, _) => 2.0 * radius_of(&index, c),
        _ => {
            return Err(SolveError::Invalid(
                "measure: the constraint is not a dimension".to_string(),
            ))
        }
    })
}
