//! Explicit curve geometry for the kernel, and measurement of dimensions.

use std::f64::consts::{PI, TAU};
use std::fmt;

use glam::DVec2;

use crate::conic;
use crate::ellipse;
use crate::eq::signed_angle;
use crate::index::{Circ, Con, EntKind, Index, LineRef};
use crate::model::{Constraint, DistanceKind, Geometry, Sketch, SketchKey, TextPiece};
use crate::spline::clamped_uniform_knots;
use crate::{ResolvedCurve, SolveError};

/// The curves of a text's outline (in box units, ADR 0015) placed in the
/// box with corners `[bottom-left, bottom-right, top-right, top-left]`
/// (only the first, second and fourth are used), or `None` when the box
/// has no area. Equal outline points map to equal plane points, so the
/// contours stay exactly closed.
pub fn text_curves(
    corners: [DVec2; 4],
    aspect: f64,
    outline: &[Vec<TextPiece>],
) -> Option<Vec<ResolvedCurve>> {
    let bl = corners[0];
    let across = corners[1] - bl;
    let up = corners[3] - bl;
    let turn = across.perp_dot(up);
    if !(turn.is_finite() && turn != 0.0) {
        return None;
    }
    let map = |q: DVec2| bl + across * (q.x / aspect) + up * q.y;
    Some(
        outline
            .iter()
            .flatten()
            .map(|piece| {
                let pts: Vec<DVec2> = piece.points.iter().map(|&q| map(q)).collect();
                if piece.degree == 1 {
                    ResolvedCurve::Line {
                        start: pts[0],
                        end: pts[1],
                    }
                } else {
                    let degree = piece.degree as usize;
                    ResolvedCurve::Spline {
                        degree,
                        knots: clamped_uniform_knots(pts.len(), degree),
                        points: pts,
                    }
                }
            })
            .collect(),
    )
}

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
                if radius.is_nan() || radius <= 0.0 || de == DVec2::ZERO {
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
            EntKind::Ellipse(el) => {
                let center = pos(el.c);
                let (major, minor) = (pos(el.m) - center, pos(el.n) - center);
                let turn = major.perp_dot(minor);
                if major == DVec2::ZERO || minor == DVec2::ZERO {
                    return Err(SolveError::Degenerate(format!(
                        "ellipse {id:?} has an axis of zero length"
                    )));
                }
                if turn == 0.0 || turn.is_nan() {
                    return Err(SolveError::Degenerate(format!(
                        "ellipse {id:?} has parallel axes"
                    )));
                }
                match el.arc {
                    None => ResolvedCurve::Ellipse {
                        center,
                        major,
                        minor,
                    },
                    Some((s, e)) => {
                        let t0 = ellipse::param_of(center, major, minor, pos(s));
                        let t1 = ellipse::param_of(center, major, minor, pos(e));
                        // Counter-clockwise in the plane is increasing θ
                        // when the minor axis is left of the major one.
                        let sweep = if turn > 0.0 {
                            (t1 - t0).rem_euclid(TAU)
                        } else {
                            -(t0 - t1).rem_euclid(TAU)
                        };
                        if !(sweep.abs() > 0.0 && sweep.abs() < TAU) {
                            return Err(SolveError::Degenerate(format!(
                                "elliptical arc {id:?} ends where it starts"
                            )));
                        }
                        ResolvedCurve::EllipticalArc {
                            center,
                            major,
                            minor,
                            start_param: t0,
                            sweep_param: sweep,
                        }
                    }
                }
            }
            EntKind::Conic(cn) => {
                let (start, control, end) = (pos(cn.p0), pos(cn.p1), pos(cn.p2));
                if start == end {
                    return Err(SolveError::Degenerate(format!(
                        "conic {id:?} ends where it starts"
                    )));
                }
                let turn = (end - start).perp_dot(control - start);
                if turn == 0.0 || turn.is_nan() {
                    return Err(SolveError::Degenerate(format!(
                        "conic {id:?} is flat: its control point is in line with its ends"
                    )));
                }
                ResolvedCurve::Conic {
                    start,
                    control,
                    end,
                    weight: conic::weight(index.rhos[cn.k]),
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
            EntKind::Text(t) => {
                let Geometry::Text { outline, .. } = &e.geometry else {
                    unreachable!("text kinds come from text geometry")
                };
                let text = index.texts[t];
                let corners = text.c.map(pos);
                let curves = text_curves(corners, text.aspect, outline).ok_or_else(|| {
                    SolveError::Degenerate(format!("text {id:?} has a box of zero size"))
                })?;
                out.extend(curves.into_iter().map(|c| (id, e.construction, c)));
                continue;
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
    // The stored value is irrelevant here; validate with a placeholder in
    // range (a ρ lies in (0, 1)).
    let mut probe = constraint.clone();
    probe.set_value(if matches!(constraint, Constraint::Rho { .. }) {
        0.5
    } else {
        1.0
    });
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
        Con::Rho(cn, _) => index.rhos[cn.k],
        _ => {
            return Err(SolveError::Invalid(
                "measure: the constraint is not a dimension".to_string(),
            ))
        }
    })
}
