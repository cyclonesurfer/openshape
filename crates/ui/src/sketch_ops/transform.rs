//! Move, rotate and scale (design note 0007 §13, ADR 0014): the selection,
//! or copies of it, mapped by a similarity.
//!
//! A similarity keeps every relation among the moved entities and every
//! angle, and multiplies lengths by its factor, so constraints among the
//! selection are kept with their dimensions scaled; constraints to the
//! sketch frame that a rotation breaks are swapped, converted or removed;
//! dimensions to the rest of the sketch are measured anew; relations to the
//! rest are kept (the rest follows, as in a drag) or removed. The moved
//! geometry is pinned where the similarity puts it while the rest of the
//! sketch is solved, so the selection lands exactly there or the transform
//! is refused.

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::FRAC_PI_2;

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, ConstraintEntry, DistanceKind, Geometry};
use openshape_core_model::{SketchData, SketchEdit, SketchId};
use openshape_sketch_solver::{
    measure, resolve_curves, solve, ResolvedCurve, SolveError, SolveOptions, SolveStatus,
};

use super::trim::{Plan, Preview};
use super::{Batch, SAME_POINT};
use crate::geometry_edits;

/// Turns within this of a whole number of quarter turns are exactly that
/// many (the model's angular tolerance, ADR 0006).
pub const QUARTER_SNAP: f64 = 1e-9;

/// A similarity of the plane with a positive factor: `p ↦ centre +
/// scale·R(angle)(p − centre) + shift`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Similarity {
    /// The point the turn and the factor are about.
    pub centre: DVec2,
    /// The turn (radians, counter-clockwise).
    pub angle: f64,
    /// The factor (> 0).
    pub scale: f64,
    /// The translation applied after the turn and the factor.
    pub shift: DVec2,
}

impl Similarity {
    /// The translation by `v`.
    pub fn translation(v: DVec2) -> Similarity {
        Similarity {
            centre: DVec2::ZERO,
            angle: 0.0,
            scale: 1.0,
            shift: v,
        }
    }

    /// The turn by `angle` about `centre`.
    pub fn rotation(centre: DVec2, angle: f64) -> Similarity {
        Similarity {
            centre,
            angle,
            scale: 1.0,
            shift: DVec2::ZERO,
        }
    }

    /// The scaling by `factor` about `centre`.
    pub fn scaling(centre: DVec2, factor: f64) -> Similarity {
        Similarity {
            centre,
            angle: 0.0,
            scale: factor,
            shift: DVec2::ZERO,
        }
    }

    /// The turn as a whole number of quarter turns (0 to 3) when it is
    /// within [`QUARTER_SNAP`] of one.
    pub fn quarter_turns(&self) -> Option<u8> {
        let k = (self.angle / FRAC_PI_2).round();
        ((self.angle - k * FRAC_PI_2).abs() <= QUARTER_SNAP).then(|| k.rem_euclid(4.0) as u8)
    }

    /// Cosine and sine of the turn, exact for quarter turns.
    fn cos_sin(&self) -> (f64, f64) {
        match self.quarter_turns() {
            Some(0) => (1.0, 0.0),
            Some(1) => (0.0, 1.0),
            Some(2) => (-1.0, 0.0),
            Some(3) => (0.0, -1.0),
            _ => (self.angle.cos(), self.angle.sin()),
        }
    }

    /// The turn, exactly a multiple of a quarter turn when it snaps to one.
    fn turn(&self) -> f64 {
        match self.quarter_turns() {
            Some(k) => f64::from(k) * FRAC_PI_2,
            None => self.angle,
        }
    }

    /// The image of a vector (turned and scaled, not shifted).
    pub fn vector(&self, v: DVec2) -> DVec2 {
        let (c, s) = self.cos_sin();
        DVec2::new(c * v.x - s * v.y, s * v.x + c * v.y) * self.scale
    }

    /// The image of a point (exactly `p + shift` for a translation).
    pub fn point(&self, p: DVec2) -> DVec2 {
        if self.quarter_turns() == Some(0) && self.scale == 1.0 {
            return p + self.shift;
        }
        self.centre + self.vector(p - self.centre) + self.shift
    }

    /// True when nothing moves.
    pub fn is_identity(&self) -> bool {
        self.quarter_turns() == Some(0) && self.scale == 1.0 && self.shift == DVec2::ZERO
    }

    /// Whether the parameters are usable: finite, with a positive factor.
    fn check(&self) -> Result<(), String> {
        if !(self.centre.is_finite() && self.angle.is_finite() && self.shift.is_finite()) {
            return Err("The transform is not finite".into());
        }
        if !(self.scale.is_finite() && self.scale > 0.0) {
            return Err("A scale factor is a positive number".into());
        }
        Ok(())
    }

    /// The image of a solved curve. Similarities keep orientation, so arcs
    /// keep their start and end and ellipses their eccentric anomalies.
    pub fn curve(&self, c: &ResolvedCurve) -> ResolvedCurve {
        match c {
            ResolvedCurve::Line { start, end } => ResolvedCurve::Line {
                start: self.point(*start),
                end: self.point(*end),
            },
            ResolvedCurve::Circle { center, radius } => ResolvedCurve::Circle {
                center: self.point(*center),
                radius: radius * self.scale,
            },
            ResolvedCurve::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => ResolvedCurve::Arc {
                center: self.point(*center),
                radius: radius * self.scale,
                start_angle: start_angle + self.turn(),
                sweep_angle: *sweep_angle,
            },
            ResolvedCurve::Ellipse {
                center,
                major,
                minor,
            } => ResolvedCurve::Ellipse {
                center: self.point(*center),
                major: self.vector(*major),
                minor: self.vector(*minor),
            },
            ResolvedCurve::EllipticalArc {
                center,
                major,
                minor,
                start_param,
                sweep_param,
            } => ResolvedCurve::EllipticalArc {
                center: self.point(*center),
                major: self.vector(*major),
                minor: self.vector(*minor),
                start_param: *start_param,
                sweep_param: *sweep_param,
            },
            ResolvedCurve::Spline {
                degree,
                knots,
                points,
            } => ResolvedCurve::Spline {
                degree: *degree,
                knots: knots.clone(),
                points: points.iter().map(|p| self.point(*p)).collect(),
            },
            ResolvedCurve::Conic {
                start,
                control,
                end,
                weight,
            } => ResolvedCurve::Conic {
                start: self.point(*start),
                control: self.point(*control),
                end: self.point(*end),
                weight: *weight,
            },
        }
    }
}

/// How a transform treats the rest of the sketch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Options {
    /// Keep the originals and transform copies of them.
    pub copy: bool,
    /// Keep relations between the selection and the rest of the sketch
    /// (the rest follows); otherwise they are removed. Copies never get
    /// them.
    pub keep: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            copy: false,
            keep: true,
        }
    }
}

/// A transform ready to commit.
#[derive(Clone, Debug, PartialEq)]
pub struct Transformed {
    /// The batch and its preview.
    pub plan: Plan,
    /// The transformed entities: the sources, or their copies
    /// (batch-local ids), in source order.
    pub result: Vec<SketchId>,
    /// Items transformed (the sources).
    pub items: usize,
    /// True when the originals stay and copies were made.
    pub copy: bool,
    /// Relations to the rest of the sketch removed (Keep constraints off).
    pub detached: usize,
    /// Constraints among the selection removed (or, for copies, not
    /// copied) because the transform cannot keep them: horizontal and
    /// vertical point pairs and distances a turn breaks, and dimensions
    /// that would measure nothing.
    pub lost: usize,
    /// True when a turn removed the first horizontal or vertical line's
    /// constraint, leaving the angle free.
    pub freed: bool,
    /// Horizontal and vertical lines made parallel or perpendicular to
    /// that line.
    pub converted: usize,
}

impl Transformed {
    /// The status-bar message after committing, starting with `verb`
    /// ("Moved", "Rotated", …).
    pub fn message(&self, verb: &str) -> String {
        let mut text = match self.items {
            1 => format!("{verb} 1 item"),
            n => format!("{verb} {n} items"),
        };
        let plural = |n: usize| if n == 1 { "" } else { "s" };
        if self.detached > 0 {
            text.push_str(&format!(
                "; removed {} constraint{} to the rest of the sketch",
                self.detached,
                plural(self.detached)
            ));
        }
        match (self.freed, self.converted) {
            (false, _) => {}
            (true, 0) => text.push_str(
                "; the angle is now free (a horizontal or vertical constraint was removed)",
            ),
            (true, n) => text.push_str(&format!(
                "; the angle is now free, with {n} horizontal or vertical constraint{} made \
                 parallel or perpendicular",
                plural(n)
            )),
        }
        if self.lost > 0 {
            let n = self.lost;
            text.push_str(&if self.copy {
                format!(
                    "; {n} constraint{} not copied, the transform could not keep them",
                    plural(n)
                )
            } else {
                format!("; removed {n} constraint{} it could not keep", plural(n))
            });
        }
        text
    }
}

/// What happens to one constraint.
#[derive(Clone, Debug, PartialEq)]
enum Fate {
    /// Unchanged.
    Keep,
    /// A dimension among the selection with a new value.
    Value(f64),
    /// Replaced by another constraint (deleted and added).
    Replace(Constraint<SketchId>),
    /// Removed.
    Remove(Why),
    /// A dimension between the selection and the rest, measured anew.
    Measure,
}

/// Why a constraint is removed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Why {
    /// It tied the selection to the rest of the sketch (Keep constraints
    /// off).
    Detached,
    /// The transform cannot keep it.
    Lost,
    /// A horizontal or vertical line's constraint, removed to free the
    /// rotation (the others become parallel or perpendicular to it).
    Freed,
}

/// Which way a constraint ties to the moved entities.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Class {
    /// Everything it references moves.
    Internal,
    /// Nothing it references moves.
    External,
    /// Some of what it references moves.
    Boundary,
}

fn class_of(keys: &[SketchId], moved: &BTreeSet<SketchId>) -> Class {
    let n = keys.iter().filter(|k| moved.contains(k)).count();
    if n == 0 {
        Class::External
    } else if n == keys.len() {
        Class::Internal
    } else {
        Class::Boundary
    }
}

/// A derived constraint's anchor (mirror line, direction line, centre) and
/// units (pairs, pieces, copies), as key lists.
fn units(c: &Constraint<SketchId>) -> Option<(Option<SketchId>, Vec<[SketchId; 2]>)> {
    Some(match c {
        Constraint::Mirror { axis, pairs } => {
            (Some(*axis), pairs.iter().map(|&(s, i)| [s, i]).collect())
        }
        Constraint::Offset { pieces, .. } => {
            (None, pieces.iter().map(|&(s, i, _)| [s, i]).collect())
        }
        Constraint::LinearPattern {
            direction: anchor,
            copies,
        }
        | Constraint::CircularPattern {
            centre: anchor,
            copies,
            ..
        } => (
            Some(*anchor),
            copies.iter().map(|&(s, i, _)| [s, i]).collect(),
        ),
        _ => return None,
    })
}

/// The derived constraint keeping only the units for which `keep` holds
/// (an offset also only the smooth points of kept images), or `None` when
/// none is left.
fn keep_units(
    c: &Constraint<SketchId>,
    sketch: &SketchData,
    keep: impl Fn(usize) -> bool,
) -> Option<Constraint<SketchId>> {
    let pick = |n: usize| (0..n).filter(|&k| keep(k)).collect::<Vec<usize>>();
    let out = match c {
        Constraint::Mirror { axis, pairs } => Constraint::Mirror {
            axis: *axis,
            pairs: pick(pairs.len()).into_iter().map(|k| pairs[k]).collect(),
        },
        Constraint::Offset {
            distance,
            pieces,
            smooth,
        } => {
            let pieces: Vec<(SketchId, SketchId, bool)> =
                pick(pieces.len()).into_iter().map(|k| pieces[k]).collect();
            let image_points: BTreeSet<SketchId> = pieces
                .iter()
                .filter_map(|(_, i, _)| sketch.entity(*i))
                .flat_map(|e| e.geometry.referenced_points())
                .collect();
            Constraint::Offset {
                distance: *distance,
                smooth: smooth
                    .iter()
                    .copied()
                    .filter(|p| image_points.contains(p))
                    .collect(),
                pieces,
            }
        }
        Constraint::LinearPattern { direction, copies } => Constraint::LinearPattern {
            direction: *direction,
            copies: pick(copies.len()).into_iter().map(|k| copies[k]).collect(),
        },
        Constraint::CircularPattern {
            centre,
            angle,
            copies,
        } => Constraint::CircularPattern {
            centre: *centre,
            angle: *angle,
            copies: pick(copies.len()).into_iter().map(|k| copies[k]).collect(),
        },
        _ => return None,
    };
    units(&out)
        .is_some_and(|(_, u)| !u.is_empty())
        .then_some(out)
}

/// The constraints decided per constraint, with the rotation's reference
/// line (see [`internal_fate`]).
struct Decider<'a> {
    sketch: &'a SketchData,
    moved: &'a BTreeSet<SketchId>,
    t: Similarity,
    keep: bool,
    /// Under a turn that is not a whole number of quarter turns: the first
    /// horizontal or vertical line among the selection (constraint id,
    /// line, true if horizontal).
    reference: Option<(SketchId, SketchId, bool)>,
}

impl Decider<'_> {
    fn new<'a>(
        sketch: &'a SketchData,
        moved: &'a BTreeSet<SketchId>,
        t: Similarity,
        keep: bool,
    ) -> Decider<'a> {
        let reference = if t.quarter_turns().is_some() {
            None
        } else {
            sketch.constraints.iter().find_map(|c| match c.constraint {
                Constraint::Horizontal { line } if moved.contains(&line) => {
                    Some((c.id, line, true))
                }
                Constraint::Vertical { line } if moved.contains(&line) => Some((c.id, line, false)),
                _ => None,
            })
        };
        Decider {
            sketch,
            moved,
            t,
            keep,
            reference,
        }
    }

    /// A length dimension's fate under the factor.
    fn scaled(&self, v: f64) -> Fate {
        if self.t.scale == 1.0 {
            Fate::Keep
        } else {
            Fate::Value(v * self.t.scale)
        }
    }

    /// What happens to a constraint all of whose entities move (not a
    /// derived one).
    fn internal_fate(&self, id: SketchId, c: &Constraint<SketchId>) -> Fate {
        let quarter = self.t.quarter_turns();
        let swap = matches!(quarter, Some(1 | 3));
        let level = |a: SketchId, b: SketchId, horizontal: bool| match (
            self.sketch.point(a),
            self.sketch.point(b),
        ) {
            (Some(p), Some(q)) => {
                let d = q - p;
                (if horizontal { d.y } else { d.x }).abs() <= SAME_POINT
            }
            _ => false,
        };
        match *c {
            Constraint::Horizontal { line } | Constraint::Vertical { line } => {
                let horizontal = matches!(c, Constraint::Horizontal { .. });
                match (quarter, self.reference) {
                    (Some(_), _) if !swap => Fate::Keep,
                    (Some(_), _) => Fate::Replace(if horizontal {
                        Constraint::Vertical { line }
                    } else {
                        Constraint::Horizontal { line }
                    }),
                    (None, Some((rid, _, _))) if rid == id => Fate::Remove(Why::Freed),
                    (None, Some((_, r, _))) if r == line => Fate::Remove(Why::Lost),
                    (None, Some((_, r, rh))) => Fate::Replace(if rh == horizontal {
                        Constraint::Parallel { a: line, b: r }
                    } else {
                        Constraint::Perpendicular { a: line, b: r }
                    }),
                    (None, None) => Fate::Remove(Why::Lost),
                }
            }
            Constraint::HorizontalPoints { a, b } | Constraint::VerticalPoints { a, b } => {
                let horizontal = matches!(c, Constraint::HorizontalPoints { .. });
                match quarter {
                    Some(_) if !swap => Fate::Keep,
                    Some(_) => Fate::Replace(if horizontal {
                        Constraint::VerticalPoints { a, b }
                    } else {
                        Constraint::HorizontalPoints { a, b }
                    }),
                    None => Fate::Remove(Why::Lost),
                }
            }
            Constraint::Distance { a, b, kind, value } => {
                let s = self.t.scale;
                match (kind, quarter) {
                    (DistanceKind::Aligned, _) => self.scaled(value),
                    (_, Some(_)) if !swap => self.scaled(value),
                    (_, Some(_)) => Fate::Replace(Constraint::Distance {
                        a,
                        b,
                        kind: if kind == DistanceKind::Horizontal {
                            DistanceKind::Vertical
                        } else {
                            DistanceKind::Horizontal
                        },
                        value: value * s,
                    }),
                    (_, None) if level(a, b, kind == DistanceKind::Horizontal) => {
                        Fate::Replace(Constraint::Distance {
                            a,
                            b,
                            kind: DistanceKind::Aligned,
                            value: value * s,
                        })
                    }
                    (_, None) => Fate::Remove(Why::Lost),
                }
            }
            Constraint::Length { value, .. }
            | Constraint::Radius { value, .. }
            | Constraint::Diameter { value, .. } => self.scaled(value),
            Constraint::Fix { point, at } => {
                let to = self.t.point(at);
                if to == at {
                    Fate::Keep
                } else {
                    Fate::Replace(Constraint::Fix { point, at: to })
                }
            }
            _ => Fate::Keep,
        }
    }

    /// What happens to a constraint, in place.
    fn fate(&self, entry: &ConstraintEntry<SketchId>) -> Fate {
        let c = &entry.constraint;
        if let Some((anchor, us)) = units(c) {
            let classes: Vec<Class> = us
                .iter()
                .map(|u| {
                    let keys: Vec<SketchId> = anchor.into_iter().chain(u.iter().copied()).collect();
                    class_of(&keys, self.moved)
                })
                .collect();
            if classes.iter().all(|k| *k == Class::External) {
                return Fate::Keep;
            }
            if classes.iter().all(|k| *k == Class::Internal) {
                return match c {
                    Constraint::Offset { distance, .. } => self.scaled(*distance),
                    _ => Fate::Keep,
                };
            }
            if self.keep {
                return Fate::Keep;
            }
            return match keep_units(c, self.sketch, |k| classes[k] != Class::Boundary) {
                Some(reduced) if reduced == *c => Fate::Keep,
                Some(reduced) => Fate::Replace(reduced),
                None => Fate::Remove(Why::Detached),
            };
        }
        match class_of(&c.referenced(), self.moved) {
            Class::External => Fate::Keep,
            Class::Internal => self.internal_fate(entry.id, c),
            Class::Boundary if c.value().is_some() => Fate::Measure,
            Class::Boundary if self.keep => Fate::Keep,
            Class::Boundary => Fate::Remove(Why::Detached),
        }
    }

    /// The constraint a copy gets for `entry`, with keys mapped by `map`
    /// (the copied entities), or `None`. `Err` carries why one among the
    /// selection is not copied.
    fn copied(
        &self,
        entry: &ConstraintEntry<SketchId>,
        map: &BTreeMap<SketchId, SketchId>,
    ) -> Result<Option<Constraint<SketchId>>, Why> {
        let c = &entry.constraint;
        let mapped = |c: &Constraint<SketchId>| c.map_keys(|k| map.get(&k).copied().unwrap_or(k));
        if let Some((anchor, us)) = units(c) {
            if anchor.is_some_and(|a| !self.moved.contains(&a)) {
                return Ok(None);
            }
            let internal: Vec<bool> = us
                .iter()
                .map(|u| u.iter().all(|k| self.moved.contains(k)))
                .collect();
            let Some(mut kept) = keep_units(c, self.sketch, |k| internal[k]) else {
                return Ok(None);
            };
            if let Constraint::Offset { distance, .. } = &mut kept {
                *distance *= self.t.scale;
            }
            return Ok(Some(mapped(&kept)));
        }
        if class_of(&c.referenced(), self.moved) != Class::Internal {
            return Ok(None);
        }
        match self.internal_fate(entry.id, c) {
            Fate::Keep => Ok(Some(mapped(c))),
            Fate::Value(v) => {
                let mut c = mapped(c);
                c.set_value(v);
                Ok(Some(c))
            }
            Fate::Replace(r) => Ok(Some(mapped(&r))),
            Fate::Remove(why) => Err(why),
            Fate::Measure => Ok(None),
        }
    }
}

/// The entities a transform moves for `sources` (curves and points): the
/// sources, the points of source curves, and every other curve all of
/// whose points move.
pub fn moved_set(sketch: &SketchData, sources: &[SketchId]) -> BTreeSet<SketchId> {
    let mut moved = BTreeSet::new();
    for &id in sources {
        if let Some(e) = sketch.entity(id) {
            moved.insert(id);
            moved.extend(e.geometry.referenced_points());
        }
    }
    for e in &sketch.entities {
        let pts = e.geometry.referenced_points();
        if e.geometry.is_curve() && !pts.is_empty() && pts.iter().all(|p| moved.contains(p)) {
            moved.insert(e.id);
        }
    }
    moved
}

/// Transforms `sources` (curves and points, as `mirror::sources` picks
/// them) by `t`: one batch that moves them (or adds transformed copies),
/// with their constraints adjusted as design note 0007 §13 says.
///
/// `sketch` is the solved sketch, `labels` its placed dimension labels and
/// `curves` its solved curves (for the preview). In place, the moved
/// geometry is pinned at its image while the rest of the sketch is solved;
/// the transform is refused, with the reason, when the rest cannot follow.
pub fn transform(
    sketch: &SketchData,
    labels: &[(SketchId, DVec2)],
    curves: &[(SketchId, ResolvedCurve)],
    sources: &[SketchId],
    t: Similarity,
    options: Options,
) -> Result<Transformed, String> {
    t.check()?;
    let sources: Vec<SketchId> = sources
        .iter()
        .copied()
        .filter(|id| sketch.entity(*id).is_some())
        .collect();
    if sources.is_empty() {
        return Err("Nothing to transform: select curves or points".into());
    }
    if t.is_identity() {
        return Err(if options.copy {
            "The copy would lie on the original".into()
        } else {
            "That leaves the selection where it is".into()
        });
    }
    let moved = moved_set(sketch, &sources);
    let decider = Decider::new(sketch, &moved, t, options.keep);
    if options.copy {
        copy(sketch, labels, curves, &sources, &moved, &decider)
    } else {
        in_place(sketch, labels, &sources, &moved, &decider)
    }
}

/// Lone points (used by no curve) among `ids`.
fn lone_points(sketch: &SketchData, ids: &BTreeSet<SketchId>) -> Vec<SketchId> {
    let used: BTreeSet<SketchId> = sketch
        .entities
        .iter()
        .flat_map(|e| e.geometry.referenced_points())
        .collect();
    sketch
        .entities
        .iter()
        .filter(|e| ids.contains(&e.id) && !e.geometry.is_curve() && !used.contains(&e.id))
        .map(|e| e.id)
        .collect()
}

fn in_place(
    sketch: &SketchData,
    labels: &[(SketchId, DVec2)],
    sources: &[SketchId],
    moved: &BTreeSet<SketchId>,
    decider: &Decider<'_>,
) -> Result<Transformed, String> {
    let t = decider.t;
    // The trial: moved geometry at its image, adjusted constraints, and the
    // moved unknowns pinned there (fixes, radii, ρ), with the boundary
    // dimensions left out until the result is known.
    let mut trial = sketch.clone();
    let mut pins = Vec::new();
    for e in &mut trial.entities {
        if !moved.contains(&e.id) {
            continue;
        }
        match &mut e.geometry {
            Geometry::Point { at } => {
                *at = t.point(*at);
                pins.push(Constraint::Fix {
                    point: e.id,
                    at: *at,
                });
            }
            Geometry::Circle { radius, .. } => {
                *radius *= t.scale;
                pins.push(Constraint::Radius {
                    curve: e.id,
                    value: *radius,
                });
            }
            Geometry::Conic { rho, .. } => pins.push(Constraint::Rho {
                curve: e.id,
                value: *rho,
            }),
            _ => {}
        }
    }
    let mut out = Transformed {
        plan: Plan {
            edits: Vec::new(),
            removed: 0,
            preview: Preview::Shapes {
                curves: Vec::new(),
                points: Vec::new(),
            },
        },
        result: sources.to_vec(),
        items: sources.len(),
        copy: false,
        detached: 0,
        lost: 0,
        freed: false,
        converted: 0,
    };
    let mut deletes = Vec::new();
    let mut adds = Vec::new();
    let mut values = Vec::new();
    let mut measured = Vec::new();
    let mut kept_dims = BTreeSet::new();
    trial.constraints.clear();
    for entry in &sketch.constraints {
        let fate = decider.fate(entry);
        let mut c = entry.constraint.clone();
        match fate {
            Fate::Keep => {
                if c.value().is_some() && class_of(&c.referenced(), moved) == Class::Internal {
                    kept_dims.insert(entry.id);
                }
            }
            Fate::Value(v) => {
                c.set_value(v);
                values.push((entry.id, v));
                kept_dims.insert(entry.id);
            }
            Fate::Replace(r) => {
                if matches!(
                    (&entry.constraint, &r),
                    (
                        Constraint::Horizontal { .. } | Constraint::Vertical { .. },
                        Constraint::Parallel { .. } | Constraint::Perpendicular { .. }
                    )
                ) {
                    out.converted += 1;
                } else if units(&r).is_some() {
                    out.detached += 1;
                }
                deletes.push(entry.id);
                adds.push(r.clone());
                c = r;
            }
            Fate::Remove(why) => {
                match why {
                    Why::Detached => out.detached += 1,
                    Why::Lost => out.lost += 1,
                    Why::Freed => out.freed = true,
                }
                deletes.push(entry.id);
                continue;
            }
            Fate::Measure => {
                measured.push(entry.clone());
                continue;
            }
        }
        trial.constraints.push(ConstraintEntry {
            id: entry.id,
            constraint: c,
        });
    }
    for (k, c) in pins.into_iter().enumerate() {
        trial.constraints.push(ConstraintEntry {
            id: SketchId::new_in_batch(k as u32),
            constraint: c,
        });
    }
    let solution = solve(&trial, &SolveOptions::default());
    if let SolveStatus::Failed(e) = &solution.status {
        return Err(match e {
            SolveError::Conflicting | SolveError::NotConverged if decider.keep => {
                "Constraints to the rest of the sketch hold the selection: turn off Keep \
                 constraints to detach it"
                    .to_string()
            }
            SolveError::Conflicting | SolveError::NotConverged => {
                "Constraints on geometry sharing points with the selection hold it".to_string()
            }
            // A similarity keeps the selection's own curves as they are, so
            // what collapses is a curve reaching from it to the rest.
            SolveError::Degenerate(_) => {
                "That would collapse a curve joined to the selection to nothing".to_string()
            }
            e => format!("The result would not solve: {e}"),
        });
    }
    let solved = solution.sketch;
    // Boundary dimensions take what they measure on the result.
    for entry in measured {
        match measure(&solved, &entry.constraint) {
            Ok(v) if v.is_finite() && meaningful(&entry.constraint, v) => {
                if entry.constraint.value() != Some(v) {
                    values.push((entry.id, v));
                }
            }
            _ => {
                out.lost += 1;
                deletes.push(entry.id);
            }
        }
    }
    let mut edits = Vec::new();
    if !deletes.is_empty() {
        edits.push(SketchEdit::Delete { items: deletes });
    }
    let mut batch = Batch::default();
    for c in adds {
        batch.constraint(c);
    }
    edits.extend(batch.edits);
    edits.extend(
        values
            .into_iter()
            .map(|(constraint, value)| SketchEdit::SetValue { constraint, value }),
    );
    edits.extend(geometry_edits(sketch, &solved));
    edits.extend(turned_labels(labels, &kept_dims, &t, |id| id));
    // Preview: what moved, and what the solve moved with it.
    let changed: BTreeSet<SketchId> = sketch
        .entities
        .iter()
        .filter(
            |e| match (&e.geometry, solved.entity(e.id).map(|s| &s.geometry)) {
                (Geometry::Point { at }, Some(Geometry::Point { at: to })) => at != to,
                _ => false,
            },
        )
        .map(|e| e.id)
        .collect();
    let shown: Vec<ResolvedCurve> = resolve_curves(&solved)
        .unwrap_or_default()
        .into_iter()
        .filter(|(id, _, _)| {
            moved.contains(id)
                || solved.entity(*id).is_some_and(|e| {
                    e.geometry
                        .referenced_points()
                        .iter()
                        .any(|p| changed.contains(p))
                })
        })
        .map(|(_, _, c)| c)
        .collect();
    let points: Vec<DVec2> = lone_points(sketch, moved)
        .into_iter()
        .filter_map(|p| solved.point(p))
        .collect();
    out.plan = Plan {
        edits,
        removed: out.detached + out.lost + usize::from(out.freed),
        preview: Preview::Shapes {
            curves: shown,
            points,
        },
    };
    Ok(out)
}

/// Whether a dimension measured anew still measures something: lengths of
/// more than [`SAME_POINT`]; any angle.
fn meaningful(c: &Constraint<SketchId>, v: f64) -> bool {
    match c {
        Constraint::Angle { .. } => true,
        _ => v > SAME_POINT,
    }
}

/// Label offsets of the dimensions in `dims`, turned and scaled with them
/// (none for a translation), for the dimension ids `id` maps them to.
fn turned_labels(
    labels: &[(SketchId, DVec2)],
    dims: &BTreeSet<SketchId>,
    t: &Similarity,
    id: impl Fn(SketchId) -> SketchId,
) -> Vec<SketchEdit> {
    if t.quarter_turns() == Some(0) && t.scale == 1.0 {
        return Vec::new();
    }
    labels
        .iter()
        .filter(|(d, _)| dims.contains(d))
        .map(|&(d, offset)| SketchEdit::PlaceLabel {
            dimension: id(d),
            offset: Some(t.vector(offset)),
        })
        .collect()
}

fn copy(
    sketch: &SketchData,
    labels: &[(SketchId, DVec2)],
    curves: &[(SketchId, ResolvedCurve)],
    sources: &[SketchId],
    moved: &BTreeSet<SketchId>,
    decider: &Decider<'_>,
) -> Result<Transformed, String> {
    let t = decider.t;
    let mut batch = Batch::default();
    let mut map = BTreeMap::new();
    for e in &sketch.entities {
        if let (true, Geometry::Point { at }) = (moved.contains(&e.id), &e.geometry) {
            let id = batch.entity(Geometry::Point { at: t.point(*at) }, e.construction);
            map.insert(e.id, id);
        }
    }
    let mut preview_curves = Vec::new();
    for e in &sketch.entities {
        if !(moved.contains(&e.id) && e.geometry.is_curve()) {
            continue;
        }
        let mut g = e.geometry.map_keys(|k| map[&k]);
        if let Geometry::Circle { radius, .. } = &mut g {
            *radius *= t.scale;
        }
        let id = batch.entity(g, e.construction);
        map.insert(e.id, id);
        if let Some((_, c)) = curves.iter().find(|(c, _)| *c == e.id) {
            preview_curves.push(t.curve(c));
        }
    }
    let (mut lost, mut freed, mut converted) = (0, false, 0);
    let mut new_dims = BTreeMap::new();
    for entry in &sketch.constraints {
        match decider.copied(entry, &map) {
            Ok(Some(c)) => {
                if matches!(
                    (&entry.constraint, &c),
                    (
                        Constraint::Horizontal { .. } | Constraint::Vertical { .. },
                        Constraint::Parallel { .. } | Constraint::Perpendicular { .. }
                    )
                ) {
                    converted += 1;
                }
                let dim = c.value().is_some();
                let id = batch.constraint(c);
                if dim {
                    new_dims.insert(entry.id, id);
                }
            }
            Ok(None) => {}
            Err(Why::Freed) => freed = true,
            Err(_) => lost += 1,
        }
    }
    let dims: BTreeSet<SketchId> = new_dims.keys().copied().collect();
    let moved_labels: Vec<(SketchId, DVec2)> = labels
        .iter()
        .filter(|(d, _)| dims.contains(d))
        .map(|&(d, offset)| (d, offset))
        .collect();
    // Copies keep their labels' places relative to them, turned and
    // scaled (a translation keeps the offsets as they are).
    for (d, offset) in moved_labels {
        let offset = if t.quarter_turns() == Some(0) && t.scale == 1.0 {
            offset
        } else {
            t.vector(offset)
        };
        batch.edits.push(SketchEdit::PlaceLabel {
            dimension: new_dims[&d],
            offset: Some(offset),
        });
    }
    let points: Vec<DVec2> = lone_points(sketch, moved)
        .into_iter()
        .filter_map(|p| sketch.point(p).map(|at| t.point(at)))
        .collect();
    let result = sources.iter().map(|s| map[s]).collect();
    Ok(Transformed {
        plan: Plan {
            edits: batch.edits,
            removed: lost + usize::from(freed),
            preview: Preview::Shapes {
                curves: preview_curves,
                points,
            },
        },
        result,
        items: sources.len(),
        copy: true,
        detached: 0,
        lost,
        freed,
        converted,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn quarter_turns_are_exact() {
        let t = Similarity::rotation(DVec2::new(1.0, 2.0), 90f64.to_radians());
        assert_eq!(t.quarter_turns(), Some(1));
        assert_eq!(t.point(DVec2::new(4.0, 2.0)), DVec2::new(1.0, 5.0));
        let t = Similarity::rotation(DVec2::ZERO, -PI / 2.0);
        assert_eq!(t.quarter_turns(), Some(3));
        assert_eq!(t.point(DVec2::new(0.0, 3.0)), DVec2::new(3.0, 0.0));
        assert_eq!(
            Similarity::rotation(DVec2::ZERO, 2.0 * PI).quarter_turns(),
            Some(0)
        );
        assert_eq!(Similarity::rotation(DVec2::ZERO, 0.3).quarter_turns(), None);
        assert!(Similarity::rotation(DVec2::ZERO, 2.0 * PI).is_identity());
    }

    #[test]
    fn similarities_map_curves_consistently() {
        let t = Similarity {
            centre: DVec2::new(0.5, -1.0),
            angle: 0.7,
            scale: 2.5,
            shift: DVec2::new(0.1, 0.2),
        };
        let arc = ResolvedCurve::Arc {
            center: DVec2::new(1.0, 1.0),
            radius: 0.5,
            start_angle: 0.3,
            sweep_angle: 1.2,
        };
        let ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } = t.curve(&arc)
        else {
            unreachable!()
        };
        let start = DVec2::new(1.0, 1.0) + DVec2::from_angle(0.3) * 0.5;
        assert!(
            (center + DVec2::from_angle(start_angle) * radius - t.point(start)).length() < 1e-12
        );
        assert_eq!(sweep_angle, 1.2);
        assert!((radius - 1.25).abs() < 1e-15);
        // Distances scale by the factor.
        let (p, q) = (DVec2::new(3.0, -2.0), DVec2::new(-1.0, 4.0));
        let d = (t.point(p) - t.point(q)).length();
        assert!((d - 2.5 * (p - q).length()).abs() < 1e-12);
    }

    #[test]
    fn translations_are_exact() {
        let t = Similarity::translation(DVec2::new(0.025, -0.01));
        assert_eq!(
            t.point(DVec2::new(0.1, 0.2)),
            DVec2::new(0.1, 0.2) + t.shift
        );
        assert!(!t.is_identity());
        assert!(Similarity::translation(DVec2::ZERO).is_identity());
    }
}
