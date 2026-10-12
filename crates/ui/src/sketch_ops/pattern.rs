//! Patterns (design note 0007 §11): copies of curves and points along a
//! line or about a centre, held by one derived pattern constraint per use
//! of a pattern tool.
//!
//! As with mirrors, the copies get new points and no constraints of their
//! own: the pattern fixes every copy point from its source, so copies add
//! no degrees of freedom and a dimension on a copy would be redundant. A
//! linear pattern's step is a construction line the tool adds between the
//! user's two clicks (dimensioning it sets the spacing). A circular
//! pattern's copies are turned about a centre point by multiples of its
//! angle, a dimension; source points at the centre (the centre itself, or
//! coincident with it) are shared by every copy, so spokes drawn from a hub
//! stay joined there.

use std::collections::BTreeMap;
use std::f64::consts::TAU;

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{SketchData, SketchId};
use openshape_sketch_solver::ResolvedCurve;

use super::trim::{Plan, Preview};
use super::{coincidence_classes, Batch, SAME_POINT};
use crate::sketch_mode::{Snap, SnapTarget};

/// The fewest instances (the sources included) a pattern tool makes.
pub const MIN_COUNT: u32 = 2;
/// The most instances a pattern tool makes.
pub const MAX_COUNT: u32 = 100;

/// Where a pattern tool puts the copies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Anchor {
    /// Along a new construction line between two clicks: copy k is moved
    /// k times from `from` to `to`.
    Linear {
        /// The line's start.
        from: Snap,
        /// The line's end.
        to: Snap,
    },
    /// About a centre (the snapped point itself when it snaps to one, else
    /// a new point there): copy k is turned k times a whole turn divided
    /// by the count.
    Circular {
        /// The centre.
        centre: Snap,
    },
}

/// How a pattern moves its sources to copy k.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Motion {
    /// Moved k times by a step.
    Translate(DVec2),
    /// Turned k times by `angle` (radians, counter-clockwise) about
    /// `centre`.
    Rotate {
        /// The centre.
        centre: DVec2,
        /// One step's angle.
        angle: f64,
    },
}

impl Motion {
    /// Where copy `k` of `p` is.
    pub fn point(&self, p: DVec2, k: u32) -> DVec2 {
        match *self {
            Motion::Translate(v) => p + v * f64::from(k),
            Motion::Rotate { centre, angle } => {
                centre + DVec2::from_angle(angle * f64::from(k)).rotate(p - centre)
            }
        }
    }

    /// Copy `k` of a vector (a direction, not a position).
    fn vector(&self, v: DVec2, k: u32) -> DVec2 {
        match *self {
            Motion::Translate(_) => v,
            Motion::Rotate { angle, .. } => DVec2::from_angle(angle * f64::from(k)).rotate(v),
        }
    }

    /// Copy `k` of a solved curve (arcs keep their start and end;
    /// ellipses their eccentric anomalies, as motions keep handedness;
    /// conics their weight).
    pub fn curve(&self, c: &ResolvedCurve, k: u32) -> ResolvedCurve {
        let turn = match *self {
            Motion::Translate(_) => 0.0,
            Motion::Rotate { angle, .. } => angle * f64::from(k),
        };
        match c {
            ResolvedCurve::Line { start, end } => ResolvedCurve::Line {
                start: self.point(*start, k),
                end: self.point(*end, k),
            },
            ResolvedCurve::Circle { center, radius } => ResolvedCurve::Circle {
                center: self.point(*center, k),
                radius: *radius,
            },
            ResolvedCurve::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => ResolvedCurve::Arc {
                center: self.point(*center, k),
                radius: *radius,
                start_angle: start_angle + turn,
                sweep_angle: *sweep_angle,
            },
            ResolvedCurve::Ellipse {
                center,
                major,
                minor,
            } => ResolvedCurve::Ellipse {
                center: self.point(*center, k),
                major: self.vector(*major, k),
                minor: self.vector(*minor, k),
            },
            ResolvedCurve::EllipticalArc {
                center,
                major,
                minor,
                start_param,
                sweep_param,
            } => ResolvedCurve::EllipticalArc {
                center: self.point(*center, k),
                major: self.vector(*major, k),
                minor: self.vector(*minor, k),
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
                points: points.iter().map(|p| self.point(*p, k)).collect(),
            },
            ResolvedCurve::Conic {
                start,
                control,
                end,
                weight,
            } => ResolvedCurve::Conic {
                start: self.point(*start, k),
                control: self.point(*control, k),
                end: self.point(*end, k),
                weight: *weight,
            },
        }
    }
}

/// A pattern ready to commit.
#[derive(Clone, Debug, PartialEq)]
pub struct Patterned {
    /// The batch and its preview.
    pub plan: Plan,
    /// The pattern constraint (batch-local id).
    pub constraint: SketchId,
    /// The direction line of a linear pattern (batch-local id).
    pub line: Option<SketchId>,
    /// Sources patterned.
    pub made: usize,
    /// Sources skipped: on the centre of a circular pattern, so every copy
    /// would lie on them.
    pub skipped: usize,
    /// Instances, the sources included.
    pub count: u32,
}

impl Patterned {
    /// The status-bar message after committing.
    pub fn message(&self) -> String {
        let made = match self.made {
            1 => format!("Patterned 1 item into {} instances", self.count),
            n => format!("Patterned {n} items into {} instances", self.count),
        };
        match self.skipped {
            0 => made,
            n => format!("{made}; skipped {n} lying on the centre"),
        }
    }
}

fn point(sketch: &SketchData, id: SketchId) -> Result<DVec2, String> {
    sketch
        .point(id)
        .ok_or_else(|| format!("point {id:?} does not exist"))
}

/// Patterns `sources` (curves and points, [`super::mirror::sources`]) into
/// `count` instances: one batch adding the copies, the anchor's new
/// entities (a linear pattern's construction line, a circular pattern's
/// centre when it snapped to no point) and one pattern constraint.
///
/// `curves` are the sketch's solved curves, for the preview. Sources that
/// a circular pattern would turn onto themselves (a point at the centre, a
/// circle about it) are skipped.
pub fn pattern(
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    sources: &[SketchId],
    anchor: Anchor,
    count: u32,
) -> Result<Patterned, String> {
    if !(MIN_COUNT..=MAX_COUNT).contains(&count) {
        return Err(format!(
            "A pattern has {MIN_COUNT} to {MAX_COUNT} instances, not {count}"
        ));
    }
    if sources.is_empty() {
        return Err("Nothing to pattern: select curves or points".into());
    }
    let mut batch = Batch::default();
    let classes = coincidence_classes(sketch);
    let (motion, line, centre) = match anchor {
        Anchor::Linear { from, to } => {
            let v = to.at - from.at;
            if !(v.length().is_finite() && v.length() > SAME_POINT) {
                return Err("The step between copies has no length".into());
            }
            let start = batch.snapped_point(&from);
            let end = batch.snapped_point(&to);
            let line = batch.entity(Geometry::Line { start, end }, true);
            (Motion::Translate(v), Some(line), None)
        }
        Anchor::Circular { centre } => {
            let c = match centre.target {
                SnapTarget::Point(q) => q,
                _ => batch.snapped_point(&centre),
            };
            let angle = TAU / f64::from(count);
            (
                Motion::Rotate {
                    centre: centre.at,
                    angle,
                },
                None,
                Some((c, centre.at)),
            )
        }
    };
    // Source points every copy shares: the centre and points coincident
    // with it.
    let shared = |p: SketchId| {
        centre.is_some_and(|(c, _)| {
            p == c || classes.get(&p).is_some_and(|k| classes.get(&c) == Some(k))
        })
    };
    let at_centre = |p: DVec2| centre.is_some_and(|(_, at)| (p - at).length() <= SAME_POINT);
    // Copy k's point for each source point, k = 1..count.
    let mut images: Vec<BTreeMap<SketchId, SketchId>> = vec![BTreeMap::new(); count as usize];
    let mut copies = Vec::new();
    let mut preview_curves = Vec::new();
    let mut preview_points = Vec::new();
    let (mut made, mut skipped) = (0, 0);
    for &src in sources {
        let Some(e) = sketch.entity(src) else {
            continue;
        };
        let own_image = match &e.geometry {
            Geometry::Point { at } => shared(src) || at_centre(*at),
            Geometry::Circle { center, .. } => {
                shared(*center) || at_centre(point(sketch, *center)?)
            }
            _ => false,
        };
        if own_image {
            skipped += 1;
            continue;
        }
        let solved = curves.iter().find(|(id, _)| *id == src).map(|(_, c)| c);
        for k in 1..count {
            let map = &mut images[k as usize];
            for p in e.geometry.referenced_points() {
                if !shared(p) && !map.contains_key(&p) {
                    let q = batch.point(motion.point(point(sketch, p)?, k));
                    map.insert(p, q);
                }
            }
            let image = match &e.geometry {
                Geometry::Point { at } => {
                    preview_points.push(motion.point(*at, k));
                    Geometry::Point {
                        at: motion.point(*at, k),
                    }
                }
                g => g.map_keys(|p| map.get(&p).copied().unwrap_or(p)),
            };
            if let Some(c) = solved {
                preview_curves.push(motion.curve(c, k));
            }
            let id = batch.entity(image, e.construction);
            copies.push((src, id, k));
        }
        made += 1;
    }
    if made == 0 {
        return Err("The selection lies on the centre: turning it leaves it where it is".into());
    }
    let constraint = batch.constraint(match (motion, line, centre) {
        (Motion::Rotate { angle, .. }, _, Some((c, _))) => Constraint::CircularPattern {
            centre: c,
            angle,
            copies,
        },
        (_, Some(direction), _) => Constraint::LinearPattern { direction, copies },
        _ => unreachable!("a linear pattern has a line and a circular one a centre"),
    });
    Ok(Patterned {
        plan: Plan {
            edits: batch.edits,
            removed: 0,
            preview: Preview::Shapes {
                curves: preview_curves,
                points: preview_points,
            },
        },
        constraint,
        line,
        made,
        skipped,
        count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn motions_move_points_and_keep_arcs_running_the_same_way() {
        let step = Motion::Translate(DVec2::new(2.0, -1.0));
        assert_eq!(step.point(DVec2::new(1.0, 1.0), 3), DVec2::new(7.0, -2.0));
        let turn = Motion::Rotate {
            centre: DVec2::new(1.0, 0.0),
            angle: TAU / 4.0,
        };
        let p = turn.point(DVec2::new(3.0, 0.0), 1);
        assert!((p - DVec2::new(1.0, 2.0)).length() < 1e-12);
        let arc = ResolvedCurve::Arc {
            center: DVec2::new(3.0, 0.0),
            radius: 1.0,
            start_angle: 0.25,
            sweep_angle: 1.0,
        };
        let ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } = turn.curve(&arc, 2)
        else {
            unreachable!()
        };
        assert!((center - DVec2::new(-1.0, 0.0)).length() < 1e-12);
        assert_eq!((radius, sweep_angle), (1.0, 1.0));
        // The copy's start is the copy of the source's start.
        let src_start = DVec2::new(3.0, 0.0) + DVec2::from_angle(0.25);
        let img_start = center + DVec2::from_angle(start_angle);
        assert!((img_start - turn.point(src_start, 2)).length() < 1e-12);
    }
}
