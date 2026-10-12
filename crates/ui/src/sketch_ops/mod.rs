//! Sketch editing operations (design note 0007): pure functions from a
//! sketch and the user's input to one [`SketchEdit`] batch. The editor in
//! [`crate::sketch_mode`] decides *what* the user asked for; these functions
//! compute the exact geometry and the constraints that keep it.

pub mod shapes;

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, ConstraintEntry, Entity, Geometry, TangentKind};
use openshape_core_model::{SketchData, SketchEdit, SketchId};
use openshape_sketch_solver::ResolvedCurve;

use crate::sketch_mode::{Snap, SnapTarget};

/// Positions closer than this are the same point (metres).
pub(crate) const SAME_POINT: f64 = 1e-7;

/// Builds one sketch edit batch with batch-local ids.
#[derive(Default)]
pub struct Batch {
    /// The edits so far.
    pub edits: Vec<SketchEdit>,
    next: u32,
}

impl Batch {
    /// Adds an entity.
    pub fn entity(&mut self, geometry: Geometry<SketchId>, construction: bool) -> SketchId {
        let id = SketchId::new_in_batch(self.next);
        self.next += 1;
        self.edits.push(SketchEdit::AddEntity(Entity {
            id,
            geometry,
            construction,
        }));
        id
    }

    /// Adds a point.
    pub fn point(&mut self, at: DVec2) -> SketchId {
        self.entity(Geometry::Point { at }, false)
    }

    /// Adds a line with two new end points.
    pub fn line(&mut self, a: DVec2, b: DVec2, construction: bool) -> Ends {
        let start = self.point(a);
        let end = self.point(b);
        let curve = self.entity(Geometry::Line { start, end }, construction);
        Ends { curve, start, end }
    }

    /// Adds a counter-clockwise arc about `center` from `a` to `b` with new
    /// points (`a` and `b` should be equally far from `center`).
    pub fn arc(&mut self, center: DVec2, a: DVec2, b: DVec2, construction: bool) -> Arc {
        let c = self.point(center);
        let start = self.point(a);
        let end = self.point(b);
        let curve = self.entity(
            Geometry::Arc {
                center: c,
                start,
                end,
            },
            construction,
        );
        Arc {
            curve,
            center: c,
            start,
            end,
        }
    }

    /// Adds a constraint.
    pub fn constraint(&mut self, constraint: Constraint<SketchId>) -> SketchId {
        let id = SketchId::new_in_batch(self.next);
        self.next += 1;
        self.edits.push(SketchEdit::AddConstraint(ConstraintEntry {
            id,
            constraint,
        }));
        id
    }

    /// Makes two points coincident.
    pub fn coincident(&mut self, a: SketchId, b: SketchId) -> SketchId {
        self.constraint(Constraint::Coincident { a, b })
    }

    /// Makes two curves tangent where they meet at `point` (an end point of
    /// both, or coincident with one of each).
    pub fn tangent_at(&mut self, a: SketchId, b: SketchId, point: SketchId) -> SketchId {
        self.constraint(Constraint::Tangent {
            a,
            b,
            kind: TangentKind::AtPoint { point },
        })
    }

    /// A new point at a snap, constrained to what it snapped to.
    pub fn snapped_point(&mut self, s: &Snap) -> SketchId {
        let p = self.point(s.at);
        self.attach(p, s);
        p
    }

    /// Constrains an existing new point to what `s` snapped to.
    pub fn attach(&mut self, p: SketchId, s: &Snap) {
        match s.target {
            SnapTarget::Free => {}
            SnapTarget::Point(q) => {
                self.coincident(q, p);
            }
            SnapTarget::Curve(c) => {
                self.constraint(Constraint::PointOnCurve { point: p, curve: c });
            }
            SnapTarget::Origin => {
                self.constraint(Constraint::Fix {
                    point: p,
                    at: DVec2::ZERO,
                });
            }
        }
    }
}

/// A line (or another curve with two ends) added to a [`Batch`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ends {
    /// The curve.
    pub curve: SketchId,
    /// Its start point.
    pub start: SketchId,
    /// Its end point.
    pub end: SketchId,
}

/// An arc added to a [`Batch`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Arc {
    /// The arc.
    pub curve: SketchId,
    /// Its centre point.
    pub center: SketchId,
    /// Its start point.
    pub start: SketchId,
    /// Its end point.
    pub end: SketchId,
}

/// The left normal of `v` (turned a quarter counter-clockwise).
pub fn perp(v: DVec2) -> DVec2 {
    DVec2::new(-v.y, v.x)
}

/// The curve that has `point` as its start or end point, with `true` when
/// it is the end: the curve a new tangent curve continues.
pub fn curve_ending_at(sketch: &SketchData, point: SketchId) -> Option<(SketchId, bool)> {
    sketch.entities.iter().find_map(|e| {
        let (first, last) = match &e.geometry {
            Geometry::Line { start, end } | Geometry::Arc { start, end, .. } => (*start, *end),
            Geometry::ControlSpline { control: p, .. }
            | Geometry::InterpolatedSpline { through: p, .. } => (*p.first()?, *p.last()?),
            _ => return None,
        };
        if last == point {
            Some((e.id, true))
        } else if first == point {
            Some((e.id, false))
        } else {
            None
        }
    })
}

/// The unit direction in which a curve leaves its end (`at_end`) or start
/// (`!at_end`): the direction a tangent continuation starts in.
pub fn leaving_direction(curve: &ResolvedCurve, at_end: bool) -> Option<DVec2> {
    let d = match curve {
        ResolvedCurve::Line { start, end } => {
            if at_end {
                *end - *start
            } else {
                *start - *end
            }
        }
        ResolvedCurve::Arc {
            start_angle,
            sweep_angle,
            ..
        } => {
            // Counter-clockwise tangent at an angle a is (−sin a, cos a).
            if at_end {
                perp(DVec2::from_angle(start_angle + sweep_angle))
            } else {
                -perp(DVec2::from_angle(*start_angle))
            }
        }
        ResolvedCurve::Spline { points, .. } => {
            // A clamped B-spline is tangent to its first and last legs.
            let n = points.len();
            if n < 2 {
                return None;
            }
            if at_end {
                points[n - 1] - points[n - 2]
            } else {
                points[0] - points[1]
            }
        }
        ResolvedCurve::Circle { .. } => return None,
    };
    let len = d.length();
    (len.is_finite() && len > 0.0).then(|| d / len)
}
