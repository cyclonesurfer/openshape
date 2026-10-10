//! Sketch data: entities, constraints and dimensions.
//!
//! These types are pure data shared with the document model
//! (`openshape-core-model` stores them with its own key type, ADR 0011).
//! They are generic over the key `K` that identifies entities and
//! constraints. Coordinates are in the sketch plane, in metres; angles are in
//! radians, counter-clockwise positive.

use std::fmt::Debug;
use std::hash::Hash;

use glam::DVec2;
use serde::{Deserialize, Serialize};

/// Requirements on entity and constraint keys.
pub trait SketchKey: Copy + Eq + Ord + Hash + Debug {}
impl<T: Copy + Eq + Ord + Hash + Debug> SketchKey for T {}

/// The geometry of a sketch entity. Curves reference point entities by key;
/// connections between curves are expressed with constraints
/// (e.g. [`Constraint::Coincident`] between end points).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Geometry<K> {
    /// A point.
    Point {
        /// Position.
        at: DVec2,
    },
    /// A line segment from `start` to `end` (point entities).
    Line {
        /// Start point entity.
        start: K,
        /// End point entity.
        end: K,
    },
    /// A full circle about a point entity.
    Circle {
        /// Centre point entity.
        center: K,
        /// Radius (metres, > 0).
        radius: f64,
    },
    /// A circular arc about `center`, running counter-clockwise from `start`
    /// to `end`. Its radius is |start − center|; the solver keeps
    /// |end − center| equal to it.
    Arc {
        /// Centre point entity.
        center: K,
        /// Start point entity.
        start: K,
        /// End point entity.
        end: K,
    },
    /// A clamped B-spline with uniform interior knots on [0, 1], defined by
    /// its control points (point entities, at least `degree + 1`).
    ControlSpline {
        /// Control point entities, in order.
        control: Vec<K>,
        /// Degree, 1 to 5.
        degree: u32,
    },
    /// A cubic B-spline interpolating `through` (at least 2 point entities)
    /// at the fixed, strictly increasing parameters `params` (first 0, last
    /// 1). Its control points are derived (P&T A9.1, knots by averaging), so
    /// they add no degrees of freedom.
    InterpolatedSpline {
        /// Points the curve passes through, in order.
        through: Vec<K>,
        /// Curve parameter at each through-point.
        params: Vec<f64>,
    },
}

impl<K: SketchKey> Geometry<K> {
    /// The point entities this geometry references, in a fixed order.
    pub fn referenced_points(&self) -> Vec<K> {
        match self {
            Geometry::Point { .. } => Vec::new(),
            Geometry::Line { start, end } => vec![*start, *end],
            Geometry::Circle { center, .. } => vec![*center],
            Geometry::Arc { center, start, end } => vec![*center, *start, *end],
            Geometry::ControlSpline { control, .. } => control.clone(),
            Geometry::InterpolatedSpline { through, .. } => through.clone(),
        }
    }

    /// True for curves (everything but points).
    pub fn is_curve(&self) -> bool {
        !matches!(self, Geometry::Point { .. })
    }
}

/// A sketch entity.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Entity<K> {
    /// Key.
    pub id: K,
    /// Geometry.
    pub geometry: Geometry<K>,
    /// Construction geometry is constrained like any other but never bounds
    /// a region.
    pub construction: bool,
}

/// Which way two curves touch, for [`Constraint::Tangent`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TangentKind<K> {
    /// The curves share the end point `point` (their ends are coincident) and
    /// have parallel tangents there (G1). Works for any pair of curves.
    AtPoint {
        /// The shared point entity (an end point of both curves, or one of
        /// two coincident end points).
        point: K,
    },
    /// A line and a circle or arc touch; the centre stays on the recorded
    /// side of the line (left when looking from start to end).
    LineCircle {
        /// True if the centre is on the left of the line.
        center_left: bool,
    },
    /// Two circles or arcs touch, from outside or inside.
    CircleCircle {
        /// True for internal tangency (one inside the other).
        internal: bool,
    },
}

/// How a distance is measured, for [`Constraint::Distance`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum DistanceKind {
    /// True (shortest) distance.
    Aligned,
    /// Along the sketch x axis.
    Horizontal,
    /// Along the sketch y axis.
    Vertical,
}

/// A geometric constraint or a dimension. "Line" means a
/// [`Geometry::Line`] entity, "circular" a circle or arc, "point" a point
/// entity (including curve end points and centres).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Constraint<K> {
    /// Two points coincide.
    Coincident {
        /// First point.
        a: K,
        /// Second point.
        b: K,
    },
    /// A point lies on a curve (on a line's infinite extension, on a
    /// circular's full circle, on a spline within its domain).
    PointOnCurve {
        /// The point.
        point: K,
        /// The curve.
        curve: K,
    },
    /// A line is parallel to the sketch x axis.
    Horizontal {
        /// The line.
        line: K,
    },
    /// A line is parallel to the sketch y axis.
    Vertical {
        /// The line.
        line: K,
    },
    /// Two points have the same y.
    HorizontalPoints {
        /// First point.
        a: K,
        /// Second point.
        b: K,
    },
    /// Two points have the same x.
    VerticalPoints {
        /// First point.
        a: K,
        /// Second point.
        b: K,
    },
    /// Two lines are parallel.
    Parallel {
        /// First line.
        a: K,
        /// Second line.
        b: K,
    },
    /// Two lines are perpendicular.
    Perpendicular {
        /// First line.
        a: K,
        /// Second line.
        b: K,
    },
    /// Two lines have equal length, or two circulars equal radius.
    Equal {
        /// First curve.
        a: K,
        /// Second curve.
        b: K,
    },
    /// A point is the midpoint of a line.
    Midpoint {
        /// The point.
        point: K,
        /// The line.
        line: K,
    },
    /// Two circulars share a centre.
    Concentric {
        /// First circular.
        a: K,
        /// Second circular.
        b: K,
    },
    /// Two points are mirror images about a line.
    Symmetric {
        /// First point.
        a: K,
        /// Second point.
        b: K,
        /// The mirror line.
        axis: K,
    },
    /// Two curves are tangent.
    Tangent {
        /// First curve.
        a: K,
        /// Second curve.
        b: K,
        /// How they touch.
        kind: TangentKind<K>,
    },
    /// A line is normal to a circular (passes through its centre).
    Normal {
        /// The line.
        line: K,
        /// The circular.
        curve: K,
    },
    /// A point is fixed at a position.
    Fix {
        /// The point.
        point: K,
        /// Its fixed position.
        at: DVec2,
    },
    /// Dimension: distance between two points, a point and a line, or two
    /// parallel lines (metres, > 0).
    Distance {
        /// First point or line.
        a: K,
        /// Second point or line.
        b: K,
        /// How the distance is measured (only points support horizontal and
        /// vertical).
        kind: DistanceKind,
        /// Value.
        value: f64,
    },
    /// Dimension: length of a line (metres, > 0).
    Length {
        /// The line.
        line: K,
        /// Value.
        value: f64,
    },
    /// Dimension: signed angle from line `a`'s direction (start to end) to
    /// line `b`'s, in (−π, π] radians.
    Angle {
        /// First line.
        a: K,
        /// Second line.
        b: K,
        /// Value.
        value: f64,
    },
    /// Dimension: radius of a circular (metres, > 0).
    Radius {
        /// The circular.
        curve: K,
        /// Value.
        value: f64,
    },
    /// Dimension: diameter of a circular (metres, > 0).
    Diameter {
        /// The circular.
        curve: K,
        /// Value.
        value: f64,
    },
}

impl<K: SketchKey> Constraint<K> {
    /// The dimension value, for dimensions.
    pub fn value(&self) -> Option<f64> {
        match self {
            Constraint::Distance { value, .. }
            | Constraint::Length { value, .. }
            | Constraint::Angle { value, .. }
            | Constraint::Radius { value, .. }
            | Constraint::Diameter { value, .. } => Some(*value),
            _ => None,
        }
    }

    /// Sets the dimension value; returns false for non-dimensions.
    pub fn set_value(&mut self, v: f64) -> bool {
        match self {
            Constraint::Distance { value, .. }
            | Constraint::Length { value, .. }
            | Constraint::Angle { value, .. }
            | Constraint::Radius { value, .. }
            | Constraint::Diameter { value, .. } => {
                *value = v;
                true
            }
            _ => false,
        }
    }

    /// Every entity key this constraint references, in a fixed order.
    pub fn referenced(&self) -> Vec<K> {
        match self {
            Constraint::Coincident { a, b }
            | Constraint::HorizontalPoints { a, b }
            | Constraint::VerticalPoints { a, b }
            | Constraint::Parallel { a, b }
            | Constraint::Perpendicular { a, b }
            | Constraint::Equal { a, b }
            | Constraint::Concentric { a, b }
            | Constraint::Distance { a, b, .. }
            | Constraint::Angle { a, b, .. } => vec![*a, *b],
            Constraint::PointOnCurve { point, curve } => vec![*point, *curve],
            Constraint::Horizontal { line }
            | Constraint::Vertical { line }
            | Constraint::Length { line, .. } => vec![*line],
            Constraint::Midpoint { point, line } => vec![*point, *line],
            Constraint::Symmetric { a, b, axis } => vec![*a, *b, *axis],
            Constraint::Tangent { a, b, kind } => match kind {
                TangentKind::AtPoint { point } => vec![*a, *b, *point],
                _ => vec![*a, *b],
            },
            Constraint::Normal { line, curve } => vec![*line, *curve],
            Constraint::Fix { point, .. } => vec![*point],
            Constraint::Radius { curve, .. } | Constraint::Diameter { curve, .. } => vec![*curve],
        }
    }
}

/// A constraint with its key.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConstraintEntry<K> {
    /// Key.
    pub id: K,
    /// The constraint.
    pub constraint: Constraint<K>,
}

/// A sketch: entities and constraints, in insertion order. Order matters only
/// for diagnostics (earlier constraints win when deciding which of a
/// redundant set to report).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sketch<K> {
    /// Entities.
    pub entities: Vec<Entity<K>>,
    /// Constraints and dimensions.
    pub constraints: Vec<ConstraintEntry<K>>,
}

impl<K> Default for Sketch<K> {
    fn default() -> Self {
        Sketch {
            entities: Vec::new(),
            constraints: Vec::new(),
        }
    }
}

impl<K: SketchKey> Sketch<K> {
    /// The entity with key `id`.
    pub fn entity(&self, id: K) -> Option<&Entity<K>> {
        self.entities.iter().find(|e| e.id == id)
    }

    /// The position of point entity `id`.
    pub fn point(&self, id: K) -> Option<DVec2> {
        match self.entity(id)?.geometry {
            Geometry::Point { at } => Some(at),
            _ => None,
        }
    }

    /// The constraint with key `id`.
    pub fn constraint(&self, id: K) -> Option<&ConstraintEntry<K>> {
        self.constraints.iter().find(|c| c.id == id)
    }
}
