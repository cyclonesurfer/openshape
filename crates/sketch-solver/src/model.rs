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

/// Smallest accepted conic shape value ρ; the largest is 1 − `RHO_MIN`.
/// Beyond them a conic lies within a millionth of its triangle's size of
/// its chord or its control polygon (design note 0009 §4.6).
pub const RHO_MIN: f64 = 1e-6;

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
    ///
    /// When `periodic`, the spline is closed smoothly (ADR 0016): `control`
    /// ends on its first point again, and the curve is the uniform periodic
    /// B-spline on the closed control polygon (at least `degree + 1`
    /// distinct points, degree 2 to 5). It does not pass through its control
    /// points and has no end points.
    ControlSpline {
        /// Control point entities, in order.
        control: Vec<K>,
        /// Degree, 1 to 5.
        degree: u32,
        /// Closed smoothly on its first point (absent in documents written
        /// before ADR 0016, which read as `false`).
        #[serde(default, skip_serializing_if = "is_false")]
        periodic: bool,
    },
    /// A full ellipse `E(θ) = c + cos θ·A + sin θ·B` with semi-axes
    /// `A = major − center` and `B = minor − center` (design note 0008 §3).
    /// The solver keeps `A ⟂ B`; nothing forces `‖A‖ ≥ ‖B‖`.
    Ellipse {
        /// Centre point entity.
        center: K,
        /// Point entity at the end of the first (major) semi-axis.
        major: K,
        /// Point entity at the end of the second (minor) semi-axis.
        minor: K,
    },
    /// An arc of the ellipse on `center`, `major` and `minor` (as
    /// [`Geometry::Ellipse`]), running counter-clockwise in the sketch plane
    /// from `start` to `end`; the solver keeps both on the ellipse.
    EllipticalArc {
        /// Centre point entity.
        center: K,
        /// Point entity at the end of the first (major) semi-axis.
        major: K,
        /// Point entity at the end of the second (minor) semi-axis.
        minor: K,
        /// Start point entity.
        start: K,
        /// End point entity.
        end: K,
    },
    /// A conic arc (design note 0009 §3): the rational quadratic Bézier
    /// from `start` to `end` whose end tangents meet at `control`, with
    /// middle weight `rho / (1 − rho)`. `rho` is stored like a circle's
    /// radius: the solver treats it as an unknown and writes it back.
    Conic {
        /// Start point entity (on the curve).
        start: K,
        /// Point entity where the end tangents meet.
        control: K,
        /// End point entity (on the curve).
        end: K,
        /// Shape value in (0, 1): below ½ elliptic, ½ parabolic, above ½
        /// hyperbolic.
        rho: f64,
    },
    /// A cubic B-spline interpolating `through` (at least 2 point entities)
    /// at the fixed, strictly increasing parameters `params` (first 0, last
    /// 1). Its control points are derived (P&T A9.1, knots by averaging), so
    /// they add no degrees of freedom.
    ///
    /// When `periodic`, the spline is closed smoothly (ADR 0016): `through`
    /// ends on its first point again (at parameter 1), and the curve is the
    /// periodic cubic spline with knots at the parameters, twice
    /// continuously differentiable everywhere, the start included; at least
    /// three distinct points.
    InterpolatedSpline {
        /// Points the curve passes through, in order.
        through: Vec<K>,
        /// Curve parameter at each through-point.
        params: Vec<f64>,
        /// Closed smoothly on its first point (absent in documents written
        /// before ADR 0016, which read as `false`).
        #[serde(default, skip_serializing_if = "is_false")]
        periodic: bool,
    },
    /// A line of text (ADR 0015) filling the box on its four corner
    /// points. The solver keeps the corners a rectangle `aspect` times as
    /// wide as it is high, read the right way round (a mirror image when
    /// `mirrored`): with `R` the quarter turn clockwise (anticlockwise when
    /// mirrored), `bottom_right − bottom_left = aspect·R(top_left −
    /// bottom_left)` and `top_right − top_left = bottom_right −
    /// bottom_left`. The letters are `outline`, laid out once when the text
    /// was typed and stored, so they never change unless the text is
    /// edited.
    Text {
        /// Corner point entities in the text's own frame: bottom-left,
        /// bottom-right, top-right, top-left.
        corners: [K; 4],
        /// The characters.
        content: String,
        /// The bundled font they were laid out in.
        font: String,
        /// Box width over height (> 0).
        aspect: f64,
        /// True for a mirror image of the text.
        mirrored: bool,
        /// Closed contours in box units: a box point `(u, v)`, `u` in
        /// `[0, aspect]` and `v` in `[0, 1]`, lies at `bottom_left + (u /
        /// aspect)·(bottom_right − bottom_left) + v·(top_left −
        /// bottom_left)`.
        outline: Vec<Vec<TextPiece>>,
    },
}

/// The largest accepted text aspect ratio (box width over height): about a
/// thousand wide characters.
pub const MAX_TEXT_ASPECT: f64 = 2000.0;

/// Structural checks of a text's data that need no other entity (ADR
/// 0015): four distinct corners, an aspect ratio in (0,
/// [`MAX_TEXT_ASPECT`]], pieces of degree 1 (two points) to 3 (at least
/// four) with finite coordinates within reach of the box, none ending where
/// it starts, and contours of at least two pieces that join and close
/// exactly.
pub fn check_text<K: SketchKey>(
    corners: &[K; 4],
    aspect: f64,
    outline: &[Vec<TextPiece>],
) -> Result<(), String> {
    if (0..4).any(|a| (a + 1..4).any(|b| corners[a] == corners[b])) {
        return Err("a text box needs four distinct corner points".into());
    }
    if !(aspect.is_finite() && aspect > 0.0 && aspect <= MAX_TEXT_ASPECT) {
        return Err(format!(
            "a text's aspect ratio must lie in (0, {MAX_TEXT_ASPECT}] (got {aspect})"
        ));
    }
    let reach = 4.0 * MAX_TEXT_ASPECT;
    for (i, contour) in outline.iter().enumerate() {
        if contour.len() < 2 {
            return Err(format!(
                "text outline contour {i} has fewer than two pieces"
            ));
        }
        for (j, piece) in contour.iter().enumerate() {
            let n = piece.points.len();
            let ok = match piece.degree {
                1 => n == 2,
                2 | 3 => n > piece.degree as usize,
                _ => false,
            };
            if !ok {
                return Err(format!(
                    "text outline contour {i} piece {j}: degree {} with {n} points",
                    piece.degree
                ));
            }
            if piece
                .points
                .iter()
                .any(|p| !(p.is_finite() && p.x.abs() <= reach && p.y.abs() <= reach))
            {
                return Err(format!(
                    "text outline contour {i} piece {j}: a point is out of range"
                ));
            }
            if piece.points[0] == piece.points[n - 1] {
                return Err(format!(
                    "text outline contour {i} piece {j} ends where it starts"
                ));
            }
            let next = &contour[(j + 1) % contour.len()];
            if piece.points[n - 1] != next.points[0] {
                return Err(format!(
                    "text outline contour {i} piece {j} does not meet the next"
                ));
            }
        }
    }
    Ok(())
}

/// One piece of a text outline (ADR 0015), in box units: a line (degree 1,
/// two points) or a clamped B-spline of degree 2 or 3 with uniform interior
/// knots on [0, 1] (`points` are its control points). Each piece of a
/// contour starts exactly where the previous one ends, and the last ends
/// exactly where the first starts.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TextPiece {
    /// 1, 2 or 3.
    pub degree: u32,
    /// The line's ends or the spline's control points.
    pub points: Vec<DVec2>,
}

impl<K: SketchKey> Geometry<K> {
    /// The point entities this geometry references, in a fixed order.
    pub fn referenced_points(&self) -> Vec<K> {
        match self {
            Geometry::Point { .. } => Vec::new(),
            Geometry::Line { start, end } => vec![*start, *end],
            Geometry::Circle { center, .. } => vec![*center],
            Geometry::Arc { center, start, end } => vec![*center, *start, *end],
            Geometry::Ellipse {
                center,
                major,
                minor,
            } => vec![*center, *major, *minor],
            Geometry::EllipticalArc {
                center,
                major,
                minor,
                start,
                end,
            } => vec![*center, *major, *minor, *start, *end],
            Geometry::Conic {
                start,
                control,
                end,
                ..
            } => vec![*start, *control, *end],
            Geometry::ControlSpline { control, .. } => control.clone(),
            Geometry::InterpolatedSpline { through, .. } => through.clone(),
            Geometry::Text { corners, .. } => corners.to_vec(),
        }
    }

    /// The start and end points of an open curve: a line, an arc, an
    /// elliptical arc, a conic or a spline that is not smooth closed. `None`
    /// for points, circles, ellipses, smooth closed splines (ADR 0016) and
    /// text, which have no ends.
    pub fn end_points(&self) -> Option<(K, K)> {
        match self {
            Geometry::Line { start, end }
            | Geometry::Arc { start, end, .. }
            | Geometry::EllipticalArc { start, end, .. }
            | Geometry::Conic { start, end, .. } => Some((*start, *end)),
            Geometry::ControlSpline {
                control: p,
                periodic: false,
                ..
            }
            | Geometry::InterpolatedSpline {
                through: p,
                periodic: false,
                ..
            } => Some((*p.first()?, *p.last()?)),
            _ => None,
        }
    }

    /// True for a smooth closed spline (ADR 0016).
    pub fn is_closed_spline(&self) -> bool {
        matches!(
            self,
            Geometry::ControlSpline { periodic: true, .. }
                | Geometry::InterpolatedSpline { periodic: true, .. }
        )
    }

    /// True for curves (everything but points).
    pub fn is_curve(&self) -> bool {
        !matches!(self, Geometry::Point { .. })
    }

    /// The same geometry with every referenced key passed through `f`.
    pub fn map_keys<L>(&self, f: impl Fn(K) -> L) -> Geometry<L> {
        match self {
            Geometry::Point { at } => Geometry::Point { at: *at },
            Geometry::Line { start, end } => Geometry::Line {
                start: f(*start),
                end: f(*end),
            },
            Geometry::Circle { center, radius } => Geometry::Circle {
                center: f(*center),
                radius: *radius,
            },
            Geometry::Arc { center, start, end } => Geometry::Arc {
                center: f(*center),
                start: f(*start),
                end: f(*end),
            },
            Geometry::Ellipse {
                center,
                major,
                minor,
            } => Geometry::Ellipse {
                center: f(*center),
                major: f(*major),
                minor: f(*minor),
            },
            Geometry::EllipticalArc {
                center,
                major,
                minor,
                start,
                end,
            } => Geometry::EllipticalArc {
                center: f(*center),
                major: f(*major),
                minor: f(*minor),
                start: f(*start),
                end: f(*end),
            },
            Geometry::Conic {
                start,
                control,
                end,
                rho,
            } => Geometry::Conic {
                start: f(*start),
                control: f(*control),
                end: f(*end),
                rho: *rho,
            },
            Geometry::ControlSpline {
                control,
                degree,
                periodic,
            } => Geometry::ControlSpline {
                control: control.iter().map(|k| f(*k)).collect(),
                degree: *degree,
                periodic: *periodic,
            },
            Geometry::InterpolatedSpline {
                through,
                params,
                periodic,
            } => Geometry::InterpolatedSpline {
                through: through.iter().map(|k| f(*k)).collect(),
                params: params.clone(),
                periodic: *periodic,
            },
            Geometry::Text {
                corners,
                content,
                font,
                aspect,
                mirrored,
                outline,
            } => Geometry::Text {
                corners: corners.map(&f),
                content: content.clone(),
                font: font.clone(),
                aspect: *aspect,
                mirrored: *mirrored,
                outline: outline.clone(),
            },
        }
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
    /// A line and a circle, arc, ellipse or elliptical arc touch; the centre
    /// stays on the recorded side of the line (left when looking from start
    /// to end).
    LineCircle {
        /// True if the centre is on the left of the line.
        center_left: bool,
    },
    /// Two circles or arcs touch, from outside or inside (not ellipses yet:
    /// design note 0008 §4.5).
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
    /// circular's full circle, on an ellipse or elliptical arc's full
    /// ellipse, on a spline within its domain).
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
    /// Two curves meet at `point` with one tangent line and one curvature
    /// (G2, design note 0010): a spline (degree ≥ 2) and a line, arc,
    /// elliptical arc, conic arc or another such spline, sharing the end
    /// point `point` as [`TangentKind::AtPoint`] does. Implies tangency.
    Curvature {
        /// First curve.
        a: K,
        /// Second curve.
        b: K,
        /// The shared point entity (an end point of both curves, or one of
        /// two coincident end points).
        point: K,
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
    /// Dimension: a conic's shape value ρ (design note 0009 §4.4), a plain
    /// number in [`RHO_MIN`], 1 − [`RHO_MIN`]].
    Rho {
        /// The conic.
        curve: K,
        /// Value.
        value: f64,
    },
    /// Derived geometry (design note 0007 §9): each image entity is the
    /// mirror image of its source about the line `axis`, and moves with it.
    ///
    /// A pair is two points, lines, circles, arcs, ellipses or elliptical
    /// arcs, or two splines of the same form (control splines of one
    /// degree, or interpolated splines with the same parameters) with as
    /// many points. Each image point is the reflection of the corresponding
    /// source point, except that a reflection reverses the sense of
    /// rotation, so an image arc (circular or elliptical) starts at
    /// the reflection of its source's end and ends at the reflection of its
    /// start; an image circle has its source's radius. A source point on
    /// the axis may be shared with the image (the point is then held on the
    /// axis). Images add no degrees of freedom: a dimension on an image of a
    /// dimensioned curve is redundant. Deleting a source or an image removes
    /// its pair ([`Constraint::without`]).
    Mirror {
        /// The mirror line.
        axis: K,
        /// `(source, image)` entity pairs.
        pairs: Vec<(K, K)>,
    },
    /// Derived geometry and a dimension (design note 0007 §10): each image
    /// curve is its source offset by `distance` (metres, > 0) to one side,
    /// and moves with it.
    ///
    /// A piece `(source, image, left)` is two lines, two arcs or two
    /// circles. The image lies `distance` to the source's left (`true`) or
    /// right in the source's own direction, so towards an arc's or circle's
    /// centre when `left`; image arcs and circles share their source's
    /// centre point and keep its direction. Images that share an end point
    /// meet there where their offset curves cross, or, at the points in
    /// `smooth` (joints whose sources meet tangentially), on the joint's
    /// normal; an image end shared with no other piece lies on the normal
    /// at its source's end. Images add no degrees of freedom.
    ///
    /// After every solve the offset must still be valid: every image end on
    /// its offset curve (at a smooth joint that holds only while the
    /// sources stay tangent), no line image reversed and no arc image
    /// turned back, and no two images crossing or touching except
    /// consecutive ones at their joint. Otherwise the solve fails with an
    /// error naming the offset. Deleting a source or an image removes its
    /// piece ([`Constraint::without`]).
    Offset {
        /// Distance (metres, > 0).
        distance: f64,
        /// `(source, image, left)` pieces.
        pieces: Vec<(K, K, bool)>,
        /// Image points at smooth joints.
        smooth: Vec<K>,
    },
    /// Derived geometry (design note 0007 §11): each image entity is its
    /// source moved by `k` times the vector of the line `direction` (from
    /// its start to its end), and moves with it.
    ///
    /// A copy `(source, image, k)`, k ≥ 1, is two points, lines, circles,
    /// arcs, ellipses or elliptical arcs, or two splines of the same form
    /// with as many points; each
    /// image point is its source point moved, an image arc keeps its
    /// source's start and end, an image circle its radius. Images add no
    /// degrees of freedom. Deleting a source or an image removes its copy;
    /// deleting the line removes the constraint.
    LinearPattern {
        /// The line whose vector is one step.
        direction: K,
        /// `(source, image, k)` copies.
        copies: Vec<(K, K, u32)>,
    },
    /// Derived geometry and a dimension (design note 0007 §11): each image
    /// entity is its source turned about the point `centre` by `k` times
    /// `angle` (radians, counter-clockwise, 0 < |angle| < 2π), and moves
    /// with it.
    ///
    /// Copies are as for [`Constraint::LinearPattern`]; a source point that
    /// is the centre (or coincident with it) may be its own image. No copy
    /// may be turned by a whole number of turns (it would lie on its
    /// source). Deleting the centre removes the constraint.
    CircularPattern {
        /// The centre point.
        centre: K,
        /// The angle between consecutive copies (radians).
        angle: f64,
        /// `(source, image, k)` copies.
        copies: Vec<(K, K, u32)>,
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
            | Constraint::Diameter { value, .. }
            | Constraint::Rho { value, .. } => Some(*value),
            Constraint::Offset { distance, .. } => Some(*distance),
            Constraint::CircularPattern { angle, .. } => Some(*angle),
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
            | Constraint::Diameter { value, .. }
            | Constraint::Rho { value, .. }
            | Constraint::Offset {
                distance: value, ..
            }
            | Constraint::CircularPattern { angle: value, .. } => {
                *value = v;
                true
            }
            _ => false,
        }
    }

    /// The same constraint with every referenced key passed through `f`
    /// (values and flags are kept).
    pub fn map_keys<L>(&self, f: impl Fn(K) -> L) -> Constraint<L> {
        use Constraint as C;
        match self {
            C::Coincident { a, b } => C::Coincident { a: f(*a), b: f(*b) },
            C::PointOnCurve { point, curve } => C::PointOnCurve {
                point: f(*point),
                curve: f(*curve),
            },
            C::Horizontal { line } => C::Horizontal { line: f(*line) },
            C::Vertical { line } => C::Vertical { line: f(*line) },
            C::HorizontalPoints { a, b } => C::HorizontalPoints { a: f(*a), b: f(*b) },
            C::VerticalPoints { a, b } => C::VerticalPoints { a: f(*a), b: f(*b) },
            C::Parallel { a, b } => C::Parallel { a: f(*a), b: f(*b) },
            C::Perpendicular { a, b } => C::Perpendicular { a: f(*a), b: f(*b) },
            C::Equal { a, b } => C::Equal { a: f(*a), b: f(*b) },
            C::Midpoint { point, line } => C::Midpoint {
                point: f(*point),
                line: f(*line),
            },
            C::Concentric { a, b } => C::Concentric { a: f(*a), b: f(*b) },
            C::Symmetric { a, b, axis } => C::Symmetric {
                a: f(*a),
                b: f(*b),
                axis: f(*axis),
            },
            C::Tangent { a, b, kind } => C::Tangent {
                a: f(*a),
                b: f(*b),
                kind: match kind {
                    TangentKind::AtPoint { point } => TangentKind::AtPoint { point: f(*point) },
                    TangentKind::LineCircle { center_left } => TangentKind::LineCircle {
                        center_left: *center_left,
                    },
                    TangentKind::CircleCircle { internal } => TangentKind::CircleCircle {
                        internal: *internal,
                    },
                },
            },
            C::Curvature { a, b, point } => C::Curvature {
                a: f(*a),
                b: f(*b),
                point: f(*point),
            },
            C::Normal { line, curve } => C::Normal {
                line: f(*line),
                curve: f(*curve),
            },
            C::Fix { point, at } => C::Fix {
                point: f(*point),
                at: *at,
            },
            C::Distance { a, b, kind, value } => C::Distance {
                a: f(*a),
                b: f(*b),
                kind: *kind,
                value: *value,
            },
            C::Length { line, value } => C::Length {
                line: f(*line),
                value: *value,
            },
            C::Angle { a, b, value } => C::Angle {
                a: f(*a),
                b: f(*b),
                value: *value,
            },
            C::Radius { curve, value } => C::Radius {
                curve: f(*curve),
                value: *value,
            },
            C::Diameter { curve, value } => C::Diameter {
                curve: f(*curve),
                value: *value,
            },
            C::Rho { curve, value } => C::Rho {
                curve: f(*curve),
                value: *value,
            },
            C::Mirror { axis, pairs } => C::Mirror {
                axis: f(*axis),
                pairs: pairs.iter().map(|(s, i)| (f(*s), f(*i))).collect(),
            },
            C::Offset {
                distance,
                pieces,
                smooth,
            } => C::Offset {
                distance: *distance,
                pieces: pieces.iter().map(|(s, i, l)| (f(*s), f(*i), *l)).collect(),
                smooth: smooth.iter().map(|p| f(*p)).collect(),
            },
            C::LinearPattern { direction, copies } => C::LinearPattern {
                direction: f(*direction),
                copies: copies.iter().map(|(s, i, k)| (f(*s), f(*i), *k)).collect(),
            },
            C::CircularPattern {
                centre,
                angle,
                copies,
            } => C::CircularPattern {
                centre: f(*centre),
                angle: *angle,
                copies: copies.iter().map(|(s, i, k)| (f(*s), f(*i), *k)).collect(),
            },
        }
    }

    /// The constraint left after deleting the entities for which `gone` is
    /// true, or `None` when nothing of it remains. A derived constraint
    /// ([`Constraint::Mirror`], [`Constraint::Offset`] and the patterns)
    /// loses the pairs, pieces or copies that reference a deleted entity
    /// (an offset also its deleted smooth points), and goes with a mirror's
    /// axis, a pattern's line or centre, or its last pair, piece or copy;
    /// any other constraint goes when an entity it references does.
    pub fn without(&self, gone: impl Fn(K) -> bool) -> Option<Constraint<K>> {
        match self {
            Constraint::Mirror { axis, pairs } => {
                if gone(*axis) {
                    return None;
                }
                let kept: Vec<(K, K)> = pairs
                    .iter()
                    .copied()
                    .filter(|(s, i)| !gone(*s) && !gone(*i))
                    .collect();
                (!kept.is_empty()).then_some(Constraint::Mirror {
                    axis: *axis,
                    pairs: kept,
                })
            }
            Constraint::Offset {
                distance,
                pieces,
                smooth,
            } => {
                let kept: Vec<(K, K, bool)> = pieces
                    .iter()
                    .copied()
                    .filter(|(s, i, _)| !gone(*s) && !gone(*i))
                    .collect();
                let smooth: Vec<K> = smooth.iter().copied().filter(|p| !gone(*p)).collect();
                (!kept.is_empty()).then_some(Constraint::Offset {
                    distance: *distance,
                    pieces: kept,
                    smooth,
                })
            }
            Constraint::LinearPattern { direction, copies } => {
                if gone(*direction) {
                    return None;
                }
                let kept = kept_copies(copies, &gone);
                (!kept.is_empty()).then_some(Constraint::LinearPattern {
                    direction: *direction,
                    copies: kept,
                })
            }
            Constraint::CircularPattern {
                centre,
                angle,
                copies,
            } => {
                if gone(*centre) {
                    return None;
                }
                let kept = kept_copies(copies, &gone);
                (!kept.is_empty()).then_some(Constraint::CircularPattern {
                    centre: *centre,
                    angle: *angle,
                    copies: kept,
                })
            }
            _ => (!self.referenced().into_iter().any(gone)).then(|| self.clone()),
        }
    }

    /// For a constraint holding two curves at a shared end point (a
    /// tangency at a point, curvature continuity): `(a, b, point)`.
    pub fn joint(&self) -> Option<(K, K, K)> {
        match *self {
            Constraint::Tangent {
                a,
                b,
                kind: TangentKind::AtPoint { point },
            }
            | Constraint::Curvature { a, b, point } => Some((a, b, point)),
            _ => None,
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
            Constraint::Curvature { a, b, point } => vec![*a, *b, *point],
            Constraint::Normal { line, curve } => vec![*line, *curve],
            Constraint::Fix { point, .. } => vec![*point],
            Constraint::Radius { curve, .. }
            | Constraint::Diameter { curve, .. }
            | Constraint::Rho { curve, .. } => vec![*curve],
            Constraint::Mirror { axis, pairs } => std::iter::once(*axis)
                .chain(pairs.iter().flat_map(|(s, i)| [*s, *i]))
                .collect(),
            Constraint::Offset { pieces, smooth, .. } => pieces
                .iter()
                .flat_map(|(s, i, _)| [*s, *i])
                .chain(smooth.iter().copied())
                .collect(),
            Constraint::LinearPattern {
                direction: first,
                copies,
            }
            | Constraint::CircularPattern {
                centre: first,
                copies,
                ..
            } => std::iter::once(*first)
                .chain(copies.iter().flat_map(|(s, i, _)| [*s, *i]))
                .collect(),
        }
    }
}

/// Serde's test for leaving a `false` flag out of a document.
fn is_false(b: &bool) -> bool {
    !*b
}

/// The copies of a pattern that reference no deleted entity.
fn kept_copies<K: SketchKey>(
    copies: &[(K, K, u32)],
    gone: &impl Fn(K) -> bool,
) -> Vec<(K, K, u32)> {
    copies
        .iter()
        .copied()
        .filter(|(s, i, _)| !gone(*s) && !gone(*i))
        .collect()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn map_keys_maps_exactly_the_referenced_keys() {
        let all: Vec<Constraint<u32>> = vec![
            Constraint::Coincident { a: 1, b: 2 },
            Constraint::PointOnCurve { point: 1, curve: 2 },
            Constraint::Horizontal { line: 1 },
            Constraint::Vertical { line: 1 },
            Constraint::HorizontalPoints { a: 1, b: 2 },
            Constraint::VerticalPoints { a: 1, b: 2 },
            Constraint::Parallel { a: 1, b: 2 },
            Constraint::Perpendicular { a: 1, b: 2 },
            Constraint::Equal { a: 1, b: 2 },
            Constraint::Midpoint { point: 1, line: 2 },
            Constraint::Concentric { a: 1, b: 2 },
            Constraint::Symmetric {
                a: 1,
                b: 2,
                axis: 3,
            },
            Constraint::Tangent {
                a: 1,
                b: 2,
                kind: TangentKind::AtPoint { point: 3 },
            },
            Constraint::Tangent {
                a: 1,
                b: 2,
                kind: TangentKind::LineCircle { center_left: true },
            },
            Constraint::Tangent {
                a: 1,
                b: 2,
                kind: TangentKind::CircleCircle { internal: true },
            },
            Constraint::Curvature {
                a: 1,
                b: 2,
                point: 3,
            },
            Constraint::Normal { line: 1, curve: 2 },
            Constraint::Fix {
                point: 1,
                at: DVec2::new(1.0, 2.0),
            },
            Constraint::Distance {
                a: 1,
                b: 2,
                kind: DistanceKind::Vertical,
                value: 0.5,
            },
            Constraint::Length {
                line: 1,
                value: 0.5,
            },
            Constraint::Angle {
                a: 1,
                b: 2,
                value: 0.5,
            },
            Constraint::Radius {
                curve: 1,
                value: 0.5,
            },
            Constraint::Diameter {
                curve: 1,
                value: 0.5,
            },
            Constraint::Rho {
                curve: 1,
                value: 0.5,
            },
            Constraint::Mirror {
                axis: 1,
                pairs: vec![(2, 3), (4, 5)],
            },
            Constraint::Offset {
                distance: 0.5,
                pieces: vec![(2, 3, true), (4, 5, false)],
                smooth: vec![6],
            },
            Constraint::LinearPattern {
                direction: 1,
                copies: vec![(2, 3, 1), (2, 4, 2)],
            },
            Constraint::CircularPattern {
                centre: 1,
                angle: 0.5,
                copies: vec![(2, 3, 1), (4, 5, 1)],
            },
        ];
        for c in &all {
            let m = c.map_keys(|k| k + 10);
            let expected: Vec<u32> = c.referenced().iter().map(|k| k + 10).collect();
            assert_eq!(m.referenced(), expected, "{c:?}");
            assert_eq!(m.value(), c.value(), "{c:?}");
            // Mapping back gives the original, flags and values included.
            assert_eq!(m.map_keys(|k| k - 10), *c);
        }
        let geometry: Vec<Geometry<u32>> = vec![
            Geometry::Point { at: DVec2::X },
            Geometry::Line { start: 1, end: 2 },
            Geometry::Circle {
                center: 1,
                radius: 0.5,
            },
            Geometry::Arc {
                center: 1,
                start: 2,
                end: 3,
            },
            Geometry::ControlSpline {
                control: vec![1, 2, 3],
                degree: 2,
                periodic: false,
            },
            Geometry::InterpolatedSpline {
                through: vec![1, 2],
                params: vec![0.0, 1.0],
                periodic: false,
            },
            Geometry::Ellipse {
                center: 1,
                major: 2,
                minor: 3,
            },
            Geometry::EllipticalArc {
                center: 1,
                major: 2,
                minor: 3,
                start: 4,
                end: 5,
            },
            Geometry::Conic {
                start: 1,
                control: 2,
                end: 3,
                rho: 0.3,
            },
            Geometry::Text {
                corners: [1, 2, 3, 4],
                content: "A".into(),
                font: "Liberation Sans".into(),
                aspect: 0.6,
                mirrored: true,
                outline: vec![vec![
                    TextPiece {
                        degree: 1,
                        points: vec![DVec2::ZERO, DVec2::X],
                    },
                    TextPiece {
                        degree: 2,
                        points: vec![DVec2::X, DVec2::ONE, DVec2::ZERO],
                    },
                ]],
            },
        ];
        for g in &geometry {
            let m = g.map_keys(|k| k + 10);
            let expected: Vec<u32> = g.referenced_points().iter().map(|k| k + 10).collect();
            assert_eq!(m.referenced_points(), expected);
            assert_eq!(m.map_keys(|k| k - 10), *g);
        }
    }
}
