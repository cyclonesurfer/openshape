//! Error type shared by every fallible constructor in the crate.

use thiserror::Error;

/// Error returned when a NURBS entity cannot be built from the given input.
///
/// Every constructor in this crate validates its input and reports problems
/// through this type instead of panicking.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum NurbsError {
    /// The degree must be at least 1.
    #[error("degree must be at least 1")]
    ZeroDegree,

    /// A knot vector needs at least two values.
    #[error("knot vector needs at least 2 values, found {found}")]
    TooFewKnots {
        /// Number of knots supplied.
        found: usize,
    },

    /// A knot value is NaN or infinite.
    #[error("knot value at index {index} is not finite")]
    NonFiniteKnot {
        /// Index of the offending knot.
        index: usize,
    },

    /// The knot sequence decreases somewhere.
    #[error("knot vector is decreasing at index {index}")]
    DecreasingKnots {
        /// Index of the first knot that is smaller than its predecessor.
        index: usize,
    },

    /// The number of knots does not equal `control points + degree + 1`.
    #[error(
        "expected {expected} knots for {control_points} control points of degree {degree}, found {found}"
    )]
    KnotCountMismatch {
        /// Degree of the entity.
        degree: usize,
        /// Number of control points (in the relevant direction).
        control_points: usize,
        /// Required number of knots.
        expected: usize,
        /// Number of knots supplied.
        found: usize,
    },

    /// Fewer than `degree + 1` control points were supplied.
    #[error("degree {degree} requires at least {required} control points, found {found}")]
    TooFewControlPoints {
        /// Degree of the entity.
        degree: usize,
        /// Minimum number of control points (`degree + 1`).
        required: usize,
        /// Number of control points supplied.
        found: usize,
    },

    /// A knot is repeated more often than the degree allows.
    ///
    /// Interior knots may appear at most `degree` times (otherwise the entity
    /// would be discontinuous); any knot may appear at most `degree + 1` times.
    #[error("knot {value} has multiplicity {multiplicity}, at most {max} is allowed")]
    KnotMultiplicityTooHigh {
        /// The repeated knot value.
        value: f64,
        /// Its multiplicity.
        multiplicity: usize,
        /// Largest multiplicity allowed at that position.
        max: usize,
    },

    /// The parametric domain `[u_p, u_{n+1}]` has zero length.
    #[error("parametric domain is empty")]
    EmptyDomain,

    /// The number of weights does not match the number of control points.
    #[error("expected {expected} weights, found {found}")]
    WeightCountMismatch {
        /// Number of control points.
        expected: usize,
        /// Number of weights supplied.
        found: usize,
    },

    /// A weight is zero, negative, NaN or infinite.
    ///
    /// For surfaces `index` is the row-major index `i * n_v + j`.
    #[error("weight at index {index} must be positive and finite, found {value}")]
    InvalidWeight {
        /// Index of the offending weight.
        index: usize,
        /// The offending value.
        value: f64,
    },

    /// A control point has a NaN or infinite coordinate.
    ///
    /// For surfaces `index` is the row-major index `i * n_v + j`.
    #[error("control point at index {index} is not finite")]
    NonFiniteControlPoint {
        /// Index of the offending control point.
        index: usize,
    },

    /// The rows of a control net (or its weights) do not all have the same length.
    #[error("control net row {row} has {found} entries, expected {expected}")]
    RaggedControlNet {
        /// Index of the offending row.
        row: usize,
        /// Length of the first row.
        expected: usize,
        /// Length of the offending row.
        found: usize,
    },

    /// Not enough data points for the requested operation.
    #[error("need at least {required} points, found {found}")]
    TooFewPoints {
        /// Minimum number of points.
        required: usize,
        /// Number of points supplied.
        found: usize,
    },

    /// A geometric argument is degenerate or out of range (zero radius,
    /// zero-length axis, empty sweep, coincident points, ...).
    #[error("invalid argument: {0}")]
    InvalidArgument(&'static str),

    /// A linear system that should be regular turned out to be singular.
    #[error("singular linear system")]
    SingularSystem,
}
