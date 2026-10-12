//! 2D geometric constraint solver for OpenDatum sketches.
//!
//! See `docs/design/0005-sketch-solver.md`. The data model ([`Sketch`],
//! [`Entity`], [`Geometry`], [`Constraint`]) is shared with the document
//! model; [`solve()`] moves the geometry to satisfy the constraints, nearest to
//! where it is, and diagnoses what is free, redundant or conflicting.
//! [`resolve_curves`] turns a sketch's curves into explicit geometry for the
//! kernel.

pub mod model;

mod build;
pub mod conic;
mod diagnose;
pub mod ellipse;
mod eq;
mod index;
mod newton;
mod offset;
mod reduce;
mod resolve;
mod solve;
mod spline;

pub use glam::DVec2;
pub use model::*;
pub use resolve::{control_spline_curve, interpolated_spline_curve, text_curves};
pub use spline::clamped_uniform_knots;

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// Options for [`solve()`].
#[derive(Clone, Debug, PartialEq)]
pub struct SolveOptions<K> {
    /// Points being dragged and their target positions (design note §7).
    pub drag: Vec<(K, DVec2)>,
    /// Iteration limit per component.
    pub max_iterations: usize,
}

impl<K> Default for SolveOptions<K> {
    fn default() -> Self {
        SolveOptions {
            drag: Vec::new(),
            max_iterations: 100,
        }
    }
}

/// Why a sketch could not be solved.
#[derive(Clone, Debug, PartialEq, thiserror::Error, Serialize, Deserialize)]
pub enum SolveError {
    /// A reference or value is invalid (wrong entity kind, missing entity,
    /// non-positive length, …). The message names the item.
    #[error("invalid sketch: {0}")]
    Invalid(String),
    /// Geometry is degenerate where a constraint needs it not to be (e.g. a
    /// zero-length line under a parallel constraint), or a curve collapses
    /// (see [`solve()`]) with no single constraint to blame, e.g. because it
    /// is already collapsed in the input. The message names the item.
    #[error("degenerate geometry: {0}")]
    Degenerate(String),
    /// The constraints cannot all be satisfied, or only by collapsing a
    /// curve; see [`Solution::conflicting`].
    #[error("constraints conflict")]
    Conflicting,
    /// No acceptable solution was reached (the iteration limit, or a result
    /// that fails the final check of [`solve()`]) and no conflict could be
    /// identified.
    #[error("solver did not converge")]
    NotConverged,
}

/// Outcome of [`solve()`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum SolveStatus {
    /// All constraints satisfied; none redundant.
    Solved,
    /// All constraints satisfied, but some are implied by others (see
    /// [`Solution::redundant`]).
    Redundant,
    /// Not solved; the returned geometry is the input geometry.
    Failed(SolveError),
}

/// Whether an entity can still move.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum EntityState {
    /// Some of its unknowns are not determined by the constraints.
    Free,
    /// Fully determined.
    Constrained,
}

/// Result of [`solve()`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Solution<K: Ord> {
    /// Outcome.
    pub status: SolveStatus,
    /// The sketch with solved geometry (unchanged input on failure).
    pub sketch: Sketch<K>,
    /// Remaining degrees of freedom.
    pub dof: usize,
    /// State of every entity (curves are free if any of their points or
    /// radius is).
    pub entities: BTreeMap<K, EntityState>,
    /// Constraints implied by earlier ones.
    pub redundant: Vec<K>,
    /// Constraints that cannot be satisfied together, or that collapse a
    /// curve (on failure).
    pub conflicting: Vec<K>,
    /// Largest constraint residual at the result, in metres or radians.
    pub max_residual: f64,
    /// Total solver iterations.
    pub iterations: usize,
}

/// Solves `sketch` (design note 0005).
///
/// A result is accepted only if every constraint holds to the model
/// tolerance (1e-7 m, 1e-9 rad — or, where rounding of the coordinates to
/// f64 alone exceeds that, as for very short directions far from the origin
/// or across long levers, to a bound on that rounding), every angle is on
/// its branch, no direction a constraint needs is degenerate and no curve
/// has collapsed at the result's own scale (so solving the result again
/// does not report it degenerate), and coordinates and radii stay within
/// ±1e6 m. Otherwise the solve fails as `NotConverged` (or `Conflicting`),
/// unless nothing is dragged and the input itself passes this check: then
/// the input is the solution, so a solved sketch stays solved.
///
/// A curve has collapsed when a line's end points meet, a circle's or
/// arc's radius vanishes, or an arc's end comes onto its start (where it
/// is not so in the input), all within 8 × 1e-10 of the sketch size. That
/// is never a solution, even when every constraint holds. It is blamed
/// like a conflict: the constraints of the collapsed curve's connected part
/// of the sketch that hold it collapsed (to first order) are tried in
/// reverse order (latest first, at most 32), and the first whose removal
/// lets that part solve with nothing collapsed is reported in
/// [`Solution::conflicting`] with status `Conflicting` (one constraint per
/// collapsed part), e.g. the vertical constraint of a line made horizontal
/// and then vertical. If no single removal helps, as for a curve already
/// collapsed in the input, the solve fails as `Degenerate` naming the
/// curve. A solve that fails while heading into a collapse (a radius
/// cannot reach zero) is blamed the same way.
///
/// A conic is collapsed when its ends meet, its control point comes within
/// that threshold of the line through its ends, or its ρ leaves
/// [`RHO_MIN`], 1 − [`RHO_MIN`]] (design note 0009 §4.6).
///
/// On failure the returned geometry is the input. For `Invalid` and
/// `Degenerate` failures nothing is solved or diagnosed: `dof` counts every
/// unknown (2 per point, 1 per circle or conic), every entity is `Free`, the lists
/// are empty and `max_residual` is 0. For `Conflicting` and `NotConverged`,
/// the diagnosis describes the least-squares state the solver reached, and
/// `max_residual` is measured on the (unchanged) input geometry.
pub fn solve<K: model::SketchKey>(sketch: &Sketch<K>, options: &SolveOptions<K>) -> Solution<K> {
    solve::solve(sketch, options)
}

/// A sketch curve with explicit geometry, for the kernel.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ResolvedCurve {
    /// Line segment.
    Line {
        /// Start.
        start: DVec2,
        /// End.
        end: DVec2,
    },
    /// Full circle.
    Circle {
        /// Centre.
        center: DVec2,
        /// Radius.
        radius: f64,
    },
    /// Counter-clockwise arc.
    Arc {
        /// Centre.
        center: DVec2,
        /// Radius.
        radius: f64,
        /// Angle of the start point.
        start_angle: f64,
        /// Counter-clockwise sweep, in (0, 2π).
        sweep_angle: f64,
    },
    /// Full ellipse `center + cos θ·major + sin θ·minor`, θ ∈ [0, 2π)
    /// (design note 0008 §3): `major` and `minor` are the semi-axis vectors.
    Ellipse {
        /// Centre.
        center: DVec2,
        /// First (major) semi-axis.
        major: DVec2,
        /// Second (minor) semi-axis.
        minor: DVec2,
    },
    /// Elliptical arc: the points `center + cos θ·major + sin θ·minor` for
    /// θ from `start_param` to `start_param + sweep_param`, which run
    /// counter-clockwise in the plane (`sweep_param` is negative when
    /// `major × minor < 0`).
    EllipticalArc {
        /// Centre.
        center: DVec2,
        /// First (major) semi-axis.
        major: DVec2,
        /// Second (minor) semi-axis.
        minor: DVec2,
        /// Eccentric anomaly of the start point.
        start_param: f64,
        /// Signed extent of the eccentric anomaly, 0 < |sweep| < 2π.
        sweep_param: f64,
    },
    /// Conic arc (design note 0009 §3): the rational quadratic Bézier
    /// from `start` to `end` with control point `control` and weights 1,
    /// `weight`, 1.
    Conic {
        /// Start.
        start: DVec2,
        /// Where the end tangents meet.
        control: DVec2,
        /// End.
        end: DVec2,
        /// Middle weight `ρ/(1 − ρ)`, positive.
        weight: f64,
    },
    /// Non-rational clamped B-spline.
    Spline {
        /// Degree.
        degree: usize,
        /// Knot vector.
        knots: Vec<f64>,
        /// Control points.
        points: Vec<DVec2>,
    },
}

/// The explicit geometry of every curve entity of `sketch` (points are
/// skipped), with its key and construction flag, in entity order.
///
/// Arcs run counter-clockwise from the start point's angle to the end
/// point's angle, with the radius of the start point. Elliptical arcs run
/// counter-clockwise between their end points' eccentric anomalies. Conics
/// carry their middle weight. Control-point splines
/// get clamped uniform knots on `[0, 1]`; interpolated splines are cubic
/// (quadratic through 3 points, a straight segment through 2) with knots by
/// averaging (Piegl & Tiller eq. 9.8, A9.1). A smooth closed spline (ADR
/// 0016) is resolved to the clamped spline it equals on one period, closed
/// exactly: a control-point one on clamped uniform knots, a through-point
/// one on the knots `0⁴, t₁ … t_{n−1}, 1⁴`. Zero-length lines, zero-radius
/// arcs, arcs that end where they start and conics whose ends meet or whose
/// three points are in line are `Degenerate`; bad references are `Invalid`.
pub fn resolve_curves<K: model::SketchKey>(
    sketch: &Sketch<K>,
) -> Result<Vec<(K, bool, ResolvedCurve)>, SolveError> {
    resolve::resolve_curves(sketch)
}

/// The current value of a dimension's measured quantity (what its value
/// would be if it were added now), for creating dimensions from geometry.
///
/// The conventions are those of the solver's equations: distances and
/// lengths in metres (point–line distances to the infinite line, line–line
/// distances from the first line's start to the second line), angles as the
/// signed angle from `a`'s direction to `b`'s in `(−π, π]`, a conic's ρ as
/// stored. The constraint's
/// stored value is ignored. Non-dimensions and bad references are
/// `Invalid`; zero-length lines where a direction is needed are
/// `Degenerate`.
pub fn measure<K: model::SketchKey>(
    sketch: &Sketch<K>,
    constraint: &Constraint<K>,
) -> Result<f64, SolveError> {
    resolve::measure(sketch, constraint)
}
