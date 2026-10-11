//! Planar arrangements of sketch curves: the regions users pick for extrude
//! and revolve (design note 0006).
//!
//! Curve parameters used throughout this module:
//! * [`ProfileCurve::Line`]: t ∈ [0, 1] from start to end.
//! * [`ProfileCurve::Arc`]: t ∈ [0, |sweep|], the angle travelled from the
//!   start point along the arc's own direction.
//! * [`ProfileCurve::Circle`]: t ∈ [0, 2π), the counter-clockwise angle
//!   from the +x direction.
//! * [`ProfileCurve::Spline`]: the B-spline's own knot parameter.
//!
//! Tolerance semantics: points closer than the linear tolerance are one
//! vertex; a curve passing within the tolerance of a vertex passes through
//! it; two pieces joining the same vertices and lying within the tolerance
//! of each other are one edge (an overlap, kept on the lower curve index).
//! Configurations that cannot be resolved consistently at the tolerance
//! (for example a region thinner than the tolerance) are reported as
//! [`KernelError::InvalidInput`] naming the location, never as a wrong
//! region. See `topology` for the algorithm and its consistency checks.

use glam::DVec2;
use serde::{Deserialize, Serialize};

use crate::{KernelError, ProfileCurve, Tolerance};

mod bezier;
mod curve;
mod intersect;
mod measure;
mod output;
mod topology;

/// A point where curves end or meet.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ArrangementVertex {
    /// Position: the mean of the most precise points that define the vertex
    /// (curve end points, else transversal crossings, else other contacts).
    pub point: DVec2,
    /// Indices of the input curves through this point, sorted, no duplicates.
    pub curves: Vec<usize>,
    /// The vertex's parameter on each curve in `curves` (same order and
    /// length), in the module's parameter conventions; for a circle in
    /// [0, 2π). Where a curve passes through the vertex more than once (a
    /// self-intersection, or the coinciding ends of a closed curve), the
    /// smallest of its parameters there.
    pub params: Vec<f64>,
}

/// A maximal piece of one input curve between two vertices, with no other
/// vertex in its interior.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ArrangementEdge {
    /// Index of the input curve.
    pub curve: usize,
    /// Parameter interval on that curve, `t0 < t1` (see the module docs).
    /// For a circle the interval may extend past 2π when it wraps.
    pub t0: f64,
    /// End of the parameter interval.
    pub t1: f64,
    /// Vertex at `t0`.
    pub start: usize,
    /// Vertex at `t1` (equal to `start` for a closed curve with no other
    /// vertices).
    pub end: usize,
    /// Region on the left when travelling from `t0` to `t1`, if bounded.
    pub left: Option<usize>,
    /// Region on the right, if bounded.
    pub right: Option<usize>,
}

/// A closed loop of edges.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RegionLoop {
    /// Edges in order, with `true` when traversed from `t0` to `t1`.
    pub edges: Vec<(usize, bool)>,
}

/// A bounded face of the arrangement.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Region {
    /// Outer boundary, counter-clockwise (region on the left).
    pub outer: RegionLoop,
    /// Holes, clockwise (region on the left).
    pub holes: Vec<RegionLoop>,
    /// Area (outer minus holes), > 0: Green's theorem on the exact curve
    /// pieces, with straight connectors across the tolerance-sized gaps
    /// between consecutive pieces (so defined to about the tolerance times
    /// the perimeter when curves only meet within the tolerance).
    pub area: f64,
    /// A point strictly inside, at least the tolerance from the boundary.
    pub interior_point: DVec2,
}

/// The planar arrangement of a set of curves.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SketchArrangement {
    /// Vertices.
    pub vertices: Vec<ArrangementVertex>,
    /// Edges, excluding dangling edges (which bound no region).
    pub edges: Vec<ArrangementEdge>,
    /// Bounded regions.
    pub regions: Vec<Region>,
}

/// Builds the arrangement of `curves` (design note 0006, §4).
///
/// Degenerate curves (zero length or radius, non-finite values, coordinates
/// beyond 10¹³ tolerances) are [`KernelError::InvalidInput`] naming the
/// curve index; inputs whose configuration cannot be resolved at the
/// tolerance are `InvalidInput` naming the location; internal
/// inconsistencies are [`KernelError::InvalidResult`].
pub fn arrange(curves: &[ProfileCurve], tol: &Tolerance) -> Result<SketchArrangement, KernelError> {
    if !(tol.linear.is_finite() && tol.linear > 0.0) {
        return Err(KernelError::InvalidInput(
            "the linear tolerance must be positive".into(),
        ));
    }
    topology::build(curves, tol.linear)
}

/// The boundary of the union of `regions` as trimmed profile curves: outer
/// loops counter-clockwise, holes clockwise, each loop's curves in order and
/// joined end to end (design note 0006, §3), for an arrangement built with
/// the default tolerance. See [`region_profile_with_tolerance`].
pub fn region_profile(
    arrangement: &SketchArrangement,
    curves: &[ProfileCurve],
    regions: &[usize],
) -> Result<Vec<ProfileCurve>, KernelError> {
    region_profile_with_tolerance(arrangement, curves, regions, &Tolerance::default())
}

/// [`region_profile`] for an arrangement built with tolerance `tol`.
///
/// Each boundary edge (an edge with exactly one side selected) becomes one
/// curve, trimmed to its interval and oriented with the union on the left;
/// edges between two selected regions disappear, and adjacent pieces of one
/// input curve stay separate. The profile builder merges end points closer
/// than ten tolerances, so the output is adjusted to that resolution: an
/// edge whose ends are that close (for example a closed edge) is cut in two
/// (at its parameter midpoint, unless that is near its ends), and a run of boundary edges that each
/// stay within ten tolerances of their start is left out when the run stays
/// within twenty tolerances of where it starts and ends within ten of it
/// (the curves before and after it then meet within ten tolerances, and
/// the boundary moves by at most twenty). Duplicate indices in `regions`
/// are ignored. An empty or out-of-range selection, a boundary loop or a
/// longer chain of boundary edges below that resolution, or an arrangement
/// that does not match `curves`, is [`KernelError::InvalidInput`].
pub fn region_profile_with_tolerance(
    arrangement: &SketchArrangement,
    curves: &[ProfileCurve],
    regions: &[usize],
    tol: &Tolerance,
) -> Result<Vec<ProfileCurve>, KernelError> {
    output::region_profile(arrangement, curves, regions, tol)
}

/// The index of the region whose interior contains `point` more than the
/// tolerance from the region's boundary (winding number of the exact curves
/// over the outer loop and the holes), or `None` (outside every region,
/// within the tolerance of a boundary, or an arrangement that does not
/// match `curves`). A point in a region's hole is not in that region; it is
/// in the region filling the hole.
pub fn region_containing(
    arrangement: &SketchArrangement,
    curves: &[ProfileCurve],
    point: DVec2,
    tol: &Tolerance,
) -> Option<usize> {
    output::region_containing(arrangement, curves, point, tol)
}
