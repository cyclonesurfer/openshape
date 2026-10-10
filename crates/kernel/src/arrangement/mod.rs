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

use glam::DVec2;
use serde::{Deserialize, Serialize};

use crate::{KernelError, ProfileCurve, Tolerance};

/// A point where curves end or meet.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ArrangementVertex {
    /// Position.
    pub point: DVec2,
    /// Indices of the input curves through this point, sorted, no duplicates.
    pub curves: Vec<usize>,
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
    /// Area (outer minus holes), > 0.
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
pub fn arrange(curves: &[ProfileCurve], tol: &Tolerance) -> Result<SketchArrangement, KernelError> {
    let _ = (curves, tol);
    unimplemented!("design note 0006 — implemented by the regions work package")
}

/// The boundary of the union of `regions` as trimmed profile curves: outer
/// loops counter-clockwise, holes clockwise, each loop's curves in order and
/// joined end to end (design note 0006, §3).
pub fn region_profile(
    arrangement: &SketchArrangement,
    curves: &[ProfileCurve],
    regions: &[usize],
) -> Result<Vec<ProfileCurve>, KernelError> {
    let _ = (arrangement, curves, regions);
    unimplemented!("design note 0006 — implemented by the regions work package")
}
