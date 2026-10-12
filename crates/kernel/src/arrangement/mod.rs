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

pub(crate) use loops::ProfileLoops;

mod bezier;
mod curve;
mod intersect;
mod loops;
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

/// A point where two curves meet, see [`curve_contacts`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurveContact {
    /// The point on the first curve, `A(ta)`.
    pub point: DVec2,
    /// Its parameter on the first curve (see the module docs).
    pub ta: f64,
    /// The parameter of the matching point on the second curve.
    pub tb: f64,
}

/// A point of a curve with its direction and curvature, see
/// [`curvature_samples`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurvatureSample {
    /// Parameter (see the module docs).
    pub t: f64,
    /// The point.
    pub point: DVec2,
    /// Unit tangent, in the direction of increasing parameter.
    pub tangent: DVec2,
    /// Signed curvature (reciprocal metres): positive where the curve turns
    /// to the left of `tangent` (counter-clockwise).
    pub curvature: f64,
}

/// Where two curves meet (design note 0007, §6): transversal crossings,
/// tangencies, end points of either curve lying on the other, and the two
/// ends of an overlap, with the same intersection code and tolerance
/// semantics as [`arrange`], so tools and regions agree on where curves
/// meet. Parameters follow the module docs. Contacts closer than the
/// tolerance to one another are merged (the one with the smaller `ta` is
/// kept; for a closed first curve also across its seam), and the result is
/// sorted by `ta`. Every contact satisfies `|A(ta) − B(tb)| ≤ tolerance`.
///
/// Degenerate curves are [`KernelError::InvalidInput`] naming the curve
/// (0 for `a`, 1 for `b`), as in [`arrange`].
pub fn curve_contacts(
    a: &ProfileCurve,
    b: &ProfileCurve,
    tol: &Tolerance,
) -> Result<Vec<CurveContact>, KernelError> {
    if !(tol.linear.is_finite() && tol.linear > 0.0) {
        return Err(KernelError::InvalidInput(
            "the linear tolerance must be positive".into(),
        ));
    }
    let ca = curve::Crv::new(0, a, tol.linear)?;
    let cb = curve::Crv::new(1, b, tol.linear)?;
    let mut ctx = intersect::Ctx {
        tol: tol.linear,
        budget: topology::BUDGET,
    };
    let mut pairs = intersect::pair_contacts(&ca, &cb, &mut ctx)?;
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
    let mut out: Vec<CurveContact> = Vec::with_capacity(pairs.len());
    for (ta, tb) in pairs {
        let point = ca.point(ta);
        let near = |c: &CurveContact| (c.point - point).length() <= tol.linear;
        let seam = ca.periodic && out.first().is_some_and(near);
        if !(out.last().is_some_and(near) || seam) {
            out.push(CurveContact { point, ta, tb });
        }
    }
    Ok(out)
}

/// The point of `curve` at parameter `t` (see the module docs; `t` is
/// clamped to the curve's domain). Degenerate curves are
/// [`KernelError::InvalidInput`], as in [`arrange`].
pub fn curve_point(curve: &ProfileCurve, t: f64, tol: &Tolerance) -> Result<DVec2, KernelError> {
    let c = curve::Crv::new(0, curve, tol.linear)?;
    Ok(c.point(c.clamp(t)))
}

/// Samples of `curve` for a curvature comb: `n + 1` points evenly spaced
/// along its length, ends included, and at each spline knot
/// where the curvature jumps, its value on both sides (left first), all in
/// order of parameter. Curvature comes from the exact first and second
/// derivatives, `κ = (C′ × C″) / |C′|³`; points where the curve has no
/// direction (a vanishing first derivative) are left out. Degenerate curves
/// are [`KernelError::InvalidInput`], as in [`arrange`].
pub fn curvature_samples(
    curve: &ProfileCurve,
    n: usize,
    tol: &Tolerance,
) -> Result<Vec<CurvatureSample>, KernelError> {
    let c = curve::Crv::new(0, curve, tol.linear)?;
    let n = n.clamp(1, 100_000);
    let sample = |t: f64, right: bool| {
        let [point, d1, d2] = c.eval(t, right);
        let speed = d1.length();
        let curvature = d1.perp_dot(d2) / (speed * speed * speed);
        (speed > 0.0 && curvature.is_finite()).then(|| CurvatureSample {
            t,
            point,
            tangent: d1 / speed,
            curvature,
        })
    };
    // Length along the curve: Simpson's rule on the speed over a grid
    // (lines and arcs move at constant speed; each spline span is cut into
    // 64), then Newton's method for the parameter at each target length.
    let speed = |t: f64, right: bool| c.eval(t, right)[1].length();
    let length = |t0: f64, t1: f64| {
        let m = 8;
        let h = (t1 - t0) / m as f64;
        let sum: f64 = (0..=m)
            .map(|i| {
                let w = if i == 0 || i == m {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                w * speed(t0 + h * i as f64, i < m)
            })
            .sum();
        sum * h / 3.0
    };
    let mut grid = vec![c.lo];
    let spans = c.spans();
    if spans.is_empty() {
        grid.push(c.hi);
    } else {
        for s in spans {
            for k in 1..=64 {
                grid.push(s.t0 + (s.t1 - s.t0) * k as f64 / 64.0);
            }
        }
    }
    let mut len = vec![0.0];
    for w in grid.windows(2) {
        let l = len[len.len() - 1] + length(w[0], w[1]);
        len.push(l);
    }
    let total = len[len.len() - 1];
    let mut out = Vec::with_capacity(n + 1);
    let mut j = 0;
    for k in 0..=n {
        let t = if k == 0 {
            c.lo
        } else if k == n {
            c.hi
        } else {
            let s = total * k as f64 / n as f64;
            while j + 2 < len.len() && len[j + 1] < s {
                j += 1;
            }
            let (a, b) = (grid[j], grid[j + 1]);
            let f = if len[j + 1] > len[j] {
                ((s - len[j]) / (len[j + 1] - len[j])).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t = a + (b - a) * f;
            for _ in 0..3 {
                let v = speed(t, true);
                if v <= 0.0 {
                    break;
                }
                t = (t - (len[j] + length(a, t) - s) / v).clamp(a, b);
            }
            t
        };
        out.extend(sample(t, k < n));
    }
    // Both sides of the knots where the curvature jumps.
    for w in spans.windows(2) {
        let t = w[0].t1;
        if let (Some(a), Some(b)) = (sample(t, false), sample(t, true)) {
            let jump = (a.curvature - b.curvature).abs();
            if jump > 1e-9 * a.curvature.abs().max(b.curvature.abs()) {
                let at = out.partition_point(|s| s.t < t);
                let past = out.partition_point(|s| s.t <= t);
                out.splice(at..past, [a, b]);
            }
        }
    }
    Ok(out)
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
/// are ignored. An empty or out-of-range selection, a boundary loop below
/// that resolution (or within twice it of a point once its edges below it
/// are left out), a longer chain of boundary edges below it, or an
/// arrangement that does not match `curves`, is
/// [`KernelError::InvalidInput`].
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
