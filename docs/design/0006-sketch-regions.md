# Design note 0006 — Regions from sketch curves

Module: `crates/kernel/src/arrangement.rs`, exposed through
`Kernel::sketch_regions` and `Kernel::region_profile`.

## 1. Problem

Users draw curves that cross, touch, overlap and dangle. Extrude and revolve
need closed, non-self-intersecting loops. The sketch's regions are the
bounded faces of the planar arrangement of its non-construction curves; the
user picks one or more of them, and the profile is the boundary of their
union.

## 2. Literature

* Planar arrangements and the doubly-connected edge list: de Berg, Cheong,
  van Kreveld & Overmars, *Computational Geometry*, 3rd ed., ch. 2 (map
  overlay, face extraction, hole assignment).
* Curve intersection: Sederberg & Nishita, "Curve intersection using Bézier
  clipping", *CAD* 22(9), 1990; Patrikalakis & Maekawa, *Shape Interrogation
  for CAD/CAM*, 2002, ch. 5; closed forms for line and circle pairs.
* Robust predicates: Shewchuk, "Adaptive precision floating-point arithmetic
  and fast robust geometric predicates", 1997 (the `robust` crate, already
  used by the kernel).

## 3. Input and output

Input: `&[ProfileCurve]` (lines, arcs, circles, B-splines) in plane
coordinates, and the model tolerance. Output, a `SketchArrangement`:

* `vertices`: points where curves end or meet, each with the sorted list of
  input curves through it.
* `edges`: maximal curve pieces between vertices: input curve index, its
  parameter interval, start and end vertex, and the region on each side.
* `regions`: bounded faces: boundary loops (outer counter-clockwise, holes
  clockwise) as sequences of (edge, forward/backward), signed area, and an
  interior point (a point strictly inside, at least a tolerance away from
  the boundary).

`region_profile(arrangement, curves, selected)` returns the boundary of the
union of the selected regions as `ProfileCurve`s: edges with exactly one
side selected, trimmed to their intervals and oriented with the union on
the left. Those curves feed `extrude`/`revolve` unchanged.

## 4. Algorithm

1. **Intersections** for every pair of curves whose bounding boxes overlap
   (sweep over x-sorted boxes): closed form for line–line, line–circle and
   circle–circle (arcs as their circles, then filtered by angle range with
   tolerance); Bézier clipping for anything involving a spline, refined by
   Newton to 1e-12 relative. Tangential contacts are found as double roots
   and kept as one intersection (refined to the point of contact, design
   note 0007 §6). **Overlaps** (collinear segments, arcs on
   the same circle, a spline overlapping itself or another curve within
   tolerance over an interval) are detected explicitly; the overlapping
   interval is kept once, on the curve with the lower index.
2. **Self-intersections** of splines by subdividing at curvature extrema and
   intersecting the pieces.
3. **Vertices**: all end points and intersection points, merged when closer
   than the tolerance (union–find); each records its incident curves.
4. **Split** curves at their vertices' parameters into edges; drop edges
   shorter than the tolerance.
5. **Prune** dangling edges: repeatedly remove edges with an end of degree 1.
   They cannot bound a region.
6. **Order** edges around each vertex by outgoing tangent angle; ties
   (tangent edges) are broken by signed curvature, then by second-order
   offset sampling. Exact orientation tests use the `robust` predicates.
7. **Trace faces**: next half-edge = the previous one's twin's clockwise
   neighbour. Cycles with positive signed area (computed exactly from the
   curves: polygon terms plus circular-segment and spline Green's-theorem
   terms) are region outer loops; negative cycles are outer boundaries of
   connected components and become holes of the smallest region containing
   them (point-in-region with a winding number over the exact curves), or
   are discarded if unbounded.
8. **Interior point**: from the region's triangulation of a fine polygonal
   approximation, the centroid of the largest triangle whose centroid is at
   least the tolerance from the boundary.

Any inconsistency (an odd number of crossings found where an even number is
required, a face whose area doesn't match the sum of its pieces) is a
`KernelError`, never a silently wrong region.

## 5. Naming regions (sketch level)

The features crate maps curve indices to sketch entity ids and derives a
`RegionKey` per region: for every boundary edge, the entity id, the
direction, and the keys of its two end vertices (each the sorted entity ids
meeting there plus an ordinal when the same set meets more than once,
ordered along the lower entity's parameter). Two crossing circles produce
three regions with distinct keys. After an edit, a stored key resolves to
the region whose key matches exactly; otherwise to the region with the
largest key overlap that contains the stored interior point; otherwise the
feature fails with "profile region no longer exists" (never a guess).

## 6. Tests

* Rectangle with diagonals: 4 regions, areas exact.
* Two crossing circles: 3 regions; lens area by closed form.
* Circle tangent internally and externally to another; line tangent to a
  circle (tangential vertex, no spurious sliver regions).
* Overlapping collinear segments and coincident arcs: no duplicate edges.
* Dangling lines and construction-only geometry produce no regions.
* Nested shapes: a circle inside a square inside a circle (holes assigned
  correctly); disjoint islands.
* Spline crossing a line several times; a self-intersecting figure-eight
  spline (2 regions).
* Property tests: random polygons and circles; the sum of region areas
  equals the area of the union computed independently; every region's
  boundary profile extrudes to a valid solid whose volume is area × height.
* Scale invariance at 1e-4 and 1e3.
