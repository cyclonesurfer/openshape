# Design note 0001 — B-rep data model and sweep construction

## Data model

Arena-based boundary representation (Mäntylä 1988, ch. 10; Weiler's
half-edge structures): `Vertex`, `Edge` (curve + parameter interval +
end vertices), `Coedge` (directed use of an edge), `Loop`, `Face` (surface +
sense + loops, outer first), `Shell`, `Solid`. Geometry is attached to
edges (3D curves) and faces (surfaces); vertices carry points. Parameter
curves (pcurves) are not stored yet — edge points are projected onto faces
when needed; pcurves will be added with booleans (Phase 3), where trimming
needs them.

Orientation: a face's outward normal is `S_u × S_v` times its sense; loops
are counter-clockwise seen from outside (material on the left).

## Validation (every operation)

1. Index integrity. 2. Edge interval non-empty, curve ends on its vertices.
3. Loops closed (end of each coedge = start of the next). 4. Every edge
lies on each adjacent face's surface (sampled, 100× model tolerance).
5. Closed manifold shells: each edge used once forward and once reversed per
shell. 6. Euler–Poincaré `V − E + F − R = 2(S − H)` with integral `H ≥ 0`.
7. Outer shells enclose positive volume, void shells negative (coarse
tessellation). A failing check discards the result with an error.

## Profiles

Curves are converted to oriented segments; arcs over π and closed splines
are split (ADR 0008). End points are clustered within 10× tolerance into
vertices, every vertex must have degree 2 (open/branching profiles are
rejected), loops are walked, intersecting loops are rejected (below), and
nesting depth decides outer loops (even depth, CCW) and holes (odd depth,
CW), from the arrangement's exact loop areas and winding numbers. Each
outer loop with its direct holes is a region; each region becomes a solid
(lump).

Intersection check (exact since 2026-10-12; before, 48-point polylines
per segment, which refused curves closer than their chord error and could
in principle miss two crossings inside one chord). Every pair of profile
curves is intersected with the sketch arrangement's contact finder
(design note 0006 §4: closed forms for lines and circles, Bézier clipping
for splines, self-contacts of splines), at the model tolerance. A contact
is allowed only where the two curves (or the two ends of one closed
curve) end at the same profile vertex and the contact is part of that
vertex: within the vertex resolution (10 tolerances) of it, or at the end
of a run along which the two curves stay within twice the tolerance of
each other all the way from the vertex (two curves leaving a vertex
tangentially, as a horn between a line and a tangent arc). Every other
contact refuses the profile, naming where: crossing curves, a hole
touching its outer loop, a curve touching another away from a shared
vertex.

## Extrude

For each region: bottom/top copies of every profile vertex and edge,
vertical line edges, one side face per profile segment (plane for lines,
linear extrusion surface otherwise), and two planar caps. Face senses are
chosen geometrically (outward = tangent × plane normal), loops reversed to
match.

## Revolve

The profile must lie on one side of the axis (touching allowed at vertices
and along whole straight edges). The sweep angle θ is split into
k = ⌈θ/π⌉ pieces (k ≥ 2 for full turns). Off-axis vertices get one copy per
piece boundary and one circular edge per piece; on-axis vertices are shared
singular vertices. Each non-axis profile edge produces k faces on surfaces of
revolution whose frame starts at the piece boundary, so every face's angular
domain is [0, ≤π]. Straight edges on the axis produce no face (full turn) or
a single edge shared by both planar caps (partial turn).

Surface kinds are recovered by classification (plane, cylinder, cone,
sphere, torus) for downstream algorithms and exporters.
