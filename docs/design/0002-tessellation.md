# Design note 0002 — Tessellation

Goals: watertight meshes, user-controlled chord and angle tolerance, robust
on trimmed and singular faces.

1. **Edges first.** Each edge is discretized once (lines: ends; circles:
   uniform steps meeting `r(1 − cos(Δ/2)) ≤ chord` and `Δ ≤ max_angle`;
   NURBS: adaptive bisection on chord deviation at ¼, ½, ¾ and tangent
   turn). End samples are pinned to the vertex points, so shared edges give
   bit-identical boundary vertices to both faces.
2. **Parameter-space boundary.** Boundary samples are projected to the face
   surface. At singular points (`|S_u|` or `|S_v|` ≈ 0: cone apex, sphere
   pole) the ambiguous coordinate is taken from the previous and next samples
   and the point is split into two parameter-space points with the same 3D
   position (welded on output).
3. **Constrained Delaunay triangulation** (`spade`, robust predicates) of the
   boundary in coordinates scaled by `|S_u|`, `|S_v|` at the domain centre.
   Inside/outside is decided by flood fill from the outer face, toggling at
   constraint edges — purely topological, so nearly collinear boundaries
   cannot misclassify slivers.
4. **Refinement.** For non-planar faces, triangles are refined by inserting
   their parameter centroid while the surface deviates from the triangle's
   plane by more than the chord tolerance (checked at the centroid and at
   interior-edge midpoints) or vertex normals differ by more than the angle
   tolerance. Slivers are never refined (their centroid could fall on a
   boundary constraint), and a refinement that splits a boundary constraint
   is reported as an error instead of producing a crack.

Tests assert watertightness (every directed edge appears once, with its
reverse once) on the whole golden corpus and on random prisms and revolves,
and that triangle centroids of cylinders lie within the chord tolerance.

References: Shewchuk, "Delaunay refinement algorithms for triangular mesh
generation", CGTA 2002; Piegl & Richard, "Tessellating trimmed NURBS
surfaces", CAD 27(1) 1995.
