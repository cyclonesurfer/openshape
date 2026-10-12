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

5. **Local feature size** (2026-10-12). A chord tolerance chosen for the
   whole body (the validator uses 1% of its size) is too coarse where two
   boundary curves of a face come closer than about twice the chord
   error: a sliver between a line and a shallow arc, a narrow notch, a
   hole near the outer loop, two curves leaving a vertex at a small angle.
   There the independently discretized edges cross in parameter space, a
   hole's samples can all fall between the outer loop's chord and its
   arc, a two-edge loop can get only two samples, and a sliver's two sides
   can be one chord each: the mesh is invalid, or encloses no volume,
   although the face is fine. The edge discretization therefore adapts to
   the boundary, after step 1 and before triangulation, in rounds:
   * every segment of a loop with fewer than three samples is bisected;
   * the boundary segments of each face are mapped to the scaled
     parameter space of step 3, each with a *deviation*: 1.25 times the
     largest distance from its chord of the curve's points at ¼, ½ and ¾
     of the segment (on curved faces projected by Gauss–Newton steps from
     the chord's own point, so a point over a closure or seam stays on the
     segment's side; capped at half the chord, since a segment turning
     less than the angle tolerance strays less than that and more is a
     projection artefact);
   * two segments that cross or touch other than at a shared end (the
     robust `orient2d` predicate the triangulation uses), or whose chords
     are within the sum of their deviations (so the curves may come closer
     than the chords show), conflict; so does a segment sharing an end
     with another whose far end lies within its deviation of its chord
     (the curve may pass it: the inside-out triangle of a thin lens). A
     sort along one axis keeps the pairwise test near-linear. Two segments
     of one loop joined (the shorter way round, by at most 32 segments) by
     a stretch whose segments all run within 45° of the line from its
     start to its end do not conflict by their deviations alone: segments
     turn by less than the angle tolerance, so the curve along that
     stretch is a graph over that line and the two cannot meet, though
     their bands overlap across short segments between them;
   * of a conflicting pair the segment with the larger deviation is
     bisected at its parameter midpoint on the exact curve (both when
     within 10%; by chord length when neither deviates).
   Bisection shrinks a deviation fourfold, so a gap g between two curves is
   resolved once the local segment length is about √(8g/κ) (κ the
   curvature), and two curves leaving a vertex at angle α once their first
   segments turn by less than α. Refinement edits the shared edge
   parameters, so the mesh stays watertight; only faces whose edges changed
   are re-projected. No segment shorter than the model tolerance is split,
   and refinement stops after 60 rounds or once the samples grow past
   eight times the initial count plus 100 000 (bisection is then chasing
   something it cannot resolve, and the triangulation reports the
   crossing).

   **Cusps.** Where two curves of a loop leave a vertex in directions less
   than the angle tolerance apart, the conflict tests above cannot be met:
   the curves converge to the vertex, so their deviation bands overlap at
   every scale, and independently placed samples keep the two polygons in
   the curves' order only when their spacings agree to a few percent
   (radii 3 and 2.9 tangent inside each other: within 4%). Bisection then
   chased the cusp down to segments a tolerance long whose polygons crossed
   by rounding (found by the lattice stress test). So each such vertex gets
   a direction d halfway between the two curves' directions, and each side
   a *run*: its samples outward while the distance along d grows and the
   tangent stays within 60° of d (segments turn by less than the angle
   tolerance, so the curve is a graph over d along the run). Out to the end
   of the shorter run, each side gets a sample at every distance along d
   at which the other has one (solved on the exact curve by safeguarded
   regula falsi; none closer than the tolerance to an existing sample).
   Two curves that are graphs over one line and do not meet, sampled at
   the same abscissae, have polygons that never cross: their difference is
   piecewise linear with breakpoints of one sign. Within the matched runs
   two segments from opposite sides therefore conflict only if they
   actually cross. A segment in a run that must be bisected for another
   reason is split at the midpoint of its extent along d, not of its
   parameter, so that when both sides split there they split at one
   abscissa; matching is then restored after the round. (Splitting at
   parameter midpoints put the two sides' new samples a fraction of a
   millimetre apart at a 26° corner of two circles of 1 km radius;
   matching added each to the other, the long segments beside each short
   one conflicted again, and the samples doubled every round, found by
   the scale test. An edge matched at both ends along two directions
   still gets nearby samples from its two partners where the runs
   overlap; the monotone-stretch rule above keeps the segments beside
   such a short one from conflicting, found by the lattice test.)
   This is the zero-angle analogue of Ruppert's concentric shells, which
   split the segments at a small input angle at equal distances from its
   apex so that neither keeps encroaching on the other.

Tests assert watertightness (every directed edge appears once, with its
reverse once) on the whole golden corpus and on random prisms and revolves,
and that triangle centroids of cylinders lie within the chord tolerance.

References: Shewchuk, "Delaunay refinement algorithms for triangular mesh
generation", CGTA 22, 2002; Ruppert, "A Delaunay refinement algorithm for
quality 2-dimensional mesh generation", J. Algorithms 18(3), 1995
(concentric shells at small angles); Piegl & Richard, "Tessellating trimmed
NURBS surfaces", CAD 27(1) 1995.
