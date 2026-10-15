# Design note 0002 — Tessellation

Goals: watertight meshes, user-controlled chord and angle tolerance, robust
on trimmed and singular faces.

1. **Edges first.** Each edge is discretized once. Lines: their ends.
   Circles: uniform steps of at most the angle tolerance and of
   `Δ = 4 asin √(chord / 2r)`, the step whose sagitta `r(1 − cos(Δ/2))`
   is the chord tolerance (written so, without the cosine's cancellation
   when the tolerance is a tiny part of the radius). NURBS: adaptive
   bisection while the curve strays from the chord by more than the
   tolerance at ¼, ½ or ¾ of a piece, or, for a piece longer than the
   tolerance, its tangent turns by more than the angle tolerance. An edge
   needing more samples than the vertex budget (two million) is refused
   rather than left to exhaust the machine. Until 2026-10-15 a circle's
   step had a floor of 10⁻⁴ rad, so a circle of 1 km held to a chord of a
   micrometre strayed by a quarter more than that, and the angle rule kept
   splitting NURBS pieces far shorter than the tolerance when the angle
   tolerance was tiny (found by review; a 1 km circle at 10⁻⁶ now takes
   35,126 samples). End samples are pinned to the vertex points, so shared
   edges give bit-identical boundary vertices to both faces.
2. **Parameter-space boundary.** Boundary samples are projected to the face
   surface. At singular points (`|S_u|` or `|S_v|` ≈ 0: cone apex, sphere
   pole, the centre of a disc of revolution) the ambiguous coordinate is
   taken from the previous and next samples and the point becomes a run of
   parameter-space points with the same 3D position (welded on output),
   on a surface of revolution in steps of at most the angle tolerance and
   π/4 of the angle `u` (2026-10-15; two points spanning half a turn made
   fans of triangles from the apex that refinement took tens of passes to
   break up). Such a point has no normal of its own: `S_u × S_v` there is
   the rounding of a vanishing derivative and can point either way (a
   disc's centre, its radius 1.6e-17 off zero, got the normal facing in);
   each of its points in the parameters takes the surface's normal a step
   from it toward the middle of the face, and the welded vertex their
   mean (each had taken one meridian's: at a slender cone's apex, across
   a face of half a turn, that normal faced away from triangles on the
   far side, found by review). The runs, like the edges, are held to the
   vertex budget (a tiny angle tolerance asked for millions of points and
   aborted the process, found by review).
3. **Constrained Delaunay triangulation** (`spade`, robust predicates) of the
   boundary in coordinates scaled by `|S_u|`, `|S_v|` at the domain centre.
   Inside/outside is decided by flood fill from the outer face, toggling at
   constraint edges — purely topological, so nearly collinear boundaries
   cannot misclassify slivers.
   On a ruled face of revolution or extrusion (a cylinder or a cone,
   straight along `v`) the chord error of a triangle depends only on its
   width across the rulings, which the tolerance holds to
   `h = √(8 · chord · R)` (`R` the radius's scale), so `v` is scaled down
   until a triangle `h` wide spans the face along `v` (by at most 256;
   2026-10-15). Scaled isotropically, a face cut across its rulings (by
   an ellipse, a conic or a slanted line, which booleans now leave) was
   filled with triangles `h` wide each way: an obliquely cut rod meshed
   to 18,148 triangles and sixteen flattened pins to 345,000; they now
   take 624 and 9,740, following the rulings as a face bounded by its
   parameter lines does. Refinement measures the error in space, so the
   mesh meets the tolerance either way.
4. **Refinement.** For non-planar faces, a triangle is refined while it
   strays from the surface by more than the chord tolerance, faces against
   the surface's normal at its centroid (a long sliver folded over), or a
   vertex normal is more than the angle tolerance from that normal.

   Its distance from the surface is measured at thirteen points: the
   quarter points of its edges, its centroid, and the points halfway from
   the centroid to the corners. Each is measured from the surface point
   at the parameters under it, along the surface's normal there (an offset
   along the surface is the parameters' slack, not the mesh's error), and
   signed by the side of the triangle the surface is on. From an apex or a
   pole an edge runs along the other end's ruling or meridian, so the free
   parameter of the points under it is that end's (the parameters as far
   between are a point beside the chord, where the distance along the
   normal overstates the chord's: a polar cap took 10,424 triangles that
   way against 7,262).

   Between the points the error is estimated by the quartic that is zero
   at the corners and takes the measured values at the other twelve
   points of the quartic Lagrange triangle (the quarter points of the
   edges and the three inside; the centroid is a further check). This is
   the a posteriori estimate of a linear element's error from a
   higher-degree interpolant of it (Bank and Smith use the quadratic):
   exact where the surface's height over the triangle is a quartic, close
   to it where the surface is smooth at the triangle's scale. Its largest
   value is searched on a grid of twelfths of the triangle, which holds
   the thirds, where a cubic error along an edge peaks (at 4/27, above the
   9/64 of the nearest point measured), and refined by a local search,
   down to 1/384, from each of the grid's local peaks within a tenth of
   its largest. No such quartic exceeds 3.36 times its largest value at
   the twelve points (the sum of the magnitudes of their Lagrange
   functions peaks at 3.3536, near the corners), so the search is not made
   when that bound is under the tolerance, and it stops at the first point
   past it; meshing takes the same time as with the quadratic. A boundary
   segment's samples hold its chord to the tolerance (step 1), so its
   points are not tested; their values, which the parameters between its
   ends only approximate for a curved boundary, are capped at the
   tolerance in the estimate, and the search leaves the segment out (no
   point added inside moves it).

   The worst failing edge is split at its middle if it is the triangle's
   longest (in the scaled parameters), as in Rivara's longest-edge
   bisection, which keeps the triangles' angles from shrinking; otherwise
   the centroid is added, and likewise for a peak of the estimate on an
   edge. Each pass tests only the triangles the last one made (inserting a
   point re-triangulates only the triangles about it) and those that
   failed and are still there, in order of their corners, so the result
   does not depend on a hash table's order. A pass that adds no point
   while triangles fail is an error, and so is a mesh past the vertex
   budget; there is no limit on passes (refinement moves as a front
   through a thin wedge of a sphere or a torus, a row every few passes,
   and a limit of 64 refused wedges that had meshed before, found by
   review; a wedge of 20° at a chord of 10⁻⁶ now meshes to 814,902
   triangles in 6.6 s).

   A triangle that faces against the surface, or whose normals turn too
   far, beside a boundary segment that is its longest side is not split
   inside: points added inside close in on the segment's middle in the
   parameters, which can lie off its chord, and settle there folded over.
   The segment is bisected on its curve instead (Delaunay refinement
   splits an encroached segment rather than crowd it: Ruppert 1995) and
   the faces using it are meshed again, for at most 60 rounds; a face
   whose segments there are too short to bisect is refused (found by
   review: triangles folded over by a cut near a ball's pole and by an
   apex). Tiny triangles and slivers (an area under (chord/100)² in the
   scaled parameters, or thinner than 10⁻⁸ of their longest side) are not
   split, as their centroid could fall on a boundary constraint within
   rounding. One whose normals turn too far is refused. One that faces
   against the surface (three points nearly in a line in the parameters,
   the surface's bend carrying the middle one across the chord of the
   other two in space) is mended once refinement is done by flipping the
   edge it shares with a neighbour: the longest such edge for which both
   new triangles turn counter-clockwise in the parameters and pass the
   tests above (Lawson's edge flip). One that no flip mends is refused
   (found by the near-tangent apex fuzz: a frustum cut by a plane nearly
   tangent to it by its apex failed validation with "too small to split").
   A triangle whose corners lie in a line within the model tolerance (its
   height under 10⁻⁷ m) has no side to face, and is not tested for one:
   needles at a pointed cone's apex, their corners on nearly one ruling,
   faced either way by rounding, and no flip could mend them (found by the
   fuzz of planes across cones at fifty times its size).

   **The chord tolerance is estimated, not bounded.** Where the surface
   bends faster than a quartic between the points measured, the mesh can
   stray a little further. A rigorous bound needs bounds on the surface's
   second derivatives over each triangle (Filip, Magedson and Markot,
   from the control net of a NURBS's derivatives); it is not built. A
   unit test measures exact distances on a grid of 24ths of every
   triangle: a unit ball at a chord of 10⁻³ 0.9989 tolerances, two polar
   caps 0.9888 and 0.9994, a nearly closed torus at 10⁻² 0.9948, and
   solids turned from two S-shaped splines at 10⁻² 0.9607 and 0.9862. A
   probe on a grid of twelfths against closest points (2026-10-15): an
   S-spline at 10⁻³ 0.9989, a torus wedge of 4.5 rad 0.9933, the ball at
   the defaults 0.9998, a polar cap at 10⁻⁵ 1.0000 (561,598 triangles).
   With the quadratic through the edges' middles instead (the estimate
   until 2026-10-15), the S-spline strayed 5.2% and the torus 1.3%;
   without an estimate between the points, the ball 3.3% at the defaults
   (found by review).

   Until 2026-10-15 only the centroid and the middles of interior edges
   were measured, by the surface's height over the triangle's plane at
   the parameters' midpoints, and refinement stopped after twelve passes
   whether or not it had converged. Review measured meshes against the
   surfaces' closest points at seven points of every triangle (a unit
   test now does at fifteen, for two cones, a polar cap and 43 cut cones
   and cylinders at three tolerances, and a ball at two): a pointed cone
   was 16 tolerances off by its apex at the default tolerance and 159 at a
   chord of 10⁻⁵, faces cut across their rulings a quarter off along their
   untested diagonals, long slivers on frustums folded over, and a ball
   5.2% off (8.2% at 10⁻⁵), mostly 10° to 30° from its poles, with a polar
   cap 45% off at 10⁻⁵. All are now within the tolerance and face out.
   Which split takes fewest triangles depends on the body (review, on the
   code before the quartic estimate): splitting the worst failing edge at
   its middle whatever its length meshed the time guards' ring to 78,918
   triangles, the centroid to 115,432 and the longest-edge rule to
   117,050, while on the dimpled plate the longest-edge rule took fewest
   (151,950, against 174,822 and 152,966). The longest-edge rule is kept
   for the triangles' shape, not their count.

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
   the scale test.)
   This is the zero-angle analogue of Ruppert's concentric shells, which
   split the segments at a small input angle at equal distances from its
   apex so that neither keeps encroaching on the other.

   **Groups of cusps.** Where the edges of three faces leave one vertex
   within the angle tolerance of each other (a box's edge grazing a
   cone, its two faces' curves on the cone leaving the point where it
   pierces it), the vertex is a cusp of each face, and each edge belongs
   to two of them. Matched along each pair's own direction, the sample a
   cusp gives an edge lies at a distance along the next cusp's direction
   that its partner lacks; each matching added a sample for each sample
   in the runs, and while any conflict kept the rounds going the samples
   doubled every round up to the limit (a sliver a box cut from a
   frustum: 152,644 triangles and 5 s in validation's coarse mesh, found
   by the conic fuzz at fifty times its size, design note 0012 §15). So
   cusps that share a side are matched along one direction, the mean of
   the directions of all their sides: every sample added keeps its
   distance along it, so each side ends with the distances of all the
   sides it is matched with (within each pair's runs) and matching
   settles; it is repeated until it adds nothing, which takes at most
   one pass per cusp after the first. For the same reason an edge with
   cusps at both ends is matched from each end only to its parameter
   midpoint: along two directions, the samples the two ends add would
   feed each other where the runs overlapped (earlier, with matching
   done once a round, that only put nearby samples on the edge, and the
   monotone-stretch rule kept the segments beside them from conflicting,
   found by the lattice test). Beyond the midpoint such an edge is
   refined by its deviation like any other.

Tests assert watertightness (every directed edge appears once, with its
reverse once) on the whole golden corpus and on random prisms and revolves,
that triangle centroids of cylinders lie within the chord tolerance, that
meshes of cones, balls and cut faces stay within it at fifteen points of
every triangle and face out, that balls, caps, a torus and turned
S-splines stay within it on a grid of 24ths against exact distances, that
the quartic's peak search finds known peaks, that folded slivers are
flipped, and that edges of large circles meet fine tolerances or are
refused.

A ball meshes to about 2.6 times as many triangles as an even mesh at the
same tolerance would take (128,000 against about 48,400 at the defaults,
an even mesh's triangle being equilateral with its sag at the
tolerance). Review measured it in bands of equal area: 4.5 times an even
mesh within a fifth of the radius of the poles, where the sphere's
parameters stretch the parameter-plane Delaunay triangles thin, and still
1.8 times at the equator, where they do not; there the triangles' shape
(three quarters of an equilateral's area for their sag) and refinement
stopping short of the tolerance (a mean sag of 0.74 of it) cost the rest.
Removing the stretch alone would leave about 1.8 times. Not yet improved.

References: Shewchuk, "Delaunay refinement algorithms for triangular mesh
generation", CGTA 22, 2002; Ruppert, "A Delaunay refinement algorithm for
quality 2-dimensional mesh generation", J. Algorithms 18(3), 1995
(concentric shells at small angles); Piegl & Richard, "Tessellating trimmed
NURBS surfaces", CAD 27(1) 1995; Bank and Smith, "A posteriori error
estimates based on hierarchical bases", SIAM J. Numer. Anal. 30(4), 1993,
and "Mesh smoothing using a posteriori error estimates", SIAM J. Numer.
Anal. 34(3), 1997; Lawson, "Software for C¹ surface interpolation", in
Mathematical Software III, Academic Press, 1977 (edge flips); Filip,
Magedson and Markot, "Surface algorithms using bounds on derivatives",
CAGD 3(4), 1986;
Rivara, "Algorithms for refining triangular grids suitable for adaptive
and multigrid techniques", Int. J. Numer. Meth. Eng. 20(4), 1984.
