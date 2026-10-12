# Design note 0011 — Curve and surface intersections (Phase 3)

Status: design, before implementation. Booleans (design note 0012) and
every feature that cuts or joins bodies stand on this.

## 1. Problem

A boolean between two B-rep bodies needs, for every pair of faces that
can meet:

* **surface/surface** (SSI): the curves along which the two surfaces
  meet, each as a 3D curve and as a curve in the parameter space of each
  surface (a *pcurve*), restricted to the two faces;
* **curve/surface**: where each edge of one body pierces or touches each
  face of the other (these become vertices, and end the SSI curves);
* **curve/curve** (3D): where edges of the two bodies meet or overlap.

The kernel's surfaces today are planes, linear extrusions of a line,
circle or NURBS curve, surfaces of revolution of such a curve, and
rational NURBS surfaces; `Surface::kind` recognises cylinders, cones,
spheres and tori among them. Curves are lines, circles and rational
NURBS. Model tolerance is 1e-7 m linear, 1e-9 rad angular (ADR 0006).

What makes this hard is not the generic case (two surfaces crossing at a
healthy angle) but completeness and degeneracy: finding *every* branch,
including small closed loops; branches through singular points (cone
apex, sphere pole, a self-intersection of the curve); tangential contact,
where the surfaces touch along a curve or at a point and the problem is
ill-conditioned; and coincident patches (two faces on the same plane or
cylinder), which are not curves at all. Every failure must be an error,
never a missing or extra branch.

## 2. Literature

* N. M. Patrikalakis & T. Maekawa, *Shape Interrogation for Computer
  Aided Design and Manufacturing*, Springer, 2002: nonlinear polynomial
  solvers (the interval projected polyhedron algorithm) and intersection
  problems (the classification of methods into analytic, lattice,
  subdivision and marching, and the role of starting points and loop
  detection).
* C. L. Bajaj, C. M. Hoffmann, R. E. Lynch & J. E. H. Hopcroft, "Tracing
  surface intersections", *CAGD* 5(4), 1988: marching along the
  intersection with step size controlled by curvature, and handling of
  singular points.
* R. E. Barnhill & S. N. Kersey, "A marching method for parametric
  surface/surface intersection", *CAGD* 7, 1990: marching in the four
  parameters with Newton correction onto both surfaces.
* T. W. Sederberg & R. J. Meyers, "Loop detection in surface patch
  intersections", *CAGD* 5(2), 1988: if the normal cones (Gauss maps) of
  two patches satisfy a separability condition, their intersection has no
  closed loop; subdivision until it holds guarantees every branch meets a
  patch boundary.
* S. Krishnan & D. Manocha, "An efficient surface intersection algorithm
  based on lower-dimensional formulation", *ACM TOG* 16(1), 1997: the
  intersection of rational patches as an algebraic plane curve, its
  turning and singular points as starting points for tracing.
* T. A. Grandine & F. W. Klein IV, "A new approach to the surface
  intersection problem", *CAGD* 14(2), 1997: the intersection as the
  solution of a differential-algebraic system, with topology found from
  its critical points.
* T. W. Sederberg & T. Nishita, "Curve intersection using Bézier
  clipping", *CAD* 22(9), 1990: fat-line clipping for curve/curve and,
  by extension, curve/surface roots.
* X. Ye & T. Maekawa, "Differential geometry of intersection curves of
  two surfaces", *CAGD* 16(8), 1999: tangent, curvature and higher
  derivatives of an intersection curve from the two surfaces'
  derivatives, transversal and tangential (for step control).
* D. J. Jackson, "Boundary representation modelling with local
  tolerances", *Proc. 3rd ACM Symposium on Solid Modeling*, 1995:
  edges whose geometry is an approximation carry their own tolerance.
* Closed forms for plane and natural-quadric pairs: R. N. Goldman &
  J. R. Miller, "Combining algebraic rigor with geometric robustness for
  the detection and calculation of conic sections in the intersection of
  two natural quadric surfaces", *Proc. ACM Symp. Solid Modeling
  Foundations and CAD/CAM Applications*, 1991; and C.-K. Shene &
  J. K. Johnstone, "On the lower degree intersections of two natural
  quadrics", *ACM TOG* 13(4), 1994.
* Our own design notes 0004 (NURBS evaluation after Piegl & Tiller),
  0006 (the 2D arrangement's contact finding: subdivision of monotone
  Bézier pieces, Newton refinement, tangency finishing), whose techniques
  carry over.

## 3. Representation of results

```text
IntersectionCurve {
    curve: Curve3,          // exact (line, circle, conic as NURBS) or a fitted B-spline
    on_a: Pcurve,           // the same curve in face A's (u, v)
    on_b: Pcurve,           // and in face B's
    tol: f64,               // largest measured distance of curve and pcurves from both surfaces
    kind: Transversal | Tangent,
    ends: [End; 2],         // boundary point, closed loop, or singular point
}
```

* Exact where a closed form exists (§4.1); otherwise a cubic B-spline
  through marched points, refined until the measured deviation from both
  surfaces (checked at every knot span's interior Gauss points, not only
  at the fitted points) is within `tol ≤ 0.1 ×` the model tolerance.
  The edge built from it records `tol` (Jackson's local tolerance), so
  the validator checks the edge against the surfaces within its own
  tolerance and nothing downstream assumes it is exact.
* Pcurves are kept because trimming, tessellation and later operations
  work in parameter space; recomputing them by projection near a
  singularity is exactly where projection is unreliable.

## 4. Surface/surface

### 4.1 Closed forms first

Pairs with known exact answers are solved in closed form, after
classification (`Surface::kind`) and with tolerant tests for the special
positions (parallel, coaxial, tangent):

| pair | result |
|---|---|
| plane/plane | line, or coincident (§4.5) |
| plane/cylinder | 0–2 lines, a circle or an ellipse (exact rational), or tangent line |
| plane/cone | conic (circle, ellipse, parabola, hyperbola branches) or lines through the apex |
| plane/sphere, plane/torus | circle; torus: up to two circles when the plane contains or is perpendicular to the axis, else the general case (§4.2) |
| plane/extrusion | the profile (or its affine image) when the plane is not parallel to the direction; rulings (lines) where it is, from the 2D intersection of the profile with the plane's trace |
| plane/revolution | circles when the plane is perpendicular to the axis; meridian profiles when it contains the axis; else §4.2 |
| two cylinders, cones or spheres with a common axis, or two spheres | circles (or coincidence) |
| two parallel cylinders | 0–2 lines |

The extrusion and revolution reductions reuse the 2D curve/curve
machinery of design note 0006 (subdivision of monotone pieces, Newton,
tangency finishing), which already refuses or reports what it cannot
resolve. Natural-quadric pairs in general position (two cylinders at an
angle) are quartic space curves; they go to §4.2, which is exact to the
fitted tolerance. A degenerate quadric pair that factors into conics (the
cases of Goldman & Miller and Shene & Johnstone) is recognised and
returned exactly; tangential factorisations are §4.4.

### 4.2 The general case: subdivide, find starts, march

1. **Bounding.** Each face's surface, restricted to the face's parameter
   box, is split into rational Bézier patches (knot insertion, Piegl &
   Tiller ch. 5)
   with bounding boxes and normal cones. Pairs whose boxes are farther
   apart than the tolerance are discarded.
2. **Loop detection.** For each remaining pair, test the Sederberg–Meyers
   condition on the two normal cones. While it fails and the patches are
   larger than a size floor, subdivide both. Below the floor without
   separation means a near-tangency: §4.4.
3. **Starting points.** Every branch that satisfies the condition meets a
   patch boundary, so the starts are the curve/surface intersections of
   each patch's four boundary curves with the other patch (§5), plus
   singular points found during subdivision. Each start is refined by
   Newton's method on the 4-parameter system with a minimum-norm step.
4. **Marching** (Barnhill & Kersey; Bajaj et al.): from each start, step
   along `t = n_A × n_B` in 3D, mapped to both parameter spaces, with step
   length from the intersection curve's curvature (Ye & Maekawa) so that
   the chord deviation stays under a tenth of the tolerance, then Newton
   back onto both surfaces. A step that changes the tangent's direction
   by more than an angle limit, or whose Newton correction does not
   converge quadratically, is halved; one that would jump to another
   branch is caught by checking that the corrected point lies on the
   chord's side predicted by the curvature.
5. **Termination**: at a patch boundary (hand over to the neighbouring
   patch pair, which must find the same point as a start), at a face
   boundary, back at the start (a closed loop), or at a singular point
   (where `|t|` vanishes: §4.4).
6. **Bookkeeping.** Each start is consumed by exactly one traced branch
   (from one side or the other); a start left over, or a branch that ends
   anywhere else, is an internal inconsistency and an error.

### 4.3 Certification

A branch is accepted only when its fitted curve and both pcurves are
within `tol` of both surfaces at every knot span's Gauss points; when
consecutive marched points are monotone along the curve; and when its
ends match the starts and boundary points they were traced between.
This turns "the marching looked fine" into a measured claim. The
fallback when certification fails is a smaller step and a refit; after a
bounded number of refits it is an error.

### 4.4 Tangency and singularities

Where `|n_A × n_B|` falls below the angular tolerance the problem is
ill-conditioned: the curve's position across the surfaces is known only
to about `sqrt(tol / Δκ)`, as in the 2D near-tangencies of design note
0006, which stay within the tolerance over a stretch. The kernel
distinguishes:

* **isolated tangent point** (sphere on a plane): a single point found
  as a minimum of the distance function, returned as a vertex;
* **tangent curve** (cylinder resting on a plane, fillet faces against
  their supports): recognised in closed form for the analytic pairs of
  §4.1; for the general case, traced as a curve of minima of the distance
  (solving `n_A ∥ n_B` with the surfaces within tolerance), marked
  `Tangent`;
* **singular point** of a branch (two branches crossing, a cone apex on
  the other surface): found by the subdivision as a point where both the
  distance and `n_A × n_B` vanish; branches are traced out of it along
  the tangent directions of the quadratic approximation (Bajaj et al.).

Anything else ill-conditioned (two surfaces nearly tangent over an area
but not coincident within the tolerance) is refused with a message naming
where. Booleans then fail with that message; they never guess.

### 4.5 Coincident patches

Two faces on the same plane, cylinder, cone, sphere or torus (recognised
in closed form within the tolerance), or NURBS patches within the
tolerance of each other everywhere (tested by projection after the
bounding step), are not an intersection curve but a coincidence; the
result is the overlap of the two faces' parameter-space regions on the
common surface, which design note 0012 handles as an "on" classification.

## 5. Curve/surface

Closed forms for line/plane, line/quadric, circle/plane and circle/sphere
(and circle/quadric by reduction to a quartic in the angle, its real roots
isolated by Sturm sequences and polished by Newton's method);
otherwise subdivision of the curve's Bézier pieces against the surface's
patches (bounding boxes, then Bézier clipping of the curve against a fat
plane of the patch, after Sederberg & Nishita), and Newton on the
three-parameter system. Tangent contacts (a line tangent to a cylinder)
are recognised by `|c′ · n|` near zero at the root and finished as
tangencies, as in design note 0007 §6 for the 2D case. A curve lying in a
surface within the tolerance is a coincidence, returned as parameter
intervals.

## 6. Curve/curve (3D)

Lines and circles in closed form; otherwise subdivision with bounding
boxes and Newton on the two-parameter distance minimisation; overlaps
(coincident within the tolerance over an interval) returned as
intervals, as the 2D arrangement does.

## 7. Tests

* Closed-form pairs against independent formulas, under random rigid
  motions and scalings (1e-3 to 1e3, as the sketch tests do): the
  marching path (§4.2), forced on the same pairs, must reproduce the
  closed-form curves within its reported `tol`, which validates the
  general path against exact answers.
* Small loops: a sphere just piercing a plane, a bump of a NURBS patch
  barely crossing another, tori meeting in four loops; the count of
  branches is asserted (loop detection is what a marching-only method
  gets wrong).
* Singular points: a cone apex on a plane, two cylinders of equal radius
  with intersecting axes (two ellipses crossing at two points).
* Tangency: a cylinder on a plane, a sphere in a cylinder of the same
  radius, nearly tangent pairs swept through the tolerance (refusals
  must name where, never miss a branch silently).
* Property tests over random patch pairs: every reported point is on
  both surfaces within `tol`; the branch count is stable under
  perturbation far from tangency; the intersection of a surface with its
  own translate by more than the tolerance has the expected count.
* A failure-rate metric (refusals per thousand random pairs) printed in
  CI, as the booleans' will be.

## 8. Order of work

1. Data structures (§3), closed forms for plane pairs and the
   extrusion/revolution reductions (§4.1): enough for booleans of the
   bodies the kernel makes today (extrudes and revolves of sketch
   profiles meeting planar and cylindrical faces).
2. Curve/surface and curve/curve (§5–6).
3. The general path (§4.2–4.3) with loop detection, then tangency and
   singularities (§4.4), each with the tests above before it is used by a
   boolean.
