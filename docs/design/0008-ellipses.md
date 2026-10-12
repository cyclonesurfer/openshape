# Design note 0008 — Ellipses and elliptical arcs

Status: accepted and implemented (2026-10-12), except the tangencies of
§4.5 and the trim, extend and split of §8, which are refused with a
message. Extends design notes 0005 (solver), 0006 (regions) and 0007
(editing tools).

## 1. Scope

Ellipses and elliptical arcs as sketch entities: their data, the solver's
equations and diagnostics, which constraints apply to them and which are
refused, how the kernel receives them exactly, the drawing tools, how the
editing tools treat them, and the tests that prove them. Conics in general
(parabolic and hyperbolic arcs, the conic tool) come later and reuse §6.

## 2. Literature

* Piegl & Tiller, *The NURBS Book*, 2nd ed.: §7.2–7.3 (conics as rational
  quadratic Bézier curves; the weight of the middle control point), §7.5
  and A7.1 (circular arcs of up to 360° as rational quadratic B-splines,
  split into at most four segments of ≤ 90°), and the affine invariance of
  rational curves (§4.2: transforming the control points transforms the
  curve; weights are unchanged).
* D. Eberly, "Distance from a point to an ellipse, an ellipsoid, or a
  hyperellipsoid", Geometric Tools, 2013: the closest point by bisection on
  a monotone function of one variable, robust near the axes. Considered for
  §7 and not needed there (see §7).
* J. Kiefer, "Sequential minimax search for a maximum", *Proc. AMS* 4,
  1953 (golden-section search), and W. Wang, H. Pottmann & Y. Liu,
  "Fitting B-spline curves to point clouds by curvature-based squared
  distance minimization", *ACM TOG* 25(2), 2006, §3 (re-projecting a
  curve parameter onto the foot point), as the solver already uses them
  for splines (0005).
* L. N. Trefethen & J. A. C. Weideman, "The exponentially convergent
  trapezoidal rule", *SIAM Review* 56(3), 2014: the reference perimeters
  in the kernel tests (§9).
* R. T. Farouki & C. A. Neff, "Analytic properties of plane offset curves",
  *CAGD* 7, 1990: the offset of an ellipse is an algebraic curve of degree
  eight, not a conic.
* R. Schneider, *Convex Bodies: The Brunn–Minkowski Theory*, 2nd ed., §1.7
  (support functions), for §4.3.
* Design note 0005 §3 (hidden curve parameters, scaling), §9 (collapse
  checks); design note 0007 §9–§11 (derived images).

## 3. Data

```text
Geometry::Ellipse { center, major, minor }
Geometry::EllipticalArc { center, major, minor, start, end }
```

All five are point entities. With A = major − center and B = minor −
center (the semi-axes), the curve is

  E(θ) = c + cos θ · A + sin θ · B,  θ ∈ [0, 2π),

θ being the eccentric anomaly. The solver holds A ⟂ B with an implicit
equation, as it holds an arc's end on its circle (0005 §3), so A and B are
the ellipse's axes and |A|, |B| its semi-axis lengths. The names say which
axis the user drew first; nothing forces |A| ≥ |B| (the curve is the same
ellipse either way, and refusing a drag past a circle would be arbitrary).

An elliptical arc is the part of that ellipse running counter-clockwise in
the sketch plane from `start` to `end`, like a circular arc. Its defining
ellipse is not an entity of its own. `start` and `end` must lie on the
ellipse (§4.2).

Handedness: cross(A, B) > 0 when the minor end is left of the major axis.
Then increasing θ runs counter-clockwise; otherwise clockwise. Both occur
(a mirror image has the other handedness), so every consumer uses the sign
of cross(A, B) rather than assuming one. A solve cannot flip it: B would
have to pass through zero length, a collapse that 0005 §9's final check
rejects.

Degrees of freedom: an ellipse has 6 point coordinates less one implicit
equation = 5 (centre, orientation, two radii); an elliptical arc 10 − 1 −
2 = 7 (also its two end angles).

## 4. Equations

### 4.1 Axes

`dot(A, B)/(‖A‖‖B‖) = 0`, the existing perpendicular residual. A
zero-length axis is a collapse (0005 §9), like a zero radius.

### 4.2 Points on the ellipse

A point p on an ellipse gets a hidden parameter θ (one unknown) and two
residual rows

  E(θ) − p = (1 − cos θ − sin θ) c + cos θ · m + sin θ · n − p,

where m and n are the major and minor end points. For fixed θ the residual
is linear in the points, so its gradients are constant multiples of the
identity; ∂/∂θ = −sin θ · A + cos θ · B = E′(θ). Rows are lengths, scaled
as every length residual (0005 §3). θ is periodic and unbounded: unlike a
spline parameter it is never clamped. Its start value is the parameter of
the point of the ellipse closest to p: the squared distance is sampled
around the whole ellipse, each local minimum of the samples refined by
golden-section search between its neighbours and then by safeguarded
Newton on (E(θ) − p)·E′(θ) = 0, the nearest winning (as for splines,
0005). When Gauss–Newton stalls, a parameter is re-seeded the same way if
that is strictly closer (Wang, Pottmann & Liu §3), which lets it jump past
a fold. Net effect: −1 DOF, as for a point on any curve.

Where a parameter is read back from a point already on the ellipse (an
arc's ends, when the solved sketch is resolved into curves), the exact
inverse θ = atan2(B·(p − c)/‖B‖², A·(p − c)/‖A‖²) is used.

This is `PointOnCurve(p, ellipse or elliptical arc)` (on the full ellipse,
as `PointOnCurve` on an arc means its full circle), and it is also how an
elliptical arc holds its own `start` and `end`, with two hidden parameters
the arc owns. An end coincident with a point already on the same ellipse
is not held again (§8.1).

### 4.3 Line tangent to an ellipse

Let the line pass through s with direction d and unit left normal
n̂ = perp(d)/‖d‖. Over the ellipse, n̂·(E(θ) − c) = cos θ (n̂·A) + sin θ (n̂·B)
has maximum h = √((n̂·A)² + (n̂·B)²) (Cauchy–Schwarz; h is the support
function of the centred ellipse in direction n̂, Schneider §1.7). The line
touches the ellipse exactly when the centre's signed distance from it is
±h:

  cross(d, c − s)/‖d‖ − σ h = 0,  σ = +1 if the centre is left of the line.

This is the `LineCircle { center_left }` tangency with h in place of the
radius: one row, no hidden parameter, smooth wherever h > 0 (always, as A
and B are non-zero and perpendicular). As for circular arcs, an elliptical
arc is tangent where its full ellipse is.

### 4.4 Tangency at a shared end point

`Tangent { kind: AtPoint }` with an elliptical arc uses the arc's tangent at
the shared point p in the existing normalized-cross residual (0005 §4),
written without the end's hidden parameter. With w = p − c, the
coordinates of w in the axis frame are κ = w·A/‖A‖² and σ = w·B/‖B‖²,
which are cos θ and sin θ when p is on the ellipse, so

  T(p) = −σ A + κ B

equals E′(θ) there. T depends on the four points only, so the row has no
θ column, and the tangency means the same whichever equation (the arc's
own end equation or a coincidence) puts p on the curve.

### 4.5 Tangency between an ellipse and a circle or another ellipse

No closed form (the contact solves a quartic). Formulation for later: a
hidden contact parameter on each ellipse, with the circle's centre at
E(θ) ± r N(θ), N the unit normal (two rows, one hidden unknown), and for
two ellipses E₁(θ₁) = E₂(θ₂) with parallel tangents (three rows, two hidden
unknowns). Until it is built these tangencies are refused (§5), and so are
sketch fillets that would need them.

### 4.6 Derived images

Mirror and both patterns map an ellipse to an ellipse whose points are the
images of its points: reflections, translations and rotations are affine,
and the image of E(θ) is the image ellipse at the same θ (with handedness
reversed by a reflection, so an image arc starts at the image of its
source's end, as for circular arcs). An image's own implicit equations
(axes perpendicular, ends on the ellipse) are implied by its source's and
are not added, as for image arcs today (0007 §9); an image arc gets no
hidden parameters.

## 5. Constraints

Supported on ellipses and elliptical arcs: coincident (their points),
point on curve (§4.2), tangent to a line (§4.3), tangent at a shared end
(§4.4), mirror and both patterns (§4.6), and everything on their points
and on lines between them. The ellipse tool makes the two semi-axes
construction lines (§7), so horizontal, vertical, parallel, perpendicular,
angle and length apply to the axes as to any line.

Refused before solving, with an error naming the constraint: radius,
diameter, equal, concentric, normal, distance to the curve, offset, and
tangency with a circle, arc or another ellipse not at a shared end (§4.5).
Concentric with a circle can be had with a coincident constraint on the
centres.

## 6. The kernel receives exact rational quadratics

Options: (a) analytic ellipse curves in the kernel next to its lines,
circles and arcs; (b) exact rational quadratic B-splines, which the kernel
already supports in profiles, regions, extrude and revolve, with tested
exact areas and volumes (a 9-point NURBS circle is in the region tests).

Recommendation and choice: (b). An ellipse is the affine image of a
circle, and rational curves are affine invariant (P&T §4.2), so A7.1's arc
construction mapped by (x, y) ↦ c + x·A + y·B is the ellipse exactly: the
eccentric anomaly range [θ₀, θ₁] (sweep 0 < Δ ≤ 2π) splits into
k = ⌈Δ / (π/2)⌉ equal segments of δ = Δ/k; each has end control points
E(θᵢ), E(θᵢ₊₁), and middle control point
c + (cos θₘ · A + sin θₘ · B)/cos(δ/2), θₘ the segment's middle, with weight
cos(δ/2); knots 0,0,0, 1/k,1/k, …, 1,1,1. A full ellipse is four
segments. Nothing downstream changes, and nothing is approximated.
Analytic ellipse edges and elliptic-cylinder faces can be recognised from
this form when an exchange format needs them (Phase 4).

## 7. Tools

* **Ellipse**: click the centre, the end of the first (major) axis, then
  the extent of the second: the third click is projected onto the line
  through the centre perpendicular to the first axis, so it marks where
  the minor axis ends (on the side clicked, which sets the handedness).
  Refused, keeping the first two clicks, while that extent is within
  1e-7 m of zero. Makes the three points, the ellipse and two
  construction lines, centre to major end and centre to minor end, which
  reference the ellipse's points (no extra constraints, no extra
  freedom). Deleting the ellipse deletes these two lines with it (unless
  another ellipse uses them); deleting a line leaves the ellipse.
* **Elliptical arc**: the same three clicks define the ellipse (shown
  dashed with its semi-axes), then the start and the end: each click is
  taken at its eccentric anomaly on the ellipse (§4.2's exact inverse,
  i.e. along the ray from the centre in the axis frame), so it lands on
  the curve. The arc runs the way the pointer turned from start to end:
  the editor accumulates the change of anomaly as the pointer moves, as
  the centre-point arc tool accumulates angle. Arcs run counter-clockwise
  in the plane, so the ends are stored swapped when the turn's sign and
  the handedness disagree. A start or end snapped to an existing point or
  curve is constrained to it like any snapped click. The arc gets the same
  construction semi-axes.
* Picking uses the drawn polyline. A snap onto an ellipse starts at the
  polyline's nearest point and refines it by Newton's method on
  (E(θ) − p)·E′(θ) = 0 from that point's anomaly; the start is within the
  polyline's chord error, so a few steps give the closest point to
  rounding, and the snapped point is then held on the curve (§4.2). The
  global search Eberly's bisection would provide is not needed here, and
  the solver has its own (§4.2).
* Display samples the kernel's exact rational form (§6) at 16 points per
  control point, about 144 segments for a full ellipse (the same rule as
  splines). A zoom-dependent step is not used: no sketch curve has one
  yet, and this one is below a hundredth of the semi-major axis in chord
  error.

## 8. Editing tools

* Mirror and the patterns copy ellipses and elliptical arcs (§4.6).
* Trim, extend and split: every piece of an ellipse is an elliptical arc on
  the same centre and axis points, so cuts are exact; the kernel finds the
  cut points on the rational form (0007 §6), and each cut's eccentric angle
  comes from its position in the axis frame. As built in §8.1.
* Fillet: needs tangency with a circle (§4.5); refused, naming the
  elliptical arc, until that exists.
* Offset: refused, the offset is not a conic (Farouki & Neff).
* Tangent constraint tool: a line and an ellipse or elliptical arc get
  §4.3's tangency, the centre's side taken from where it is; curves
  meeting at an elliptical arc's end get §4.4's. An ellipse with a circle,
  an arc or another ellipse and no shared end is refused with a message
  (§4.5). Equal, concentric, normal, radius and diameter do not offer
  ellipses (§5).
* Line inference: a line drawn from an elliptical arc's end along its
  tangent is made tangent there, as from a circular arc's.

### 8.1 As built (trim, extend, split)

* **Parameters.** A cut's parameter on an ellipse is its eccentric
  anomaly (§4.2's exact inverse, as the cut point is on the curve)
  turned counter-clockwise in the plane from the first axis end, and on
  an elliptical arc from its start, in [0, 2π). Pieces are then ordered,
  chosen and previewed exactly as on a circle or a circular arc, for
  either handedness. The kernel's own parameter on the rational form runs
  the same way but not evenly, so it is only used to find the points.
* **Own semi-axes do not cut.** The construction lines the tool drew from
  the centre to the axis ends (§7) are part of how the ellipse is drawn
  and touch it only at its axis ends. Counted as cutters they would stop
  every trim and extension at those two points, which nobody drawing a
  line across an ellipse means. Every other curve cuts as usual,
  construction lines included, and the semi-axes of one ellipse cut
  another.
* **Trim.** An ellipse needs two crossings (else refused, naming it); it
  becomes an elliptical arc on its own three points, which takes every
  constraint on the ellipse: they all use the full ellipse (point on
  curve, tangency to a line, images). An elliptical arc's end piece moves
  that end to the cut and holds it on the cutter (0007 §7.4a). A middle
  piece leaves a second arc on the same three points from the far cut to
  the old end, which takes over the old end's constraints; the shared
  points define the ellipse, so, unlike a circular arc's piece, it needs
  no equal-radius constraint. An ellipse nothing else cuts is deleted
  with its semi-axes and points, as any uncut curve.
* **Extend.** An elliptical arc's end runs on around the rest of its
  ellipse (counter-clockwise from its end, clockwise from its start) to
  the first curve met, the arc's own other end and semi-axes not
  counted. An ellipse has no ends to extend, and an arc whose gap is
  under about 2 µm is refused as already whole.
* **Split.** An elliptical arc splits at the pointer's closest point
  (Newton's method from the eccentric anomaly) into two arcs on the same
  three points joined by a coincidence. The two ends at the joint each
  have their own equations putting them on the one ellipse, and with the
  coincidence the normal component of that is said twice: one row short
  of full rank, so the edit was refused as redundant. The solver now
  holds an elliptical arc's end on its ellipse once per coincidence
  class: an end coincident with a point already on the same ellipse (an
  axis end, or another arc's end, the ellipse being the same when its
  three points are, up to coincidence) gets no equations of its own. The
  equations dropped are implied, so degrees of freedom are unchanged; a
  user's point-on-curve constraint at such a joint is still reported
  redundant, as it says again what the arcs say. An ellipse cannot be
  split at one point (no ends): trim it.
* **Tests.** `crates/ui/tests/sketch_ellipse.rs` drives each edit through
  the editor for both handednesses: the pieces kept, their shared points,
  the constraints added and removed, degrees of freedom, no redundancy,
  the hover previews, one undo step, and exact elliptic segment areas
  (ab/2)(Δ − sin Δ) of the regions left (1e-11 relative) and of their
  extrusions; the semi-axes must not stop a trim or an extension. The
  solver test checks the joint rule (no redundancy, DOF 5 + 3, dragging
  keeps the joint on the ellipse) and fails without it. The trim fuzz
  (`crates/ui/tests/sketch_edit.rs`) now also draws random ellipses and
  elliptical arcs of either handedness, directly and through the tools'
  own builders with their semi-axes: in a 2,000-case run per property,
  1,469 ellipse edits committed (374 ellipse trims, 431 arc trims, 456
  splits, 208 extensions), every one solved, without redundancy, with
  nothing unplanned moving.

## 9. Tests

As built (2026-10-12). The solver's tests are in
`crates/sketch-solver/tests/ellipse.rs`, the kernel's in
`crates/kernel/tests/ellipse.rs` and the editor's in
`crates/ui/tests/sketch_ellipse.rs`; each bullet below is covered.

Measured solver robustness, because over-determination makes least
squares lose here more often than elsewhere: random sketches of motifs
with ellipses (points on them, tangent lines, arcs joined to other
curves) whose point coordinates are perturbed by up to a given fraction
of the sketch size, then solved. Failures, all reported as failures with
a diagnosis (never a wrong answer, in about 1.1 million cases, each
success checked by independent closed forms):

| perturbation | failures |
|---|---|
| ≤ 0.01 | 1 in 300,000 |
| ≤ 0.02 | 1 in 300,000 |
| ≤ 0.05 | about 8 in 100,000 |

The same generator without ellipses fails 0 in 100,000 at ≤ 0.05. The
failures are local minima of an over-determined subsystem and tangencies
started near an ellipse's vertex, where the contact point is a fold of
the equations (the generator avoids contacts within 0.3 rad of a vertex
for that reason). The proptest in CI runs ≤ 0.01; larger perturbations
are tested for "never lies" only.

Planned when this note was written:

* Solver: DOF 5 and 7; a point on an ellipse removes one; tangency to a
  line checked against an independent closed form (the distance from the
  centre equals the support value), on both sides; a fully constrained
  ellipse through its axis lines (DOF 0, no redundancy); images add no
  freedom and no redundancy; refusals name the constraint; a fuzz of
  random ellipses, arcs and constraints that either solve and pass an
  independent check or fail with a diagnosis.
* Kernel and features: area πab and elliptic sector (ab/2)(θ₁ − θ₀) minus
  its triangle, to 1e-12 relative; extrude volume πab·h; a half-ellipse
  revolved about its major axis is a prolate spheroid of volume
  (4/3)πab², about its minor axis an oblate one, (4/3)πa²b.
* Editor: both tools by clicking, refusals, handedness both ways (a
  mirrored arc), snapping onto the curve, dragging an axis end.
