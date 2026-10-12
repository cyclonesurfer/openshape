# ADR 0016 — Smooth closed splines

Status: accepted (Phase 2)

## Context

A spline drawn back to its first point closes there with a corner: both
sketch splines are clamped (design note 0007 §4), so the curve leaves its
first point along the first leg of its polygon and arrives along the last
one. Closed outlines (cams, knobs, lettering, organic profiles) need a
curve that closes smoothly, with no seam. The roadmap lists "smooth
closed" splines. Four questions:

1. **What closing a spline does**: a corner, a smooth closure, or a choice.
2. **What the document stores.**
3. **What the solver and the kernel see.**
4. **How the through-point spline closes**: its parameters and knots.

## Options

### 1. Closing

* **A. Keep the corner**; smooth closure through a separate option.
* **B. Closing is smooth**; a corner is still possible by ending the
  spline next to its start and joining the two points with a coincident
  constraint (that is the clamped closed form documents already hold).
* **C. Ask every time.**

### 2. Stored form

* **A. A `periodic` flag on both spline kinds.** The point list still ends
  on its first point, as a corner-closed spline's does, so "closed" keeps
  one meaning in the document; the flag says the closure is smooth.
  Documents written before the flag read as `false` (serde default), so
  existing corner closures are unchanged.
* **B. Two new geometry kinds.** Every match over geometry grows two arms
  that mostly repeat the open kinds' code.
* **C. Drop the repeated point** for periodic splines. The list then has
  no closing entry, which every consumer that walks a spline's points
  (mirror pairs, patterns, transform, display of the control polygon)
  would have to special-case.

### 3. Solver and kernel

* **A. Periodic knot vectors throughout**: the kernel and the region code
  learn unclamped, periodic B-splines.
* **B. An exact clamped equivalent.** A periodic B-spline restricted to
  one period is a clamped B-spline: inserting the domain's end knots until
  they have multiplicity `p` (Boehm insertion, Piegl & Tiller, *The NURBS
  Book*, 2nd ed., §5.2, algorithm A5.1) turns the wrapped control polygon
  into the control polygon of a clamped curve with the same shape on the
  same uniform interior knots (P&T §12.2 for the wrapped construction).
  Each clamped control point is a fixed linear combination of the
  sketch's points, which is exactly how the solver already treats an
  interpolated spline (a linear map from referenced points to the curve,
  design note 0005 §3), so the solver needs no new unknowns and the kernel
  receives a closed clamped spline, which it already handles.

### 4. Through-point closure

A closed through-point spline passes through `n` distinct points at
parameters `0 = t₀ < … < t_{n−1} < t_n = 1` (chord length over the closed
polygon, as the open tool's), returning to the first at `t_n`.

* **A. Knots by averaging**, as the open spline (P&T eq. 9.8). Averaging is
  defined for clamped ends; for a periodic curve it leaves the closing
  span unlike the others.
* **B. The periodic cubic spline with knots at the parameters**, the
  classical closed interpolating cubic (C. de Boor, *A Practical Guide to
  Splines*, rev. ed., Springer 2001, ch. IV and XIII; the periodic
  interpolation problem with knots at the data sites has a unique
  solution for any increasing parameters). As a clamped spline on the
  knots `0⁴, t₁ … t_{n−1}, 1⁴` it has `n + 3` control points, fixed by the
  `n + 1` interpolation conditions and by equal first and second
  derivatives at `t = 0` and `t = 1`: a square, nonsingular system, solved
  once per definition like the open spline's.

## Decision

1B, 2A, 3B, 4B.

* Clicking the first point again with the spline tools (after at least
  three distinct points) closes the spline smoothly. Existing documents
  keep their corners.
* A periodic control-point spline of degree `p` needs at least `p + 1`
  distinct control points and degree 2 to 5 (a degree-1 periodic spline is
  a polygon, not a smooth curve). The tool makes degree `min(3, n − 1)`
  for `n` distinct points, as for open splines.
* A periodic through-point spline is cubic and needs at least three
  distinct points (through two, the closed cubic collapses onto the chord).
* The first and last clamped control points are the same combination of
  the same point, so the curve closes exactly, not to rounding.
* A smooth closed spline has no end points: tangency and curvature
  continuity at its first point are refused with the reason. Points on it,
  dimensions to its points, mirror, patterns and transform work as for
  open splines and keep the flag.

## Consequences

* The kernel, region code and STL path see a closed clamped spline. Its
  second-to-last control point is set from the first two (P&T eq. 3.7),
  so its first derivatives at the two ends are equal exactly, not just to
  rounding; the construction already makes them equal mathematically, so
  this moves that point by rounding only.
* A profile loop cuts a closed curve in two, so an extruded or revolved
  closed spline has two faces on its side. The kernel's display test
  classifies an edge where swept faces from two profile pieces meet with
  equal first derivatives as a seam, not a crease, so the side shows as
  one smooth wall; a tangent joint between different curves (unequal
  speeds) and a corner closure still show their edge. A split a rounding
  step from a knot splits at the knot, so no piece ends on a sliver span
  whose end derivative would be meaningless.
* The extra control points of the clamped form are internal; the user sees
  and drags the closed control polygon (control-point spline) or the
  through points.
* A spline closed with a corner stays expressible (coincident ends), and
  documents that hold one are read unchanged.
