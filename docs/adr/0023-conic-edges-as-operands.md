# ADR 0023 — Conic edges as boolean operands: follow the stored curve

Status: accepted, 2026-10-15.

## Context

Booleans of planes across cones and cylinders (design note 0012 §14,
§15) leave ellipses, parabolas and branches of hyperbolas as edges of
curved faces. The kernel stores them as exact rational quadratic NURBS
curves (`Curve3::Nurbs`), since it has no curve kind for an ellipse or a
conic in space. A second boolean on such a body was refused: the stages
follow an edge's course in a face's chart with a proven band (its
distance from its chords, §12.1), and had bands for lines, circles and
the conic tracks they compute themselves, not for a NURBS curve on a
curved face. STATUS.md planned to "recognise the rational quadratic
edges the booleans write as their conics". This record changes that
plan, so it is written down.

## Options

* **A. Recognise the conic from the curve** (its control points and
  weights: Lee 1987; Farin, *Curves and Surfaces for CAGD*, ch. 12), and
  follow the recognised track. Ill-conditioned where it matters: on a
  short or flat piece the shape factor's difference from a parabola's
  shrinks with the square of the piece's turning, a nearly circular
  ellipse has no defined axis, and a parabola's centre runs to infinity.
  A recognition wrong by more than a band would let the split misplace a
  crossing while the result's edge, still the right NURBS, passes
  validation: a silent error. Making it safe needs a proven bound on the
  distance between the NURBS and the recognised track, which needs the
  same derivative bounds as option C.
* **A′. Recognise it with the face's help** (the axis from the cylinder
  or cone, the plane from a neighbouring face), then verify the track
  against the curve. Better conditioned, but two stored forms (an
  ellipse's half and a conic about its vertex) and data-driven choices
  of sign, wrap and branch, each a way to follow the wrong stretch.
* **B. Store exact conics on edges** (`Curve3::Ellipse`,
  `Curve3::Conic`), parameterised as the tracks are. Exact and tight
  bands, and exchange formats (STEP) will want exact conics. But every
  place a curve is evaluated, projected, boxed, sampled and meshed gains
  a case (about twenty in the kernel), and an ellipse needs its own
  convention for parameter ranges. A side record next to the NURBS
  instead gives two sources of truth that later edits (transforms, knot
  insertion) could make disagree.
* **C. Follow the stored NURBS** in cylinder and cone charts, with a
  band proven from its own derivatives: the image of a curve in an
  isometric chart bends by at most its acceleration (do Carmo §4-4),
  bounded over each Bézier piece by the Bernstein coefficients of
  `X'' = Q / W³`, with terms for its rounding off the surface from the
  surface's implicit function along it (design note 0012 §16.2).

## Decision

C. The NURBS stays the one record of the edge, nothing is decided from
the data, and the band is a mathematical property that tests can check
against dense samples. It also covers any rational curve of degree at
most five on a cylinder or a cone, not only the conics booleans write.
B is left for exchange formats, where exact conics are needed in any
case.

## Consequences

* A band from `|X''|` is looser than a track's where the curve's speed
  varies: measured against dense samples, 4.7 times the true distance
  from the chord at the median and up to `7.4 × 10⁵` far along a
  hyperbola's branch, which costs about ten more halvings of a chord
  there; no answer changes.
* Every bound carries an allowance for its own rounding (design note
  0012 §16.2), checked against double-double references, and with it
  its floor, the part of the bound that rounding makes and halving the
  chord does not narrow. Halving stops at four floors; where that is
  over a sixteenth of the tolerance and the chord could hide a contact,
  the boolean is refused as too far from the origin (from about 36 km
  at the default tolerance for an evenly weighted conic on a cylinder).
* Two curved segments between the same points are one when points of
  each, a quarter, half and three quarters of the way along it by its
  parameter and by its tangent's turning, lie within the tolerance of
  the other's stretch, both ways (§16.3). Matching by the middles of
  their parameters, as for circles, failed: a NURBS's parameter and its
  track's differ. No single point would do, as two different conics
  can pass through the same ends and middle.
* NURBS edges on spheres stay refused: no conic lies on a sphere.
