# Design note 0003 — Mass properties by boundary integration

Volume, area, first and second moments are computed on the exact geometry,
not on a mesh.

* Divergence theorem: `∫_V g dV = ∮ F·n dA` with `div F = g`. We use
  `F = p/3` (volume), `F = (x²/2, 0, 0)` etc. (first moments),
  `F = (x³/3, 0, 0)` (pure second moments) and `F = (x²y/2, 0, 0)` etc.
  (products).
* On a face with parameterization `S(u, v)`, `n dA = σ (S_u × S_v) du dv`.
  Green's theorem turns the parameter-domain integral into a boundary
  integral: `∬ f du dv = ∮ (∫_{u₀}^{u} f(s, v) ds) dv` (inner along u;
  used for revolution surfaces and planes) or `−∮ (∫_{v₀}^{v} f ds) du`
  (inner along v; used for linear extrusions, so the inner integral never
  crosses a NURBS knot).
* The loop traversal direction carries the face sense σ for the flux terms.
* The boundary integral runs along the exact edge curves (Gauss–Legendre,
  8 points per sub-interval, adaptive: below); the parameter velocity is
  the least-squares solution of `[S_u S_v]·(u′, v′) = C′(t)`. Singular
  points contribute explicit parameter-space connector segments.
* Inner integrals use 16-point Gauss–Legendre.
* Moments are accumulated about the vertex centroid for conditioning and
  shifted with the parallel-axis theorem.

## Adaptive boundary quadrature (2026-10-12)

The first version split each edge where a polyline would deviate from it
by the requested accuracy times the body's size and applied the 8-point
rule to every segment. The cost grew without bound as the accuracy
tightened (a 15-letter text extruded at 1e-9 took 128 s), while the rule
was already exact to rounding on a few pieces of smooth geometry.

Now each edge is split into its smooth pieces only (a line whole, a circle
in pieces of at most π/4, a NURBS curve at its distinct knots), and each
piece is integrated by recursive bisection in the manner of Gander and
Gautschi (W. Gander, W. Gautschi, "Adaptive Quadrature — Revisited", BIT
40(1), 2000, pp. 84–101): the rule on a piece is compared with the rule on
its two halves, and the halves are accepted when, for each of the eleven
integrals, the difference is at most

* the requested relative accuracy times the rule's integral of the
  integrand's magnitude over the piece (so the total error is bounded by
  that accuracy times the boundary integral of the magnitude; for the
  volume of a solid star-shaped about the vertex centroid that is the
  volume itself), or
* a rounding floor of 1e-13 times the body's size to the integral's
  dimension, shared among the coedges by parameter length (for components
  that vanish on a piece up to rounding, whose relative difference is
  noise).

Accuracies finer than 1e-13 mean 1e-13. A piece that has not converged
after 30 bisections makes the whole computation fail with an error; the
kernel never returns a value it cannot vouch for. The difference between
the rules estimates the error of the coarser one, so for smooth integrands
the accepted halves are far more accurate than asked (a strongly rational
conic asked for 1e-3 comes out within 2e-6).

Two cases found while building it are regression tests:

* A conic arc of weight 134 puts almost all of its length within 1/w of
  its middle parameter. A draft of this change stopped refining when bisection
  did not halve the difference ("stalled at the noise floor"), which
  before the asymptotic regime returned a volume 21 % too small. Stalls
  are no longer a reason to stop; only the relative and rounding tests
  are (`random_segments_have_their_exact_areas`, and a unit test of a
  peak of width 1e-5).
* Splitting a spline one rounding step from an existing knot leaves a
  knot span of length 5.6e-17, whose polynomial piece is numerically
  meaningless: sampled inside it, the integrand varies by 10 %, and
  bisection never converged. Knot spans shorter than 1e-12 of an edge's
  parameter range are merged into their neighbours (their contribution is
  below rounding) (`spline_crossing_a_line_several_times`).

## Tolerant edges (2026-10-12)

Each face integral follows an edge as projected onto that face's own
surface (the parameter velocity above). Where an edge lies on its two
surfaces only within the tolerance, the two faces trace two slightly
different curves, and the closed surface the divergence theorem needs is
open along a thin strip between them. That happens wherever profile end
points up to ten tolerances apart were merged into one vertex: the edge
curves end up to five tolerances from the vertex, so the edge swept from
the vertex lies off one of its two side faces by as much. The flux through
an open surface depends on the integrands' origin, so the volume was wrong
by about the strip's area times its distance from the origin: 1.6e-5 of
the volume of a 3 mm sliver with a 0.47 µm gap, about 1e-7 on ordinary
bodies, far outside the requested accuracy and invisible to every check
that compared the kernel with itself.

Following the tolerant-modelling view that an edge is a tube of tolerance
radius about its curve, inside which its faces may meet anywhere (D. J.
Jackson, "Boundary representation modelling with local tolerances", Proc.
3rd ACM Symposium on Solid Modeling and Applications, 1995), the strip
between the two projections is added to the boundary. For every edge used
by two different faces, with `w(t)` the vector from its projection onto
the first face's surface to its projection onto the second's, the strip is
ruled between them and its flux taken along the edge curve, `n dA ≈
(C′ × w) dt`, with the loop direction's sign; that differs from the exact
ruled strip by terms of the order of the width squared. The area integral
leaves the strip out (it is no part of the surface). Edges on both
surfaces to rounding (width at most 1e-14 of the distance from the origin
at nine points) are skipped, so exact bodies cost nothing extra; the others
use the same adaptive rule as the faces.

Regression tests: a square whose side stops (3e-7, 4e-7) short of a
corner now has the volume of the region its curves bound with the closing
chord to 4e-14 of its volume (it was 5e-8 off) (`square_with_a_gap_at_a_corner`), and
the three regions of a spline sketch with merged end points extrude to
area × height within 1e-12 (`spline_regions_with_gaps_at_their_corners`).

The golden corpus checks cuboids, cylinders (including inertia tensors),
cones, spheres, tori, partial revolves, plates with holes, spline
extrusions and spline revolves (against Pappus' theorem) to relative error
≤ 1e-9 (1e-6 where the reference itself is a 20,000-segment polyline).
