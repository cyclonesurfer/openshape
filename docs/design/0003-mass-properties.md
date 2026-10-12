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

The golden corpus checks cuboids, cylinders (including inertia tensors),
cones, spheres, tori, partial revolves, plates with holes, spline
extrusions and spline revolves (against Pappus' theorem) to relative error
≤ 1e-9 (1e-6 where the reference itself is a 20,000-segment polyline).
