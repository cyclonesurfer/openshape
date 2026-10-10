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
  8 points per sub-interval); the parameter velocity is the least-squares
  solution of `[S_u S_v]·(u′, v′) = C′(t)`. Singular points contribute
  explicit parameter-space connector segments.
* Inner integrals use 16-point Gauss–Legendre.
* Moments are accumulated about the vertex centroid for conditioning and
  shifted with the parallel-axis theorem.

The golden corpus checks cuboids, cylinders (including inertia tensors),
cones, spheres, tori, partial revolves, plates with holes, spline
extrusions and spline revolves (against Pappus' theorem) to relative error
≤ 1e-9 (1e-6 where the reference itself is a 20,000-segment polyline).
