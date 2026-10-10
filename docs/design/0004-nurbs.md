# Design note 0004 — NURBS library

`crates/nurbs` implements the algorithms of Piegl & Tiller, *The NURBS
Book* (2nd ed.), with algorithm numbers cited in the code:

* Knot spans and basis functions: A2.1, A2.2, A2.3 (derivatives).
* Curve/surface points and derivatives in homogeneous form, rational
  derivatives via A4.2 (curves) and A4.4 (surfaces).
* Knot insertion A5.1; splitting by repeated insertion.
* Exact circular arcs A7.1; surfaces of revolution A8.1.
* Global interpolation A9.1 with chord-length parameters.
* Point inversion (closest point) by dense sampling followed by Newton
  iteration on the orthogonality conditions (§6.1, eqs. 6.3–6.8), clamped to
  the domain and never worse than the best sample.

Surface/surface intersection (marching with Bézier clipping seeds, following
Patrikalakis & Maekawa, *Shape Interrogation for CAD/CAM*, 2002) is the next
step and gets its own design note before Phase 3.
