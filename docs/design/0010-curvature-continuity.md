# Design note 0010 — Curvature continuity (G2) in sketches

Status: accepted (2026-10-12). Extends design notes 0005 (solver), 0007
(editing tools), 0008 (ellipses) and 0009 (conics).

## 1. Scope

A sketch constraint that makes two curves meeting at a shared end point
curvature continuous there: one tangent line and one curvature, so the
joint is smooth to second order (G2). Its data, the solver's equations,
which pairs of curves it accepts, how the editing tools treat it, and the
tests. Publicly documented sketchers offer such a constraint between a
spline and another curve; the details below are ours.

## 2. Literature

* M. P. do Carmo, *Differential Geometry of Curves and Surfaces*,
  Prentice-Hall, 1976, §1-5: the signed curvature of a plane curve,
  κ = (x′y″ − y′x″)/(x′² + y′²)^{3/2}; it changes sign with the direction
  of travel and is unchanged by any other reparametrization.
* B. A. Barsky & T. D. DeRose, "Geometric continuity of parametric curves:
  three equivalent characterizations", *IEEE CG&A* 9(6), 1989: G2 at a
  joint is a common unit tangent and a common curvature vector.
* G. Farin, *Curves and Surfaces for CAGD*, 5th ed., Morgan Kaufmann,
  2002, the chapters on the differential geometry of curves and on
  rational Bézier curves: curvature of a (rational) Bézier curve at its
  end point, κ(0) = ((n − 1)/n)·(w₀w₂/w₁²)·h/a², with a the length of the
  first leg and h the distance of the third control point from it.
* Piegl & Tiller, *The NURBS Book*, 2nd ed., A2.3 (derivatives of the
  basis functions), from which the solver's spline weights give C′ and C″
  at the ends.
* Design note 0005 §4 (tangency at a shared point, scaling), 0008 §3–4
  (ellipse parametrization, the parameter-free tangent), 0009 §3–4 (the
  conic's weights and shape unknown).

## 3. Data

```text
Constraint::Curvature { a, b, point }
```

`a` and `b` are curves, `point` the point entity where they meet: an end
point of both, or one of two coincident end points (as
`Tangent { kind: AtPoint }`). It is a relation (no value).

Accepted pairs: at least one curve is a spline (control-point or
interpolated, degree ≥ 2); the other is a line, an arc, an elliptical arc,
a conic arc or a spline (degree ≥ 2). Refused, with the reason:

* no shared end point ("Curvature continuity needs the curves to share an
  end point"), or a circle or full ellipse, which has no end;
* neither curve is a spline: two lines always agree (zero curvature), a
  line never matches an arc, ellipse or conic, and two arcs agree only by
  being one circle. Those cases are refusals rather than conflicts the
  solver finds, because nobody means them;
* a degree-1 spline: it is straight at its ends, so the constraint would be
  a tangency or a conflict;
* both ends of one curve at the point (a closed spline): the end is
  ambiguous, as for tangency.

## 4. Equations

Each side contributes, at the joint, two vectors D₁, D₂ such that D₁ is the
tangent pointing **away** from the joint along the curve and

  κ = cross(D₁, D₂)/‖D₁‖³

is the signed curvature for travel away from the joint. Any pair
(λC′, λ²C″ + μC′) of a parametrization C gives the same κ (do Carmo), so
D₁, D₂ need not be true derivatives:

| side | D₁ | D₂ |
|---|---|---|
| line from p to q | q − p | 0 |
| arc on centre c, joint at its start p | perp(p − c) | c − p |
| arc, joint at its end p | −perp(p − c) | c − p |
| elliptical arc (centre c), joint p at its start / end | ±σ·T(p) | c − p |
| conic P₀P₁P₂, weight w = e^u, joint at P₀ | P₁ − P₀ | e^{−2u}(P₂ − P₁)/2 |
| conic, joint at P₂ | P₁ − P₂ | e^{−2u}(P₀ − P₁)/2 |
| spline, joint at its start | C′(0) | C″(0) |
| spline, joint at its end | −C′(1) | C″(1) |

perp rotates a quarter turn counter-clockwise; arcs and elliptical arcs run
counter-clockwise from start to end. For the ellipse T is the
parameter-free tangent of design note 0008 §4.4 (E′(θ) at the point's θ,
whose second derivative is E″ = c − E), and σ = sign cross(A, B) turns
the θ direction into the counter-clockwise one. For the conic,
cross(D₁, D₂)/‖D₁‖³ is Farin's end curvature with n = 2, w₀ = w₂ = 1,
w₁ = w (a quarter circle of radius R as a conic, w = 1/√2, gives 1/R).
For a spline C′ and C″ at the ends are fixed linear combinations of its
referenced points (Piegl & Tiller A2.3 through the spline map of design
note 0005 §3), so D₁ and D₂ are linear.

At a smooth joint the two away tangents are opposite, and curvature is
continuous for travel through the joint when κ_a = −κ_b (travel through
the joint is away from it on one side and towards it on the other, which
reverses that side's sign). Two rows:

* **Smooth joint.** The angle θ from D₁ of `a` to D₁ of `b` is π:

    r₁ = sin(θ − π),  required cos(θ − π) > 0.

  This is the angle row of an angle dimension between two directions
  (design note 0005 §4) with the value π, and it is solved the same way: a
  solution with cos(θ − π) ≤ 0 is the wrong branch and is not accepted,
  and a start far from π (the curves meet at a corner) is reached by
  continuation, the value moved from the start angle θ₀ to π in stages,
  each solved from the last.
* **Curvature.** r₂ = κ_a + κ_b.

The tangency row of design note 0005 §4 (the cross product of the two
tangents) is deliberately not used for the first row. It also holds at a
cusp, θ = 0, where the curves leave the joint together, and a first
version paired it with a curvature row valid at both, κ_a − (û_a · û_b)κ_b
(û = D₁/‖D₁‖). Newton then found cusps from ordinary starts where the
curves met at a corner (the fuzz of §6 caught them): a 'solution' that
satisfies every row and is never what is meant. With the angle row a cusp
is not a solution at all: the branch check refuses it, and the solver
reports a failure rather than return one.

Both residuals have analytic gradients (the chain rule through θ, κ and,
for the ellipse and conic, the Jacobian of T and the factor e^{−2u}),
checked against central differences. In the solver's scaled units (lengths
divided by the sketch size L) r₂ is a curvature times L, a length-like
quantity: it converges to the length tolerance, which bounds the position
error the mismatch causes over the sketch's size by the same 1e-12 · L; r₁
converges to the angular tolerance. A D₁ shorter than the degeneracy limit
(a spline whose first two control points meet) is reported as degenerate,
as for the tangency row.

DOF: the constraint removes two freedoms (tangent direction, curvature).
A tangency on the same two curves at the same point is redundant with it
and is reported as such (the sketcher replaces it, §5).

### 4.1 Which curve gives

The solver's steps are weighted minimum-norm (design note 0005 §5): of
all the changes that fix the linearized residuals, the one that moves the
geometry least. With every point weighing the same, that is often not the
change anyone means. When the spline leaves the joint bending the other
way from the curve it meets, the least change flattens both towards zero
curvature: in the base joints of §6 an arc's radius tripled, and in 10 %
of random elliptical joints the ratio of the ellipse's axes fell a
hundredfold. The
constraint is between a spline and a curve, and it is the spline that is
meant to adapt; so:

* **The spline gives first.** The points of a spline that a curvature
  constraint holds weigh 1e-3 in the step metric (a point shared with
  another curve keeps the weight of that curve's point, as the weights of
  coincident points add). The other curve then moves only for what the
  spline cannot absorb: in the base joints of §6 its curvature at the
  joint changes by under 1 % for arcs, 5 % for elliptical arcs and 4 %
  for conics (whose shape ρ is a free number). Between two splines both
  give alike, and their curvatures meet in the middle.
* **A trust region on light points.** A light unknown's step is larger
  than a unit-weight one's by about the inverse square root of its weight
  (≈ 30×), so its linearization is trusted over a correspondingly smaller
  region: the line search starts at the fraction of the step that moves
  no light point by more than a quarter of the sketch size (Nocedal &
  Wright §4). Without it the first step could throw a spline's points
  across the sketch into a loop, from which the solve settled on a
  flattened ellipse.
* **Then everything gives.** If the solve with the spline first finds no
  solution, the sketch is solved again with every point weighing the
  same, and that result is used if it is one. Measured (§6), the first
  attempt fails on 1.5 % of perturbed joints (mostly an interpolated
  spline meeting a sharply curved conic), the second rescues all but
  0.27 %; in the rescued joints the other curve may change shape (a conic
  flattening), which is a solution of the constraints as given and can be
  prevented by dimensioning that curve.

## 5. In the sketcher

* The constraint tool takes two curves; the point is their shared end
  (as for tangency). A tangency already holding the same two curves at the
  same point is replaced by the curvature constraint, which implies it;
  the message says so.
* The glyph is "C", drawn next to both curves like tangency's "T".
* Editing tools treat it like a tangency at a shared point: trim and split
  move it with the end it holds, deleting either curve or the point drops
  it, transform keeps it (a similarity keeps curvature continuity), mirror
  and patterns carry nothing over (their images are held by their own
  constraints).
* The kernel needs nothing: the constraint only positions geometry.

## 6. Tests

Solver (`crates/sketch-solver/tests/curvature.rs`):

* A spline (control-point and interpolated) joined to each accepted kind
  of curve at either end (line, arc at start and end, elliptical arc at
  start and end and with its minor axis on either side, conic at either
  end, quadratic, cubic and interpolated splines), from a start that is
  not G2: solved, every constraint holding, a smooth joint (away tangents
  opposite to 1e-12) and κ_a + κ_b = 0 to 1e-9, checked against
  curvatures computed independently in the test from closed forms (lines
  0, arcs ±1/r, elliptical arcs |A×B|/‖E′(θ)‖³, conics by the quotient
  rule on the rational form, splines by Piegl & Tiller eq. 3.8 applied
  twice); a curve that is not a spline keeps its curvature at the joint
  to within 6 % (measured up to 4.9 %).
* DOF drops by 2; a tangency on the same joint is reported redundant; a
  fixed spline against an arc of fixed radius is a conflict with the
  input returned untouched; every refusal of §3; the same joints solved
  alike at 1e-4 and 1e3 scale.
* A cusp is never the answer: a spline leaving back along the line it
  meets is turned round or the solve fails.
* Corners of 10° to 179° between a line and a spline (both kinds) all
  solve to a smooth joint.
* Property test: 256 random joints, every point moved by up to 15 % of
  the size: solved with every constraint holding, or failed with the
  input untouched.
* Measured statistics, deterministic: every kind of joint with every
  point moved by up to 15 % of the sketch size, solved with and without
  the constraint. Over 1,000 joints of each of the 11 kinds: 30 of
  11,000 fail (0.27 %; 24 of them a spline meeting an interpolated
  spline), all
  reported, none wrong; where the other curve is not a spline its
  curvature at the joint ends within 25 % of where it ends without the
  constraint in 98.3 % of joints. Before §4.1 the same joints failed 111
  times (1.0 %), and the ratio of an elliptical arc's axes fell a
  hundredfold in 311 of 3,000; with the spline first and no trust region
  it fell so in 25; with the trust region, in none. An elliptical arc
  drawn with its minor axis on the other side (the same curve) failed 255
  times in 2,000 against 13 before §4.1, and 9 against 9 after. The test
  runs 100 of each kind and bounds failures at 1 % and kept curves at
  95 % (measured 4 of 1,100 and 99.1 %).
* Gradients: the analytic gradient of both rows against central
  differences for random configurations of every side kind.

Editor (`crates/ui/tests/sketch_curvature.rs`): the tool joins a spline
drawn bending the wrong way to an arc (G2 checked from the solved control
points, the arc's radius within 1 %, one undo step), its refusals (no
spline, a circle, curves that do not meet, a point), replacing a
tangency (one freedom more, nothing redundant), and trim and split
handing the constraint on to the piece at the joint.
