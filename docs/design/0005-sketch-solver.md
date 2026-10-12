# Design note 0005 — Sketch constraint solver

Crate: `crates/sketch-solver` (`openshape-sketch-solver`), using
`crates/linalg` (ADR 0010). The data types (entities, constraints) are the
contract shared with `core-model` (ADR 0011) and are defined in
`src/model.rs`.

## 1. Requirements

* Solve 2D sketches of points, lines, circles, arcs and B-splines under
  geometric constraints and dimensions, to 1e-12 of the sketch size.
* Prefer the solution nearest the current geometry: changing a dimension or
  dragging a point moves as little as possible, and never jumps to another
  branch (a flipped arc, a mirrored triangle) when a nearby solution exists.
* Report, per sketch: degrees of freedom left; which entities are still free
  (under-constrained); which constraints are redundant (implied by others) or
  conflicting (unsatisfiable together); and a clear failure when the solver
  cannot converge. Never return geometry that violates a constraint while
  claiming success.
* Deterministic: the same input gives bit-identical output (ADR 0010).
* Interactive: a 50-entity sketch solves (and re-solves while dragging) well
  within a 60 Hz frame; a 300-entity sketch within 50 ms.

## 2. Literature

* Geometric constraint solving surveys: Bouma, Fudos, Hoffmann, Cai & Paige,
  "Geometric constraint solver", *CAD* 27(6), 1995; Hoffmann & Joan-Arinyo,
  "A brief on constraint solving", *CAD&A* 2(5), 2005; Jermann, Trombettoni,
  Neveu & Mathis, "Decomposition of geometric constraint systems: a survey",
  *IJCGA* 16, 2006.
* Graph decomposition: Owen, "Algebraic solution for geometry from
  dimensional constraints", ACM SMA 1991; Hoffmann, Lomonosov & Sitharam,
  "Decomposition plans for geometric constraint systems", *J. Symbolic
  Computation* 31, 2001.
* Numerical method: Madsen, Nielsen & Tingleff, *Methods for Non-Linear
  Least Squares Problems*, DTU 2004 (Levenberg–Marquardt, damping updates);
  Nocedal & Wright, *Numerical Optimization*, 2nd ed., ch. 10–11
  (Gauss–Newton, minimum-norm steps for underdetermined systems,
  line search); Golub & Van Loan §5.4–5.5 (rank-deficient least squares).
* Rank and redundancy: the numerical-rank approach of Kramer's and
  SolveSpace's solvers (Westhues, SolveSpace, `system.cpp`, GPL-3.0: studied
  as prior art for behaviour only; no code is copied).
* B-splines: Piegl & Tiller, *The NURBS Book*, A2.2, A2.3, A9.1.

## 3. Unknowns

Every point entity contributes (x, y); every circle its radius. Lines and
arcs reference point entities and add no unknowns. A control-point spline
references its control points. An interpolated spline references its
through-points and owns hidden control points (one per through-point, cubic,
knots by averaging the stored parameters, P&T eq. 9.8), determined by
interpolation equations C(tᵢ) = Pᵢ; these are linear in all unknowns and
remove exactly the hidden freedoms, so DOF counts are unaffected. A
point-on-spline constraint owns one hidden curve parameter t, clamped to the
domain after each step.

The arc's radius is |start − centre|; the solver adds the implicit equation
|end − centre| = |start − centre|.

**Scaling.** Let L be the larger of the diagonal of the sketch's bounding
box and the largest length dimension (and 1e-6 m if both are zero).
Unknowns are expressed in units of L; length residuals are divided by L,
angle residuals are dimensionless (sines or radians). This makes
tolerances, rank decisions and convergence independent of the model's
scale (golden tests run the same sketches at 1e-4 and 1e3 scale).

## 4. Equations

All residuals are smooth near their solutions. For lines, d = end − start.

| constraint | residuals |
|---|---|
| Coincident(p, q) | pₓ − qₓ, p_y − q_y |
| PointOnCurve(p, line) | cross(d, p − start) / ‖d‖ |
| PointOnCurve(p, circle or arc) | ‖p − c‖ − r |
| PointOnCurve(p, spline) | C(t) − p (2 rows, hidden t) |
| Horizontal(line) / Vertical(line) | start_y − end_y / startₓ − endₓ |
| HorizontalPoints / VerticalPoints | p_y − q_y / pₓ − qₓ |
| Parallel(a, b) | cross(d_a, d_b) / (‖d_a‖ ‖d_b‖) |
| Perpendicular(a, b) | dot(d_a, d_b) / (‖d_a‖ ‖d_b‖) |
| Equal(line, line) | ‖d_a‖ − ‖d_b‖ |
| Equal(circular, circular) | r_a − r_b |
| Midpoint(p, line) | p − (start + end)/2 (2 rows) |
| Concentric(a, b) | c_a − c_b (2 rows) |
| Symmetric(p, q, axis) | cross(d, (p+q)/2 − start)/‖d‖, dot(d, q − p)/‖d‖ |
| Tangent line–circle | cross(d, c − start)/‖d‖ − s·r, s = +1 if the centre is on the left |
| Tangent circle–circle | ‖c_a − c_b‖ − (r_a + r_b) (external) or ‖c_a − c_b‖ − \|r_a − r_b\| (internal) |
| Tangent at shared point p | line–circular: dot(d, p − c)/(‖d‖ r); circular–circular: cross(p − c_a, p − c_b)/(r_a r_b); spline end: cross(C′(end), other tangent) normalized |
| Normal(line, circular) | the centre lies on the line (as PointOnCurve) |
| Fix(p, at) | p − at (2 rows) |
| Distance(p, q) aligned | ‖p − q‖ − v |
| Distance(p, q) horizontal / vertical | \|qₓ − pₓ\| − v / \|q_y − p_y\| − v |
| Distance(p, line) | \|cross(d, p − start)\|/‖d‖ − v |
| Distance(line, line) (parallel) | distance from a's start to b, as above |
| Length(line) | ‖d‖ − v |
| Angle(a, b) | sin(θ − v) where θ is the signed angle from d_a to d_b; after convergence cos(θ − v) > 0 is checked, so the antipodal branch is rejected |
| Radius / Diameter | r − v / 2r − v |

Dimension values must be positive (lengths, radii) or in (−π, π]
(angles); zero lengths are rejected with "use a coincident constraint".
Entities referenced by the wrong kind (a radius on a line) are rejected
before solving with an error naming the constraint.

Degenerate geometry (a zero-length line under a parallel constraint, a
zero-radius circle) is reported as `Degenerate` instead of producing NaN.

## 5. Algorithm

1. **Validate** references and values.
2. **Decompose** the equation–unknown incidence graph into connected
   components; each is solved independently (Hoffmann et al. decompose
   further into rigid clusters; that refinement is future work and does not
   change results, only speed).
3. **Solve each component** by damped Gauss–Newton with minimum-norm steps:
   * J Δ = −r is solved in the weighted minimum-norm least-squares sense:
     minimize ‖W^{1/2} Δ‖ among the least-squares solutions, using QR with
     column pivoting of J W^{−1/2} and a complete orthogonal decomposition.
     Weights are 1, except for dragged points (§7). Images of derived
     constraints (design note 0007 §9) move first: when moving them alone
     solves the linearized system (to 1e-9 relative, or to a tenth of the
     tightest convergence tolerance), that is the step, the limit of the
     weighted step as their weight goes to zero. So a change only images
     need absorb (an offset's distance) leaves every other point exactly
     where it was. A dragged image moves last, like any dragged point.
   * Step acceptance by backtracking line search on ½‖r‖² (Armijo, c = 1e-4,
     halving down to 2⁻²⁰); if the line search fails, Levenberg–Marquardt
     damping (Madsen et al. §3.1 update rule) takes over for that iteration.
   * Converged when ‖r‖∞ ≤ 1e-12 (scaled), checked again after unscaling
     against the model tolerance. Limit: 100 iterations.
   * Hidden spline parameters are clamped to the curve domain after each step.
4. **Diagnose** at the final point with the scaled Jacobian:
   * Numerical rank ρ from the pivoted QR (|R_kk| ≤ 1e-10 · |R_11| is zero),
     confirmed by Jacobi SVD for components of up to 200 unknowns.
   * DOF = unknowns − ρ − (hidden unknowns); summed over components.
   * Free unknowns: unknown i is free if its row of an orthonormal null-space
     basis N of J has norm > 1e-8. An entity is *free* if any of its unknowns
     is; otherwise *constrained*. (Points of a fully constrained rectangle are
     constrained; a lone line's points are free.)
   * Redundancy: if equations > ρ, QR with column pivoting of Jᵀ picks a
     maximal independent set of equations in priority order (earlier
     constraints first, implicit equations first of all); the constraints
     owning the remaining equations are *redundant*.
   * Conflict: if the solve did not converge, compute the least-squares
     residual r*; the constraints whose equations carry |r*ᵢ| > 1e-9 form the
     *conflicting* set. Geometry is then left at the input (the seed), and the
     status is `Conflicting`.
5. **Write back** solved positions (in metres) into a copy of the sketch.

Redundant but consistent constraints (e.g. a horizontal line between two
points already constrained to equal y) solve, but the status is
`Redundant` and the UI shows them as over-defined. Users can delete them.

## 6. Branch stability

Minimum-norm steps from the current geometry pick the nearest solution. A
dimension change that fails to converge from the seed is retried by
continuation: the change is applied in 2, 4, … 64 equal steps, each solved
from the previous result. Orientation-sensitive constraints (tangent sides,
angle signs) store their branch in the constraint, so regeneration cannot
silently flip them.

## 7. Dragging

`SolveOptions::drag` lists points and target positions. The solver moves
those points to their targets in the starting state, then projects onto the
constraints with weight 1e6 on the dragged unknowns (they move last). A free
point ends exactly on the cursor; a constrained one ends at the nearest
feasible position. The UI calls `solve` per pointer move from the last
solved state and commits the final state as one edit.

## 8. Tests

* Unit tests for each equation's Jacobian against central finite
  differences at random configurations.
* Known DOF counts: lone point 2, line 4, circle 3, arc 5, triangle with
  coincident corners 6, constrained rectangle (H/V + coincident) 4, with
  two dimensions and a fixed corner 0.
* Constructive property tests (the fuzzing strategy): generate random
  geometry, derive constraints that the geometry satisfies (dimensions
  measured from it), perturb the positions, solve, and require every
  constraint satisfied to tolerance; the same at scales 1e-4 to 1e3.
* Redundancy: append a constraint implied by others; it is reported and the
  others are not. Conflict: append an inconsistent dimension; status is
  `Conflicting`, the set contains it, and the geometry is unchanged.
* Branch stability: change a triangle's side length by ±20% repeatedly; the
  triangle never flips orientation. Dimension changes return to the start
  exactly when undone.
* Dragging: a free point lands on the target; a point on a fixed circle
  lands on the circle nearest the target.
* Determinism: solving twice gives bit-identical output.
* Timing smoke tests for 50- and 300-entity sketches (reported, with
  generous limits in debug builds).

## 9. As built: differences from the plan above

The implementation follows §§3–8 with these changes, each made because a
test or the fuzzers showed the plan was not enough.

* **Substitution presolve** (`reduce.rs`). Before the numerical solve,
  equalities of two unknowns (coincident points, horizontal/vertical,
  concentric, equal radii) and unknowns fixed to constants (fixes, radius
  and diameter dimensions) are eliminated exactly, as SolveSpace does. Rank
  and null space transfer exactly to the full system. An elimination that
  contradicts an earlier one stays in the numerical system, where it is
  diagnosed as a conflict. Interpolated-spline control points are eliminated
  exactly by the interpolation equations.
* **DOF** is computed from the rows of the null-space basis that belong to
  visible unknowns, so hidden spline parameters never count.
* **Tolerances** account for lever length and f64 rounding noise. The final
  check (in metres, against the model tolerance) also allows the rounding
  bound for very short directions far from the origin.
* **Final acceptance check.** A result is accepted only if every constraint
  holds and nothing is degenerate at the result's own scale (threshold × 8
  margin): no direction a constraint needs, and no collapsed curve. That
  every line, circle and arc keeps its extent is an implicit requirement: a
  line whose end points meet, a circle or arc of zero radius, or an arc
  whose end has come onto its start (unless it is so in the input) is a
  collapse, and the input must pass the same check before the fallback below
  may return it, so a zero-length input line never solves. A collapse is
  attributed like a conflict: the constraints of the collapsed curve's
  connected part that hold it collapsed are tried in reverse priority order
  (latest first), each trial solving that part alone without the
  constraint, and the first whose removal gives a solution with nothing
  collapsed is reported as `Conflicting` (one per collapsed part), geometry
  left at the input — horizontal then vertical on one line blames the
  vertical. "Hold" is decided to first order at the collapsed state: the
  collapsed vector's (or radius's) gradient is written as the minimum-norm
  combination of the equation gradients (complete orthogonal
  decomposition), and the constraints with a nonzero coefficient hold it;
  equations that need a direction or length that has collapsed too add
  their constraint and the constraints holding that vector, found the same
  way. This keeps blaming to one or two trial solves in large sketches (at
  most 32 in any case). If no single removal helps, the solve fails as
  `Degenerate` naming the curve. A converged collapse is never answered by
  the fallback, and a solve that fails heading into a collapse (a radius
  cannot reach zero) is attributed the same way. Collapses that the eliminations already produce before
  solving (an arc start made coincident with its centre, a needed direction
  of zero length) are reported `Degenerate` by the start checks, as before.
  Otherwise a result failing the check is `NotConverged`; the solver never
  reports success on geometry a second solve would reject.
* **Input fallback.** If a solve without drags fails but the input already
  passes the final check, the input is returned as the solution, so
  re-solving a settled sketch always succeeds. Duplicates that agree within
  the model tolerance are reported as `Redundant`, not `Conflicting`.
* **Continuation** is adaptive and also interpolates fixed values. A
  turn-limited retry (steps that rotate a direction by more than π/4 are
  shortened) and foot-point restarts for point-on-curve constraints handle
  seeds far from the solution.
* **Rank** is first bounded rigorously from the pivoted QR's singular-value
  bounds; the Jacobi SVD runs only when the bounds are inconclusive. Both
  give the same rank. This took the 300-entity chain from 352 ms to 11.6 ms
  (release).
* **Rows that cancel are zero.** A Jacobian row whose entries sum (in
  absolute value) to at most 1e-10 of its equation's own gradient terms is
  set to zero before any rank decision. Such a row is rounding noise: the
  members of an equality class pulling equally both ways, or a derivative
  that is zero analytically (an angle's with respect to a vertical line's
  length). Alone in its component, the relative rank test above would
  count it as full rank against itself, so an angle implied by horizontal
  and vertical constraints went unreported when the input was placed
  exactly and was reported when it was nudged. The mirror fuzz (design
  note 0007 §9) found this; `constraints_implied_by_eliminations_alone_are_redundant`
  guards it. Rows of legitimate equations are not that small: a row is
  smallest against its terms when a free point is measured against a short
  fixed line (its terms grow as 1/length), and every line shorter than
  8 × 1e-10 of the sketch has collapsed (above), which keeps the ratio
  above 1e-10 for geometry within the sketch's extent.
* **Spline parameters** use an active set at the domain ends; spline
  corners (repeated knots) are supported; the start point for a
  point-on-spline is found by sampling every knot span and refining each
  local minimum (golden section, then Newton).
* **Limits.** Coordinates beyond 1e6 m are rejected; iterates leaving a
  reach box of 1e3 sketch sizes fail; a radius may at most halve per step.
  An arc with zero start radius is rejected up front as `Degenerate`, and a
  zero-length line in `resolve_curves` is `Degenerate`.
* **Tangent equations** are normalized (divided by the lengths involved),
  like parallel and perpendicular.
* **Reporting.** A constraint that is both conflicting and redundant is
  reported only as conflicting. `NotConverged` also covers results that
  diverged, left the reach box or were degenerate.

**Known limitation.** Conflict detection is local: from a seed far from any
solution, the solver can report `Conflicting` although another branch
would satisfy every constraint. It never reports success on a violated
constraint.

**Measured** (final code): about 400k constructive property cases at scales
1e-4 to 1e3 and 500k garbage-input fuzz cases with zero failures (no panic,
no NaN, every success verified by an independent checker, re-solving every
output succeeds). With random perturbations of 2–5% of the sketch size every
case solves; at 10% about 1 in 40,000 fails; at 30% about 1 in 600 fails,
mostly by flipping to another branch. A failed solve always says so.
