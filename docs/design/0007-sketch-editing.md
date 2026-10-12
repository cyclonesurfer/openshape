# Design note 0007 — Sketch editing tools

Modules: `crates/ui/src/sketch_ops/` (pure functions from a sketch to a
`SketchEdit` batch), `crates/ui/src/sketch_mode.rs` (the interaction),
`crates/kernel/src/arrangement` (`curve_contacts`), and the derived
constraints of `crates/sketch-solver` (§9).

## 1. Scope

Tools that create or change sketch geometry beyond the basic curves of
Phase 2:

* shapes: centre-point rectangle, aligned (three-point) rectangle, slot,
  regular polygon, tangent arc;
* inference while drawing lines: horizontal, vertical, parallel,
  perpendicular, tangent;
* trim, extend, split;
* sketch fillet;
* mirror, offset, linear and circular patterns;
* dimension labels the user can drag.

The behaviour (what a user sees and can do) follows what parametric CAD
systems publicly document for their sketchers; the geometry and the
constraint formulations below are our own, from the literature in §2.

## 2. Literature

* Variational (constraint-driven) sketch geometry: Light & Gossard,
  "Modification of geometric models through variational geometry", *CAD*
  14(4), 1982; the solver itself is design note 0005.
* Curve intersection: design note 0006 §4 (closed forms for lines and
  circles; Sederberg & Nishita, "Curve intersection using Bézier clipping",
  *CAD* 22(9), 1990, for splines).
* Fillets: a circle of radius r tangent to two curves has its centre on
  both curves' offsets at distance r — Faux & Pratt, *Computational
  Geometry for Design and Manufacture*, 1979, ch. 1–2; the planar case of
  the rolling-ball blend of Rossignac & Requicha, "Constant-radius blending
  in solid modeling", *Computers in Mechanical Engineering* 3(1), 1984.
* Offsets of line–arc profiles, their trimming at corners, and the
  self-intersections that large distances create: Tiller & Hanson,
  "Offsets of two-dimensional profiles", *IEEE CG&A* 4(9), 1984; Held, *On
  the Computational Geometry of Pocket Machining*, LNCS 500, 1991, ch. 2–3.
* Offsets of polynomial curves are in general neither polynomial nor
  rational: Farouki & Neff, "Analytic properties of plane offset curves",
  *CAGD* 7, 1990; approximation methods compared in Elber, Lee & Kim,
  "Comparing offset curve approximation methods", *IEEE CG&A* 17(3), 1997.
  Hoschek & Lasser, *Fundamentals of CAGD*, 1993, §12 (offsets).
* B-splines: Piegl & Tiller, *The NURBS Book*, 2nd ed.

## 3. Principles

1. **One batch, one undo step.** Every tool action is one `SketchEdit`
   batch (ADR 0011): it adds, moves and deletes entities and constraints
   atomically.
2. **Exact geometry first, constraints second.** A tool computes the
   resulting geometry exactly (closed forms for lines and circles), writes
   it into the batch, and adds constraints that this geometry already
   satisfies. The solve after the edit then moves nothing (up to rounding).
3. **Checked.** After applying a batch, the tool re-solves. If the sketch
   solved before and does not solve after, the batch is undone and the tool
   reports why. A tool never leaves the sketch in a state it did not intend.
4. **No silent loss of intent.** When an edit has to drop constraints that
   no longer make sense (the length of a line that was trimmed), the status
   bar says how many were removed and why.
5. **Degrees of freedom are what the user expects.** Each tool's
   constraints are counted below; the tests assert the counts and that no
   constraint is redundant.
6. **Refuse what cannot be exact.** Splines cannot be trimmed, extended,
   filleted or offset exactly with the sketch's spline types (§7.5, §10.5);
   those tools refuse splines with a message instead of approximating.

## 4. Shapes

All shapes are built from the existing entities and constraints.

* **Centre-point rectangle** (click centre, click a corner): four lines
  joined by coincident corners, two horizontal and two vertical; two
  construction diagonals whose ends are coincident with the corners; the
  clicked centre is a point constrained as the first diagonal's midpoint
  (it then lies on the second diagonal too, so no constraint says so).
  DOF: sides 16 + diagonals 8 + centre 2 − corners 8 − diagonal ends 8 −
  horizontal/vertical 4 − midpoint 2 = 4 (centre, width, height).
* **Aligned rectangle** (click two corners along one side, then the
  width): four lines, coincident corners, `Perpendicular(s0, s1)`,
  `Parallel(s0, s2)`, `Parallel(s1, s3)`. DOF 5 (position, angle, length,
  width).
* **Slot** (click the two arc centres, then a point on the side): a
  construction line between the centres; two arcs centred on its ends; two
  lines joining the arc ends; coincident joints; tangency at the four
  joints; equal radii. DOF 5 (two centres and the radius).
* **Regular polygon** (click centre, click a vertex; the number of sides is
  a tool setting, 3 to 64, default 6): a construction circle about the
  centre through the vertex; n lines joined end to end, each start point on
  the circle; n − 1 equal-length constraints. Equal chords of a circle
  subtend equal angles, so the polygon is regular. DOF 4 (centre, radius,
  rotation).
* **Tangent arc** (click the end of a line or arc, then the arc's end):
  with start P, unit tangent d of the curve leaving P and end E, the centre
  is P + k n, n the left normal of d and k = |E − P|² / (2 n·(E − P)); the
  arc turns left (counter-clockwise from P) when k > 0. Constraints:
  coincident with the curve's end and `Tangent` at that point. Refused
  when E is (nearly) on the tangent line, where the circle degenerates.

## 5. Inference while drawing lines

While the second point of a line moves freely (no snap), its direction is
compared with (in order):

1. the sketch axes (horizontal, vertical);
2. when the line starts at the end of an arc: the arc's tangent there
   (tangent);
3. when the line starts at the end of a line: that line (perpendicular,
   then parallel);
4. the last line the pointer passed over during this tool (perpendicular,
   then parallel).

The first direction within 3° wins: the point is moved onto it (keeping its
distance along the direction) and the matching constraint is added with the
line. The canvas shows the inferred relation's glyph next to the pointer.

## 6. Curve contacts (kernel)

`Kernel::curve_contacts(a, b)` returns where two profile curves meet:
transversal crossings, tangencies, end points of one lying on the other and
the two ends of an overlap, as parameter pairs in the arrangement's
conventions (design note 0006), merged when closer than the tolerance and
sorted by the parameter on `a`. It reuses the arrangement's pair
intersection (closed forms for line and circle pairs, Bézier subdivision
with Newton refinement for splines), so tools and region detection agree on
where curves meet. Every contact satisfies |A(tₐ) − B(t_b)| ≤ tolerance.

## 7. Trim, extend, split

### 7.1 Pieces

The clicked curve C (line, arc or circle) is cut by every other curve of
the sketch, construction geometry included. Its contact parameters, with
contacts within the tolerance of C's own ends dropped, divide C into pieces;
the piece under the pointer (the closest point's parameter) is the one
acted on. For a circle the parameters are cyclic.

### 7.2 Trim

* C has no contacts: C is deleted (with points only it used).
* The piece is at the start of an open curve: C's start point moves to the
  first cut; constraints on that point are removed (that end no longer
  exists); the point is attached to the cutting curve (coincident with a
  point there, else on the curve).
* At the end: symmetric.
* In the middle: C keeps the part before the piece; a new curve D takes the
  part after it, with a new start point (attached to its cutter) and a new
  end point that takes over every constraint of C's old end point. D is
  tied to C: a line by both its ends lying on C (collinear), an arc by
  `Concentric` and `Equal`. D copies C's construction flag.
* A circle needs at least two cuts; the rest of the circle becomes an arc
  on the same centre point; every constraint on the circle moves to the arc
  (all circle constraints mean the same for an arc: they use the full
  circle).

Constraints on C itself after a trim:

| constraint on C | kept? |
|---|---|
| horizontal, vertical, parallel, perpendicular, angle | kept (D is collinear) |
| point on curve, tangent to a circle, normal, symmetric about C, point–line and line–line distance | kept (they use the infinite line or the full circle) |
| radius, diameter, equal radii, concentric | kept (the circle is unchanged) |
| tangent at a point | stays with the curve that keeps the point |
| length, equal length, midpoint | removed (the segment changed) |

### 7.3 Extend

Clicking near an end of a line or arc extends that end to the nearest
contact beyond it: along the line's ray (a proxy segment as long as four
times the sketch's size), or along the arc's circle (the rest of the
circle). The end point moves there and is attached to that curve;
constraints on the old end point, and length or equal-length constraints
of a line, are removed and reported. No contact: refused with a message.

### 7.4 Split

Clicking a line or arc splits it at the pointer's closest point into C and
D as in the middle-piece trim, with D's start coincident with C's new end.
Circles cannot be split at one point.

### 7.4a As built (trim, extend, split)

Implementing and fuzzing §7.1–7.4 changed four details; the reasons are
the solver's rank test, which reports a constraint as redundant when its
equation is implied to first order.

* **Attaching a cut end** uses a point-on-curve constraint, never a
  coincidence with a point of the cutter. At a T junction the cutter's
  end usually already lies on the cut curve; a coincidence would say that
  twice.
* **Touching is not crossing.** A cut end is attached only to a curve that
  crosses the cut curve there at more than about 0.6° (sine 0.01), and is
  not constrained tangent to it. Where curves touch (a tangency, a curve
  running along the same line or circle), point-on-curve has a double
  root: redundant when the tangency is constrained, even indirectly, and
  ill-conditioned when it is not. Such a cut still cuts; its end is left
  free along its own curve.
* **Arcs share the centre point.** The new arc of a middle trim or a split
  uses the cut arc's centre point entity rather than a concentric
  constraint: a middle trim adds equal radii, a split adds nothing (its
  joint is a coincidence of two points on the one circle).
* **Joints survive a moved end.** When an end that other curves met at
  moves, its coincidences are removed but the remaining points stay
  coincident with one another (links are added only where the other
  coincidences do not already join them), and a tangency of two other
  curves at that joint moves to a remaining point. A tangency of the cut
  curve at that end goes with the end.

Every plan is committed through a checked edit: if the sketch solved
before and would fail or gain a redundant constraint after, the edit is
refused with the reason. A property test applies random trims, extends
and splits to random lines, arcs and circles and to constrained shapes
(rectangles, slots, polygons, joined chains): every plan commits, the
sketch stays solved without redundancy, no point the plan did not move
moves, and no curve shrinks below 1 µm. Region detection can still refuse
a sliver thinner than the tolerance that an edit creates (for example a
line extended past the tangent joint of a large arc next to another edge);
that is the arrangement's documented refusal (design note 0006), not a
wrong region.

### 7.5 Splines

Neither spline type of the sketch (uniform clamped control-point splines,
and interpolated splines) can represent an arbitrary piece of a spline
exactly: knot insertion (P&T A5.1) produces a non-uniform knot vector.
Trim and extend refuse splines (an uncut spline is deleted by trim, as any
uncut curve); a spline type with explicit knots is future work.

## 8. Fillet

Click a corner point where exactly two lines or arcs end (their end points
coincident, or one shared point) and the fillet is made with the tool's
radius; the new radius dimension then opens for editing.

**Geometry.** For two lines with unit directions u₁, u₂ leaving the corner
V at angle θ: tangent points V + uᵢ r / tan(θ/2), centre on the bisector at
r / sin(θ/2). For pairs with an arc, the centre is an intersection of the
two curves' offsets at r (line: parallel at ±r; circle: radius R ± r);
among the candidates, the one whose tangent points lie on both pieces next
to V and nearest V is used. Tangent curves at V (θ = π) and radii too large
for either piece are refused.

**Edits.** The two end points move to the tangent points (an arc end stays
on its circle). A new arc joins them, with coincident ends, `Tangent` at
both points and a `Radius` dimension. If the corner carried constraints
(other than the two ends being coincident), a *virtual sharp* point is
added at V, constrained onto both curves, and those constraints move to it;
a `Length` of a filleted line becomes a `Distance` from the virtual sharp
to the line's far end. Equal-length and midpoint constraints on the
filleted lines are removed and reported. DOF are unchanged by a fillet
(the arc's five freedoms meet five equations; the virtual sharp's two meet
two).

### 8.1 As built

* The candidate centres are the intersections of the two curves' offsets
  at ±r (parallel lines; concentric circles R ± r). A candidate is kept
  when both tangent points lie strictly inside their curves (1 µm from
  either end) and the arc continues the first curve smoothly into the
  second (travelling along the first towards the corner, the arc must
  arrive on the second heading away from it); the nearest to the corner
  wins. One rule covers line–line, line–arc and arc–arc.
* The arc's ends are new points coincident with the moved curve ends
  (not the same point entities), so a later trim can move either end on
  its own.
* Construction curves ending at a corner of two real curves (the
  diagonals of a centre-point rectangle) do not count as sides; they stay
  coincident with the virtual sharp.
* A smooth joint (tangents within 0.6°) or a fold-back is refused, as is
  a radius whose tangent points would not fit on the curves.
* A constraint that would collapse onto the sharp twice (between the two
  ends) is dropped, and one putting an end on the other curve is not
  repeated (the sharp already lies on both).

A property test fillets random corners of rectangles, polygons, joined
chains and line/arc corners with random radii: every fillet commits
through the checked edit, the sketch stays solved without redundancy, the
degrees of freedom change only by what removed equal-length (one) and
midpoint (two) constraints held, and no point but the two moved ends
moves.

## 9. Derived geometry in the solver

Mirror, offset and patterns create copies that must follow their source.
Expressing that with the ordinary constraints creates redundancy: an arc's
mirror image satisfies |end − centre| = |start − centre| automatically, so
the arc's implicit equation and two symmetric constraints say the same
thing, and the solver would rightly flag a redundant constraint. Tangent
chains offset by one distance close up only because the original closes,
which again shows as a dependent equation.

So these tools use **derived constraints**: each one gives, for every image
point, exactly two equations that fix it from the source geometry, with a
Jacobian block that is regular whatever the rest of the sketch is:

* `Mirror { axis, pairs }`: entity pairs (source, image) of the same kind;
  the images' points are reflections of the sources' points (the
  symmetric-point equations of design note 0005 §4); an image arc's start
  and end are the reflections of the source's end and start (a reflection
  reverses orientation); an image circle's radius equals the source's.
* `LinearPattern { direction, copies }` and `CircularPattern { centre,
  angle, copies }`: copies (source, image, k); image points are the source
  points translated by k times the direction line's vector, or rotated
  about the centre by k times the angle — both linear in the unknowns.
* `Offset { value, pieces }`, §10.

The image arcs' implicit radius equations are dropped (they are implied),
and image circles' radii are eliminated by substitution. Images therefore
add no degrees of freedom and no redundancy; constraints the user adds on
images are judged like any others (a dimension on a mirror image of a
dimensioned curve is redundant, as it should be).

## 10. Offset

### 10.1 Selection

Clicking a curve selects its *chain*: the maximal sequence of lines, arcs
(and a lone circle) joined end to end through points where exactly two of
them meet. Moving the pointer sets the side and the distance (the pointer's
distance to the clicked curve); a click creates the offset with a distance
dimension the user can edit.

### 10.2 Pieces and sides

The chain is oriented; each piece is offset to the chain's left or right,
expressed in the piece's own direction as `left` (a line: the left of
start → end; an arc or circle, counter-clockwise: towards the centre,
radius R − d; right: R + d).

### 10.3 Joins (Tiller & Hanson)

At each joint V between consecutive pieces:

* **Tangent joint** (the originals share a tangent at V): both offset ends
  are V + d·n, n the common normal on the offset side.
* **Corner**: the offset curves, extended (lines infinitely, arcs as full
  circles), are intersected; the intersection nearest V is the join. This
  trims overlapping offsets at concave corners and extends them at convex
  ones. If the extended offsets do not meet (possible only with an arc
  offset towards its centre), the offset is refused.
* **Free end** of an open chain: the offset end lies on the original's
  normal through the original end.

### 10.4 Equations

Each offset piece's end point is its own point entity, fixed by two
equations: one "on my offset curve", and one placing it along the curve —
on the neighbour's offset curve (corner), or on the normal line at V (free
end, tangent joint). An offset arc's centre equals the original's centre;
its start gets the explicit "on my offset circle" equation and its end the
arc's own implicit equation, so no equation is implied by another. An
offset circle's radius is the original's ± d. Corners meet transversally,
so their two equations are independent; tangent joints never use two
"on offset curve" equations at one point, so no double root appears.

### 10.5 Validity

At creation, and after every later solve, the offset must be valid: every
line piece keeps the direction of its original (a reversed piece means the
distance swallowed it), every arc keeps a positive radius, and no two
pieces cross except consecutive pieces at their join (checked with
`curve_contacts`; Held ch. 3). An invalid offset is refused at creation
and, after a later edit, reported as a sketch error naming the offset;
nothing is silently trimmed away. Splines are refused: their offsets are
not splines (Farouki & Neff), and an approximating fit with a guaranteed
tolerance (Elber, Lee & Kim) would need a constraint that holds only
approximately; that is future work.

## 11. Patterns

* **Linear pattern**: select curves, then a direction line (made by the
  tool as a construction line the user can dimension) and a count; copy k
  is the source moved by k times the direction line's vector.
* **Circular pattern**: select curves, a centre point and a count; the
  angle between copies is a dimension (default 360°/count); copy k is the
  source rotated by k times the angle.

Copies keep the construction flag of their source and can be selected and
constrained like any curve; deleting a copy removes it from the pattern.

## 12. Dimension labels

A dimension's label position is stored per dimension as an offset in plane
coordinates from the dimension's anchor (the midpoint of what it measures),
in `SketchDef::labels` (`SketchEdit::PlaceLabel`). The dimension line
follows the label (for a length, the line runs parallel to the measured
segment through the label). Dimensions without a stored position are
placed automatically as before (outside, away from the sketch's centroid).

## 13. Tests

* Shapes: entity and constraint counts, DOF exactly as in §4, no redundant
  constraint, regions and extruded volumes exact.
* Inference: each relation is added (and only then) when the pointer comes
  within 3° of it.
* Contacts: line, arc, circle and spline pairs against closed forms;
  tangencies and overlaps; parameter conventions round-trip through the
  arrangement's evaluation.
* Trim/extend/split: each case of §7 (start, middle, end, uncut, circle)
  with exact resulting geometry, constraint transfer per the table, DOF,
  and that a trimmed sketch's regions extrude to the expected volume.
* Fillet: line–line, line–arc and arc–arc corners against closed forms; the
  fillet radius dimension can be edited afterwards; virtual sharps keep
  dimensions; refused cases.
* Derived constraints: images follow their sources under drags and
  dimension changes; DOF unchanged by adding images; no redundancy; image
  arcs keep their orientation.
* Offset: rectangle, slot and L-shaped chains, inward and outward; closed
  form areas; refused cases (an arc swallowed, a crossing offset, splines);
  an edit of the original that makes an offset invalid is reported.
* Property tests: random convex polygons offset outward by d have area
  A + P·d + d²·Σ tan(θᵢ/2) (θᵢ the exterior angles); random polygons
  mirrored about random lines keep area and reverse orientation.
