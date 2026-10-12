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
* move, rotate and scale (transform), with copies;
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
* **Control-point spline** (click control points; Enter or double-click
  finishes): a `ControlSpline` on the clicked points, clamped with uniform
  interior knots (P&T eq. 2.13, clamped form), degree min(3, n − 1); at
  least three points, so the tool never makes a straight line. The clamped
  form starts at the first control point, ends at the last and is tangent
  to the polygon there (P&T §3.2). Clicking the first point again (it
  snaps within the point pick radius) finishes the curve on that same
  point entity, so the curve closes: with clamped knots it meets itself
  with a corner, not smoothly. A smooth closure needs a periodic spline
  (P&T §12.2: wrapped control points on uniform unclamped knots), which
  the sketch cannot express yet. Mirror and patterns copy it exactly (a
  B-spline is affine invariant, P&T §3.2: the image of the curve is the
  curve on the images of its control points). DOF 2n, no constraints.
  While drawing, the control polygon is shown dashed with the curve to the
  pointer; afterwards the polygon stays visible, dimmed, so its points can
  be dragged.

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

### 9.1 As built (mirror)

* **Equations.** For each image point q of a source point p: the midpoint
  (p + q)/2 lies on the line (its signed distance is zero) and q − p is
  normal to it (its projection on the line's direction is zero). In q the
  two rows are the line's normal and direction, so the block is regular
  whatever else holds the sketch. An image arc's start and end are the
  images of the source's end and start; image arcs' implicit equations
  are not emitted, and image circles' radii are substituted by the
  source's.
* **Validity** (refused with a message naming the constraint): the line
  cannot be a source or an image; an entity cannot be its own image, the
  image of two entities, or an image of two mirror constraints; a pair must
  be the same kind (splines the same form and size); a point cannot be the
  image of two points, nor both a source and an image within one
  constraint; an end of the line cannot be an image; no entity may be
  derived, through any chain of mirrors, from its own image.
* **Shared points.** A source point on the line is shared with the image,
  so a mirrored half profile closes. Such a pair (p, p) gives one
  equation, p on the line, and the edit that adds the mirror deletes the
  point-on-line constraints that held p there (and those on points
  coincident with it), so nothing is said twice. A shared point
  coincident with an end of the line gets no equation (the line holds
  it), and of shared points coincident with each other only one is held
  (the fuzz found this: a slot's centre line and its arc meet at
  coincident points, and mirroring both about a line through that joint
  held the joint twice). A point held on the line by a midpoint constraint is not shared:
  that constraint already says "on the line" and more, so its image gets
  a point of its own, which the mirror keeps on the same spot. The editor
  commits with sharing and, if the checked edit refuses it (a constraint
  would become redundant), again without; in 3,000 fuzz cases (4,650
  mirrors) the fallback was not needed once, but it keeps the tool from
  failing on a case nobody foresaw.
* **Coincident sources.** Two coincident source points (a joint drawn as
  two line ends) get two image points; the mirror makes them coincident
  as a consequence, so no constraint is added and none is redundant.
* **One constraint per line.** Mirroring again about the same line
  replaces its constraint by one holding the old and new pairs, and image
  points are reused (both ways: a reflection is its own inverse), so
  curves mirrored one at a time join exactly as their sources do. Sources
  already in a pair about the line, and entities that are their own image
  (a circle centred on the line, a line along it or symmetric across it, a
  point on it), are skipped and counted in the message.
* **Deleting.** Deleting a source or an image removes only that pair
  (`Constraint::without`); deleting the line removes the constraint, and
  the images stay as free geometry.
* **Tool.** With a selection, the first line clicked mirrors it in one
  edit; without one, the first click chooses the line and each later
  click mirrors the curve or lone point under the pointer. The tool stays
  on the line; Esc chooses another, a second Esc leaves the tool. Hovering
  previews the images.

Two property tests guard it. `crates/sketch-solver/tests/mirror.rs`
mirrors random curves of random constrained sketches about a fixed or
free line: the mirror adds no freedom and no redundancy and every
equation holds after the solve. `crates/ui/tests/sketch_mirror.rs`
mirrors random rectangles, polygons, slots, circles, arcs and points in
random batches about random lines, often through one of their points,
through the editor's planner: every mirror commits through the checked
edit, the sketch stays solved without redundancy, the degrees of freedom
drop only by one per newly shared point, no source moves, and every image
point is its source's reflection to 1e-9. Both ran 3,000 cases clean. The fuzz also found a solver defect,
fixed separately (design note 0005 §9, "Rows that cancel are zero").

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

### 10.6 As built (offset)

* **Constraint.** `Offset { distance, pieces, smooth }`: pieces are
  `(source, image, left)`, two lines, two arcs or two circles, with the
  side in the source's own direction (§10.2); `smooth` lists the image
  points at joints whose sources meet tangentially. The distance is the
  constraint's value, so it is edited like any dimension and labelled from
  the middle of the first source across to its image.
* **Equations** (a change from §10.4). Image arcs and circles share their
  source's centre point. Each image end point gets exactly two equations
  from the joint it sits at: at a corner, on both pieces' offset curves;
  at a free end, on its piece's offset curve and on the normal at its
  source's end; at a smooth joint, on the first piece's offset curve and
  on the normal at the joint (the second piece's offset passes through
  that point as long as the sources stay tangent). An image arc's own
  equation is dropped instead of being used at one end: with both ends
  placed on the same offset circle it is implied. An image circle's radius
  is the source's less or more the distance. Images therefore add no
  freedom. A joint the planner wrongly took for a corner while it is
  tangent would put two on-curve equations at a double root, which the
  rank test reports as redundant; a test that first drew such a joint by
  mistake (a line leaving a quarter arc vertically) showed exactly that.
* **Validity, in the solver** (`crates/sketch-solver/src/offset.rs`),
  after every solve, as part of accepting it: every image end on its
  offset curve (fails at a smooth joint whose sources are no longer
  tangent), no line image reversed, no arc image turned back (its sweep
  measured from the source's, each end moved by less than half a turn),
  no arc end held on a normal lying across the centre (the normal line
  meets the offset circle twice), no image circle shrunk to nothing, and
  no two images crossing, touching or overlapping except consecutive
  ones at their shared joint. A failure is `SolveError::Degenerate` with
  a message naming the offset and the curve. The solver computes these
  contacts itself (it does not depend on the kernel); the editor fuzz
  checks every committed offset against the kernel's `curve_contacts`.
* **Images move first** (design note 0005 §3). Without this, changing the
  distance spread the change over the sources as well (a rectangle's
  corner moved by a third of the change in the first test run): the
  ordinary step is the smallest change of all points together. Now the
  sources stay exactly where they were and only the offset moves; a
  dragged image still drags its sources.
* **Planner** (`crates/ui/src/sketch_ops/offset.rs`). The chain runs
  through coincidence classes where exactly two curves of the same kind
  (construction or not) end; a spline ends a chain and cannot be offset;
  a circle is a chain of its own. The side and distance come from the
  pointer's distance to the clicked curve, rounded to a 1, 2 or 5 step of
  about 1.5 screen points. Joins follow §10.3: smooth when the tangents
  agree to a sine of 1e-7, otherwise the meeting of the offset carriers
  nearest the joint; the planner refuses when they do not meet or the
  distance exceeds an arc's or circle's radius. The edit is committed
  through the checked edit (refused if the sketch would stop solving or
  gain a redundant constraint, which the solver's validity check turns
  into a refusal for a crossing or turned-back offset), then the distance
  opens for typing. Hovering shows the chain, then the offset, already
  checked by a trial solve.
* **Deleting** a source or an image removes its piece; joints it shared
  become free ends, whose images move to the normals at their sources'
  ends. Deleting every piece removes the constraint.

Tests: `crates/sketch-solver/tests/offset.rs` (rectangles in and out,
the distance as a dimension, redundancy of a dimension on an image,
chains with smooth joints, corners and free ends, circles, refusals,
deletion, and a property test over random convex polygons) and
`crates/ui/tests/sketch_offset.rs` (the tool end to end, regions and an
extruded volume exact, and a fuzz: random rectangles, polygons, slots,
circles and line chains offset by random distances to random sides; every
offset either commits with the sketch solved, nothing redundant, no
freedom added, no source moved, every image at the distance and no two
crossing by the kernel's intersections, or is refused naming the
offset). At 3,000 cases each: the property test solved every polygon
exactly; the fuzz committed 2,738 offsets, and the solver's check
refused 204 (a piece turned back or images crossing) and the planner 58
(the distance larger than an arc's or circle's radius), with no
failure.

## 11. Patterns

* **Linear pattern**: select curves, then a direction line (made by the
  tool as a construction line the user can dimension) and a count; copy k
  is the source moved by k times the direction line's vector.
* **Circular pattern**: select curves, a centre point and a count; the
  angle between copies is a dimension (default 360°/count); copy k is the
  source rotated by k times the angle.

Copies keep the construction flag of their source and can be selected and
constrained like any curve; deleting a copy removes it from the pattern.

### 11.1 Constraints and equations

Two derived constraints (§9), with copies `(source, image, k)`, k ≥ 1:

* `LinearPattern { direction, copies }`: for each source point p and its
  image q in copy k, q − p − k·(e − s) = 0, s and e the direction line's
  start and end: two equations linear in the unknowns, the identity in q.
* `CircularPattern { centre, angle, copies }`: q − c − R(k·θ)(p − c) = 0,
  c the centre point and θ the angle, a dimension (its value): linear in
  the unknowns for the given angle, a rotation in q.

Both maps keep orientation, so an image arc's start and end are the
images of its source's start and end; image arcs' implicit equations are
dropped and image circles' radii are their sources' (substituted). A
source point that is the centre (or coincident with it) is its own image
in every copy of a circular pattern: the copies share it and it gets no
equation, so spokes drawn from a hub stay joined at it. Each use of a
pattern tool makes one constraint; several may share a centre point.

Validity, refused naming the constraint: the direction line and the
centre cannot be sources or images; a pair must be the same kind; an
entity cannot be its own image or the image of two entities; k·θ must
not be a whole number of turns (a copy would lie on its source), nor may
two copies of one source have the same k. As with mirrors, images move
first (design note 0005 §3): changing the angle moves only the copies.

### 11.2 Tools

The count (instances, the source included) is a field of the sketch
dialog while a pattern tool is active, like the fillet radius. With
curves selected:

* **Linear**: two clicks, snapping like a line, give the first step: the
  tool adds a construction line between them (the direction line) and
  count − 1 copies, then opens the line's length for typing, which adds
  a length dimension.
* **Circular**: one click, snapping to points, gives the centre (a new
  point if it snaps to nothing); copies are spread at 360°/count and the
  angle opens for typing.

Hovering previews the copies; edits go through the checked edit, so a
pattern that would make the sketch unsolvable or redundant is refused.

### 11.3 As built (patterns)

* **Constraints.** `LinearPattern { direction, copies }` and
  `CircularPattern { centre, angle, copies }`, copies `(source, image, k)`.
  Each copy point gets two equations linear in the unknowns with
  coefficients summing to zero (so they hold under any shift of the
  sketch): q − p − k(e − s) = 0, and q − c − R(kθ)(p − c) = 0 with the
  rotation's cosine and sine computed once from the angle. Image arcs'
  implicit equations are dropped and image circles' radii are their
  sources' (substituted), as for mirrors (§9). The angle is the circular
  pattern's value, so it is edited and labelled like any dimension (an arc
  about the centre from a point of the first source to its first copy).
* **Validity** is checked when the solver indexes the sketch and refused
  naming the constraint: the line or centre as a source or image, a line
  end as an image, mismatched kinds, an entity its own image (only a
  point in the centre's coincidence class may be), an image in two copies,
  a point both a copy and a source, k = 0, and an angle that is not
  finite, zero, a turn or more, or turns some copy k onto its source
  (k·θ within 1e-9 rad of a whole number of turns). Copies of one k are
  checked like a mirror's pairs; across k, every image point must be new.
* **Images weigh nothing** in a solver step (design note 0005 §5): copy
  points are image slots, so a new angle moves only the copies and a new
  step moves the direction line's end and the copies, not the sources.
  Copies that nothing else holds are *peeled*: computed from the solved
  sources after the solve rather than solved for (design note 0005 §5
  step 2), so large patterns cost little.
* **Found by the tests.** The first editor test typed a linear pattern's
  step and the square being copied changed width: the old images-first
  rule only applied when moving the images alone solved the step, which a
  direction line's length change does not. The rule is now the exact
  limit (above), with a solver regression test whose source runs along
  the step. The first UI fuzz took 28 minutes: a pattern of n copies was
  one dense component of n times the source's unknowns, hence peeling.
* **Planner** (`crates/ui/src/sketch_ops/pattern.rs`). Sources are the
  selection as for mirrors (curves, and points not on a selected curve).
  Copy k of each source point is a new point, shared by every source of
  that copy that uses the point; a circular pattern's source points in
  the centre's coincidence class are shared by all copies instead. A
  linear pattern's line is a construction line between the two snapped
  clicks; a circular pattern's centre is the snapped point itself, or a
  new point constrained to what it snapped to. Sources a circular pattern
  would turn onto themselves (a point on the centre, a circle about it)
  are skipped, and a pattern with nothing left is refused. The count is
  2 to 100 instances (the Instances field of the sketch dialog, default
  3); a circular pattern's angle starts at 360°/count.
* **Tools.** Both are chosen with the selection made, and keep it as the
  sources. Hovering previews the copies after a trial solve; the click
  commits through the checked edit, returns to the Select tool and opens
  the line's length (a new dimension) or the angle (the pattern's own
  value) for typing. A typed angle is refused with the reason when it is
  zero, a turn or more, or turns a copy onto its source.
* **Deleting** a source or a copy removes that copy; deleting the line or
  the centre removes the constraint and leaves the copies as free curves.

Tests: `crates/sketch-solver/tests/pattern.rs` (copies of every kind
exact for both patterns, the angle as a dimension moving only the copies,
the line driving them, spokes sharing the hub, a dimension on a copy
redundant, every malformed pattern refused, deletion, and a property test
that random patterns add no freedom and no redundancy) and
`crates/ui/tests/sketch_pattern.rs` (the tools end to end: a square
copied along a typed step and extruded, spokes about their hub and a new
angle typed, circles about the origin following their source's radius,
deletion, refusals, and a fuzz of random patterns over random shapes).

## 12. Dimension labels

A dimension's label position is stored per dimension as an offset in plane
coordinates from the dimension's anchor (the midpoint of what it measures),
in `SketchDef::labels` (`SketchEdit::PlaceLabel`). The dimension line
follows the label (for a length, the line runs parallel to the measured
segment through the label). Dimensions without a stored position are
placed automatically as before (outside, away from the sketch's centroid).

### 12.1 As built (labels)

* Storage. `SketchDef::labels` holds `(dimension, offset)` pairs sorted by
  id (pairs, not a map, so a definition serializes to formats whose map
  keys must be strings; `#[serde(default)]`, so documents saved before it
  read with none). `SketchEdit::PlaceLabel { dimension, offset }` sets an
  offset, or with `None` returns the label to automatic placement; it is
  refused for a missing id, a constraint that is not a dimension, or a
  non-finite offset. `Delete` drops the labels of every dimension it
  removes, including those removed because what they measured went.
  Concurrent placements of one label converge on the op later in replay
  order, like concurrent value edits.
* Anchors: a length's or aligned/horizontal/vertical distance's two
  points' midpoint; for a point–line or line–line distance, the middle of
  the point (or first line's midpoint) and its foot on the other line; a
  radius's or diameter's centre; an angle's vertex (the lines'
  intersection); an offset's middle of its first source–image segment; a
  circular pattern's centre. The offset is in plane coordinates, so a
  label keeps its direction from the anchor when what it measures turns.
* Drawing with a placed label at L (the label's centre is always L):
  * length and aligned distance: the dimension line is parallel to the
    measured segment through L, reaching L when it lies past an end, with
    extension lines from the measured points to a quarter of the automatic
    offset beyond it;
  * horizontal or vertical distance: the dimension line at L's y (or x),
    likewise extended, with extension lines from both points;
  * point–line distance and offset: the measured segment plus a leader
    from its nearest point to L;
  * radius: from the centre towards L, out to L when L is outside the
    circle (to the circle otherwise); a diameter likewise, from the far
    side;
  * angle: an arc of radius |L − vertex| in whichever of the two equal
    vertical angles L is on, widened to L's direction when L is outside
    it, with extension lines along the rays where the lines stop short of
    the arc (also drawn for automatic labels now);
  * circular pattern: the arc from the first source to its first copy at
    radius |L − centre|, widened likewise, with both radii.
* Editing. In the Select tool, dragging a dimension's label moves it: the
  drawing follows the pointer and the release commits one `PlaceLabel`
  (one undo step); a drag that ends where it started commits nothing.
  Constraint glyphs are not draggable (a drag on one box-selects as
  before). The command "Reset dimension label positions" returns the
  selected dimensions' labels, or all placed ones when no dimension is
  selected, to automatic placement in one step.

Tests: `crates/core-model/tests/sketch.rs` (placing, replacing, clearing,
refusals, deletion of the dimension and of what it measures, concurrent
placements, saving, documents without the field), the geometry unit tests
and property test in `crates/ui/src/sketch_mode.rs` (every kind of
dimension with a placed label is drawn at it with its lines reaching it;
dimension lines follow labels; angles pick their side and widen; automatic
positions unchanged), `crates/ui/tests/sketch_labels.rs` (the editor end to
end for every dimension kind, labels following moved geometry, reset) and
the shell test `dragging_a_dimension_label_in_the_viewport`.

## 13. Transform

Move, rotate and scale the selection, or copies of it (ADR 0014 records
the choices below and the alternatives).

* **The map** is a similarity T(p) = c + s·R(θ)(p − c) + v with s > 0:
  Move is a translation (two clicks: from, to), Rotate a turn θ about a
  centre (clicks: centre, a reference point, the target direction), Scale
  a factor s about a base point (clicks: base, a reference point, where it
  goes: s = |t − c| / |r − c|). The amounts can also be typed. A turn
  within 1e-9 rad of a whole number of quarter turns is taken as exactly
  that, with cos and sin exactly 0 or ±1, so horizontal and vertical stay
  exact. Similarities keep orientation, so arcs keep their start and end,
  tangent sides and ρ; they keep every angle and multiply every length by
  s (any text on similarity geometry, e.g. Coxeter, *Introduction to
  Geometry*, 2nd ed., 1969, ch. 5). That is why constraints among the
  selection survive it.
* **What moves.** The moved set M is the selected curves and points, the
  points of selected curves, and every other curve all of whose points are
  in M. Each moved point goes to T(p) and each moved circle's radius
  becomes s·r. A curve with only some of its points moved (it shares a
  point with the selection) is stretched by the solve below.
* **Constraints** are *internal* (all their entities in M), *external*
  (none) or *boundary*; a derived constraint (mirror, offset, pattern) is
  classified pair by pair, together with its mirror line, direction line
  or centre.
  * External: unchanged.
  * Internal: kept, dimensions of length multiplied by s exactly (angles
    and ρ unchanged). A fix's position moves with its point. Under a whole
    number of quarter turns, horizontal and vertical constraints (lines,
    point pairs, distances) swap for an odd number and stay for an even
    one. Under any other turn, the first horizontal or vertical line (in
    constraint order) loses its constraint, freeing exactly the rotation,
    and every other becomes parallel (same kind) or perpendicular (the
    other kind) to it; a horizontal or vertical distance between points
    that were level the other way (within 1e-7 m) becomes an aligned
    distance; horizontal and vertical point pairs and other such distances
    are removed. Removals are reported.
  * Boundary dimensions are measured anew on the result; one that measures
    1e-7 m or less is removed and reported.
  * Boundary relations are kept ("Keep constraints", the default: the rest
    of the sketch follows, as in a drag) or removed (detached).
* **Solve and landing.** The sketch with the moved geometry placed at T
  and the adjusted constraints (boundary dimensions left out) is solved
  with every moved point pinned at T(p), every moved circle's radius at
  s·r and every moved conic's ρ where it was (temporary fixes, radius and
  ρ dimensions that are not committed). So the selection lands exactly
  where T puts it, and only the rest of the sketch can move to satisfy
  what ties it to the selection. If that does not solve the transform is
  refused, never approximated: "the rest of the sketch holds the
  selection" (a conflict, or no convergence), or "that would collapse a
  curve joined to the selection" (a line from a moved point to an unmoved
  one shrunk to nothing). Boundary dimensions then get their measured
  values. The batch (constraint changes and the solved geometry) goes
  through the checked edit: no new redundancy, still solvable.
* **Copy** adds transformed copies of M with the internal constraints,
  adjusted the same way, and nothing tying them to the rest of the sketch;
  the originals are untouched. A copy onto its original (the identity) is
  refused.
* **Labels** of internal dimensions placed by the user turn and scale with
  them (§12).
* **Tools.** Move, Rotate and Scale share a toolbar button and act on the
  selection made before choosing them. Hovering previews the result (or
  shows why it is refused); the last click commits as one undo step and
  returns to Select with the moved entities (or the copies) selected. The
  sketch dialog holds the options (Copy, Keep constraints) and the typed
  amounts (X and Y offsets, angle, factor), applied with its Apply button
  once the centre or base point is clicked (Move needs none).

## 14. Tests

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
