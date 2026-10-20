# Performance: budgets, guards and limits

Efficiency is part of the quality bar, alongside correctness: an operation
that is right but takes minutes is not done. This document records how that
is kept true and measured.

## Rules

1. **Every user-facing operation has a budget**: a time at a reference size
   on the reference machine, and a release-build test that fails at about
   ten times the measured time. Ordinary noise and slower machines pass; an
   accidental quadratic loop, a lost index or a debug path left on fails.
2. **Scaling is tested on its own.** The same operation is timed at a size
   and at eight times that size, the runs interleaved; the ratio is
   bounded at 24, three times the factor (an n log n operation grows about
   9 to 10 times, a quadratic one 64 times). A ratio holds on any machine,
   so these checks catch a change of complexity even where the budgets are
   loose. A cost that is a small part of its operation's time can grow
   faster without failing either; where profiling finds one, a unit test
   counts its work instead (eighth pass, below).
3. **Hot code uses the standard data structures**: spatial indexes
   (`crates/kernel/src/bvh.rs`, a bounding-volume hierarchy; uniform grids)
   instead of loops over all pairs, bulk-loaded triangulations, exact
   predicates instead of tolerance bands where an exact answer exists, and
   incremental updates instead of rebuilding state in rounds.
4. **Every operation whose work grows with its input has a work limit**
   that refuses, with a message, rather than running for minutes or
   exhausting memory (table below). Limits sit far above any model a person
   builds; a refusal names the limit.
5. **Profile before recording numbers, and record what was measured.**
   Numbers here come from release builds, best of several runs, on the
   machine named with them.
6. **CI enforces all of it**: the time guards run in release on every push
   (`cargo test --release --test timing` for the kernel, the sketch solver
   and the app), each job and long step has a timeout, and a hang fails in
   minutes, not after GitHub's six-hour default.

## Reference machine

The development container: a 4-core Intel Xeon at 2.1 GHz (cloud VM), Rust
stable, release profile. CI's runners are comparable; the limits leave
room for them.

## Budgets

Release builds, best of 2–3 runs, measured 2026-10-13 on the reference
machine (the rows that mesh curved faces, and the new rows from the
tilted countersunk holes to the creased revolve's validation, on
2026-10-15, after the fifth pass). The booleans and validations of the
dimpled, countersunk, tilted and flattened plates, and the cuts one at a
time, were measured again after the sixth pass (2026-10-15): the lowest
and the median of 12 runs interleaved with main's. The union with a
plate of flush pins was measured on 2026-10-16, after the seventh pass:
the range of six runs (three of every guard, three of the booleans' and
the scaling tests alone), each test reporting its best repetition; in
those runs rows already in the table often measured above their ranges
(the flattened pins' plate validated in 4.5–8.8 ms against 3.1–3.8), so
this row may read high against the others. The two rows of a face
lying or crossing along one line, and the scaling row of the first,
were measured on 2026-10-16, after the eighth pass: the range of three
runs of the booleans' and the scaling tests, each the best of three
repetitions. The two rows of notch tips and spikes across the line, and
their scaling rows (design note 0012 §21), were measured the same way on
2026-10-19 on a quiet machine. The limit is what the test enforces.

| Operation | Measured | Limit |
|---|---|---|
| Extrude a 1024-sided star prism (build and validate) | 12 ms | 120 ms |
| Mesh it | 8 ms | 100 ms |
| Validate it | 8.6 ms | 100 ms |
| Its mass properties to 1e-9 | 70 ms | 700 ms |
| Mesh a revolved ring (85,368 triangles) | 222–274 ms | 3,000 ms |
| Mesh a unit ball (126,656 triangles) | 319–374 ms | 3,500 ms |
| Extrude a plate with 100 round holes | 12.5–23 ms | 150 ms |
| Mesh it | 198–261 ms | 2,500 ms |
| Union or subtract two 256-sided star prisms | 38.5–44 ms | 500 ms |
| Cut 64 holes into a box, one boolean each | 225–265 ms | 3,000 ms |
| Drill 64 round holes into a plate, one boolean each | 860–1,050 ms | 10,000 ms |
| Intersect or subtract a plate with 100 round holes and a box across them | 52–55 ms | 600 ms |
| Unite a box with a corner of a plate of 400 flush square pins (the result has 4,012 faces; each place of it checked against the faces near it) | 95–113 ms | 1,000 ms |
| Intersect a notched block with a prism whose face's edge of 256 segments lies on the block's top beside the line where they cross (the strip check searches the whole strip) | 11.7–13.1 ms | 120 ms |
| Intersect a block with a prism whose face crosses its top 1,024 times along one line (513 solids) | 65–71 ms | 700 ms |
| Intersect a block with a prism whose face has 512 notch tips across the line where it crosses the block's top, each within the tolerance of the top (refused at the first; each tip sorted and tested against the stretches along the line) | 9.95–10.5 ms | 100 ms |
| The same with 512 spikes across the line, each one's flap built (the result has 1,031 faces) | 42.9–44.3 ms | 450 ms |
| Dimple a plate with 16 balls, one boolean each | 73–89 ms | 900 ms |
| Mesh it (173,680 triangles) | 382–406 ms | 4,000 ms |
| Validate it | 2.2–2.5 ms | 40 ms |
| Its mass properties to 1e-9 | 14–17 ms | 150 ms |
| Cut 16 countersunk holes (a shank and a cone each) into a plate, one boolean each | 96–102 ms | 1,200 ms |
| Mesh it (9,612 triangles) | 31–37 ms | 300 ms |
| Validate it | 3.4–4.0 ms | 50 ms |
| Its mass properties to 1e-9 | 22–24 ms | 250 ms |
| Cut 16 countersunk holes leaning 10° into a plate in one boolean (ellipses) | 21–24 ms | 350 ms |
| The same holes one at a time (each cut's cones carry the ellipses of the cuts before) | 176–197 ms | 2,500 ms |
| Mesh it (12,956 triangles) | 44–59 ms | 600 ms |
| Validate it | 3.4–4.2 ms | 60 ms |
| Its mass properties to 1e-9 | 17–20 ms | 200 ms |
| Flatten 16 pins of 5° draft along their axes in one boolean (hyperbolas) | 18–22 ms | 350 ms |
| The same pins one at a time (each cut's cones carry the hyperbolas of the cuts before) | 175–197 ms | 2,500 ms |
| Mesh it (9,868 triangles) | 26–35 ms | 350 ms |
| Validate it | 3.1–3.8 ms | 60 ms |
| Its mass properties to 1e-9 | 13 ms | 150 ms |
| Revolve a rational cubic bending nearly to a crease (validation meshes it) | 12–15 ms | 150 ms |
| Mesh it (108,486 triangles) | 446–482 ms | 5,000 ms |
| Mesh it at 1e-3 and 5° (92,264 triangles) | 400–407 ms | 4,000 ms |
| Validate it | 13–16 ms | 150 ms |
| Sketch regions of 122 curves (a 41 by 41 grid and 40 circles) | 100–115 ms | 1,200 ms |
| Regenerate the demo document | 36.5–41 ms | 400 ms |

Scaling, each size against eight times it (bound: 24, three times the
factor):

| Operation | Measured growth |
|---|---|
| Extrude a star prism, 256 to 2048 sides | 9.7–12.2 times |
| Union of two regular prisms, 128 to 1024 sides | 11.3–12.6 times |
| Mesh a plate, 16 to 121 holes (7.6 times) | 6.8–8.2 times |
| Cut a box across a plate's round holes, 16 to 121 holes (7.6 times) | 7.8–8.4 times |
| Mesh a ball, 13,964 to 105,556 triangles (7.6 times) | 8.6–9.6 times |
| Unite a box with a corner of a plate of flush pins, 49 to 400 pins (8.2 times) | 8.4–9.2 times |
| Intersect a notched block with a prism whose lying edge has 128 to 1,024 segments | 7.7–7.8 times |
| Intersect a block with a prism refused at notch tips across the line, 64 to 512 tips | 8.0–8.1 times |
| Intersect a block with a prism with spikes across the line, 64 to 512 spikes | 7.8–8.6 times |

One more checks the cut of a plate's 4th hole against that of its 16th,
four times the holes counting the one cut (bound: 12, three times the
factor):

| Operation | Measured growth |
|---|---|
| Cut the 4th or the 16th countersunk hole leaning 10° into a plate | 3.1–3.6 times |

The union is scaled on regular prisms, whose result grows with their
size. Two star prisms overlap in a band where many faces of each meet
many of the other, so their result itself grows faster than their size
(5.2 times the faces for four times the sides).

Tests: `crates/kernel/tests/timing.rs`, `crates/ui/tests/timing.rs`.

The sketch solver's own guards (design note 0005 §8) hold a 50-entity solve
to 8 ms and a 300-entity solve, or one drag step of it, to 50 ms
(`crates/sketch-solver/tests/timing.rs`); until now CI ran them only as
debug builds, where the limits are not applied.

## Work limits

| Operation | Limit | Where |
|---|---|---|
| Boolean | 2,000,000 pairs of faces that meet; 8,000,000 pairs of edges near each other; 4,000,000 candidate points; 8,000,000 edge segments | `crates/kernel/src/boolean/mod.rs` (`LIMITS`) |
| Boolean, the strip checks of design note 0012 §20 | 4,194,304 (2²²) steps for all of one boolean's strips | `crates/kernel/src/boolean/mod.rs` (`LIMITS.strip_steps`, set from `imprint::STRIP_STEPS`) |
| Meshing | 2,000,000 vertices a face; edge samples grow at most eightfold (plus 100,000) while refining | `crates/kernel/src/tess.rs` |
| Sketch regions | 40,000,000 intersection steps; curves of degree at most 24 and 10,000 points | `crates/kernel/src/arrangement` |
| Mass properties | 30 levels of bisection of a boundary piece | `crates/kernel/src/mass.rs` |
| Sketch solver | 100 Newton iterations a component (default), 32 restarts | `crates/sketch-solver` |

A limit is a refusal, never a truncated answer: the kernel reports
`TooLarge` (booleans) or a meshing or measurement error, and the feature
shows it.

## What changed, 2026-10-13

### First pass

Profiling (callgrind on release builds with line tables) found these, now
fixed; times are on the reference machine.

* **Booleans: indexes instead of all-pairs loops.** A bounding-volume
  hierarchy (`crates/kernel/src/bvh.rs`) finds the faces that meet, the
  edges near each other, the points near a segment and each face's edges
  near a plane. Inside-or-outside is decided by exact ray parity
  (Shewchuk's orientation predicate on a triangulation of the other body,
  pruned by the same tree) instead of a winding sum over every triangle.
  Collinear edges are joined incrementally instead of in rounds over the
  whole result. The union of two 1024-sided star prisms went from 5.4 s
  to under 0.6 s; the boolean fuzz gives exactly the same outcomes as
  before (6,000 random pairs, 11 refusals; 1,400 tolerance-scale motions,
  2,180 refusals).
* **Classification triangles without slivers.** A convex face was split
  into a fan for the ray tests, whose slivers all reach across the face,
  so a ray near a many-sided face met most of their boxes: the union of
  two 1024-sided regular prisms spent 44% of its time there, growing 44
  times for eight times the sides. Convex faces are now split by halving
  the boundary level by level; classification fell from 598 to 148
  million instructions there, and the union from 157 ms to 110 ms.
* **Triangulations are bulk loaded.** Inserting a polygon's points one at
  a time around its boundary costs a number of Delaunay tests growing with
  the square of its size when the points lie on a circle (any regular
  polygon, any hole). Every polygon is now triangulated by the library's
  constrained bulk load, O(n log n).
* **The same output on every run.** The triangulation library's bulk
  load inserts the points it defers in the iteration order of a randomly
  seeded hash set, so a face whose boundary points share a circle (an
  annulus, a plate with holes) was meshed differently from run to run.
  The workspace uses a copy of the library patched to defer them in order
  (ADR 0022, `vendor/README.md`), and a test meshes and combines such
  models several times over with fresh seeds and requires identical
  output.
* **Hash tables use the Fx hasher** (`rustc-hash`), not the standard
  library's SipHash: SipHash was a tenth of a boolean's time. Fx is not
  seeded at random either, so no result can depend on a seed; a lint
  (`crates/kernel/clippy.toml`) keeps the standard tables out of the
  kernel. A 512-sided star union fell from 1,278 to 1,031 million
  instructions.
* **Profiles: no all-pairs loops.** Joining curve ends into vertices used a
  linear search per end, the check that a profile's curves touch only at
  their vertices tested every pair of curves, and nesting the loops tested
  every loop against every other. A grid and two tree queries replace
  them: extruding an 8192-sided polygon went from 315 ms to 92 ms, and a
  plate with 1,600 holes from 1.35 s to 0.21 s.

### Second pass

Each boolean's stages were timed by wall clock (instrumented release
builds) on the 256-sided unions and the 64 cuts above.

* **Faces clear of the other body need no test.** A face farther than the
  tolerance from the other body's box is untouched by it and wholly
  outside it, so its piece needs no point inside it and no parity test.
  Over the 64 cuts, where all but two faces of the growing body are clear
  of each pin, finding piece points fell from 45 ms to 25 ms.
* **Convex pieces get their point without a triangulation.** A piece is
  classified at a point well inside it, the incentre of the best triangle
  of a constrained triangulation; a convex piece (decided exactly) now
  uses the better of its centroid, which is at least two thirds as far
  from the boundary as the best possible point, and the incentres of its
  balanced triangles. For the 256-sided unions, finding piece points fell
  from 11.8 ms to 7.5 ms a boolean, and the union from 48 ms to 43 ms. The
  boolean fuzz at twenty times the default size refuses exactly the same
  2,191 operations for the same reasons as before; only the place a
  sliver refusal names moved, since the point inside the piece moved.
* **Convex flat faces are meshed without a search.** A flat face bounded
  by one strictly convex loop is split by halving its boundary's runs (the
  triangles of the parity test) instead of by a constrained Delaunay
  triangulation: meshing the 1024-sided prism, whose sides are 2,048
  rectangles, fell from 9 ms to 8 ms.
* **Not carried: a patched incircle test.** The `robust` crate's exact
  incircle test zeroes a 9 KB buffer on every call of its exact stages.
  An instruction count put that at nearly half of meshing a plate with 100
  holes, but a patched copy that sizes the buffer to the stage gave the
  same answers bit for bit (36 million calls against the published crate,
  16 million of them reaching the exact stages and 300,000 the last of
  them) and saved only 3–6% of the wall-clock time: zeroing memory costs
  far less per instruction than the arithmetic around it. Too little to
  carry a patched copy for (ADR 0022).

### Third pass: booleans with curved faces

The new curved stages (design note 0012 §11) were timed on one round
hole drilled into a plate extruded with 99, by callgrind and then by wall
clock: 211 ms at first.

* **A lone shell needs no volume.** Every curved result's shells were
  measured by the exact mass properties to order them into solids and
  voids; a result with one shell is one solid, and validation already
  checks that it is not inside out. With several, each shell's volume now
  comes from the coarse mesh validation uses, with a bound on its error
  (twice its area times the mesh's distance from the surfaces), and the
  exact mass properties only where that bound could reach the sign.
* **Chords meeting at a node are halved at once.** A curved segment's
  chords are halved where their bands could hide a contact; two segments
  meeting at a point are always within their bands there, so the chords
  at every node were halved one round at a time, about eleven rounds, and
  every round rebuilt the whole graph and tested every chord. Those
  halvings are now made first, all at once, and each round tests only the
  chords new in it (a pair that stood unchanged through a round was
  tested in it).
* **Large pieces are probed, not triangulated.** A piece with more than
  64 boundary points gets its interior point from probes: from the middle
  of a few sides, inward to the first side met, halfway, each checked by
  an exact winding number and measured against every side. The
  triangulation of the face with 100 holes (5,600 chords) spent 22 ms a
  face, mostly in exact incircle tests on points of one circle; the probes
  take about 1 ms. The triangulation remains for small pieces and
  wherever the probes find no point far enough inside.

The single hole fell from 211 ms to 40 ms: 18 ms in splitting the two faces
with holes into regions, 15 ms validating the result (which meshes it;
see the slow spots below), every other stage under 3 ms.

### Fourth pass: spheres and cones (item 3a)

The time guards above for balls, cones and their meshes were added with
booleans of spheres and cones (design note 0012 §13), and found two
things.

* **Validation measured every body exactly for a while.** Validation
  now decides a shell's orientation as the booleans do: from the coarse
  mesh's volume where its error bound cannot reach the sign, else from
  the shell's mass properties (a thin ring a ball keeps outside a cone
  meshed to no volume at all and was refused as inside out). The bound
  was twice the area of the whole mesh times its distance from the
  surfaces, which for a thin flat body exceeds its volume: validating the
  1024-sided prism went from 10 ms to 87 ms, its extrude from 13 ms to
  90 ms, and the prism's scaling check failed (65 times for 8). The mesh
  of a flat face lies in its plane, so only curved faces can put the mesh
  off the body: the bound now counts their area alone (plus the
  rounding of the sum), and every guard is back to its time on main, by
  interleaved runs of both builds.
* **Points added to a curved face's mesh started their search far
  away.** Refinement inserted each new point into the triangulation with
  no hint, so the library walked to it from the last point inserted,
  across the face: about √n triangles a point. The new scaling check on
  a ball found its time growing 14.5 times for 7.8 times the triangles.
  Each point now starts from a corner of the triangle it splits, and the
  growth is 11–12 times (n log n). A unit ball's default mesh fell from
  364–411 ms to 220–249 ms and the revolved ring's from 312–387 ms to
  204–247 ms; meshes of about a million triangles take a third of the
  time. The meshes keep their triangle and point counts but are not the
  same as before: the coarse ball has the same triangles in another
  order, and finer meshes of the ball and the ring differ. Inferred from
  that: each round adds its points in the order of the triangulation's
  faces, which the library numbers differently when it reaches a point
  from another side, and points on one circle have several equally good
  triangulations. Each mesh is the same on every run.

## What changed, 2026-10-15

### Fifth pass: meshes measured against their tolerance (item 3b)

Meshing of curved faces was rebuilt to hold the chord tolerance much
more closely (design note 0002 §4): each triangle is measured at
thirteen points, with the surface's derivatives at each, and the quartic
through them is bounded. Between the points this is an estimate, not a
proof: a surface that bends faster than a quartic there can stray
further. On every probe case the worst distance found is now within the
tolerance, where main measured the surface's height over the triangle's
plane at its centroid and the middles of its interior edges (up to four
points) and strayed past it (5% past on a ball, 2.4 times the tolerance
on a nearly closed torus). The first version took 692 and 726 ms in two
runs for a unit ball's default mesh, 3.6 times main's 193 ms in the
first of them. Profiled with callgrind and the stages timed in
place:

* **A hash table of the edges measured** grew to over a hundred
  thousand entries a face on the ball, tens of megabytes. It is now a
  direct-mapped cache of at most 4,096 slots (754 kB), one an edge by a
  hash of its corners, and the triangles of a pass are judged in the
  order found about each new point, so the two beside an edge come close
  together. Counted on the ball's default mesh (2026-10-15), 42% of the
  lookups hit, of the 49% that can (each edge is looked up from both
  sides), and measurements repeated come to an eighth of the edges;
  judged in the order of their corners instead, 32% hit and repeats come
  to a third (the ring, the S-spline and the creased
  revolve likewise: 41–45% against 32–36%).
* **The longest edge is measured first** and the first failure decides,
  so a failing triangle is measured only up to it.
* **No inverse cosines:** the turn of the normal comes from above as a
  tangent (0002 §4); points on lines and circles take one `sin_cos`.
* **One evaluation a vertex:** a vertex's point, normal and whether it is
  singular come from one call of the surface's derivatives.
* **NURBS basis functions in arrays.** Meshing a rational cubic revolved
  nearly to a crease spent 59% of its 6.42 billion instructions in the
  general basis routine (Algorithm A2.3, inclusive), which builds and
  moves small vectors on every call. Values and first derivatives, all
  that evaluation for meshing asks, now come from fixed arrays with the
  same arithmetic in the same order, equal bit for bit (a property test
  checks degrees 1 to 7 on random knot vectors): 3.18 billion
  instructions, half as many.

Against main, best of five runs interleaved with two of main's, on the
reference machine (2026-10-15):

| Mesh | main | now |
|---|---|---|
| Unit ball, defaults | 104,988 triangles, 179 ms, 1.052 tolerances off | 126,656 triangles, 377 ms, within (0.9998) |
| Ring (torus, radii 1 and 0.3), defaults | 228,880, 364 ms, 1.010 | 174,756, 504 ms, 0.9998 |
| Turned S-spline, 10⁻³ | 30,400, 52 ms, 0.9994 | 24,584, 101 ms, 0.9935 |
| Creased rational revolve, 10⁻³ and 5° | 85,902, 198 ms, 0.9917 | 92,264, 337 ms, 0.9887 |

(Worst distance from closest points on a grid of sixths of every
triangle, in tolerances. The benchmark runs one case at a time; the time
guards above run beside other tests and read higher, e.g. 400–407 ms for
the creased revolve.) On these four bodies a triangle now costs 1.6 to
2.4 times what it did: the measurements and the bound, and the exact
incircle tests that bisecting edges on a regular grid brings (slow spots
below).

### Sixth pass: operands with conic edges (design note 0012 §16)

A second cut of a body whose cones and cylinders carry conic edges
follows each edge's NURBS curve in its face's chart, with a band from the
curve's derivatives (design note 0012 §16.2). Profiled with callgrind
and compared with main by wall clock:

* **Edges near an axis were sampled at 65 points, with an 80-step
  search at each local least.** Validation and the booleans' last check
  ask whether an edge of a turned face passes through its axis
  (`check::through_axis`). Along a line parallel to the axis, as the
  straight edges of a turned cylinder run, the distance is the same but
  for rounding, so many samples were local leasts, each starting its own
  search: 30% of the instructions of building, cutting, validating,
  meshing and measuring the tilted countersink and flattened-pin plates
  were spent there. A line now takes the closed form (its
  least distance from the axis, clamped to its ends), a circle a test of
  its least distance, and a NURBS curve is cleared when the control
  points of each Bézier piece, whose convex hull holds the piece (Piegl
  & Tiller §4.2), all lie farther than the tolerance from the axis; only
  what these cannot clear is sampled as before. Counted on those two
  plates: 2.44 billion instructions in the check before, 2.0 million
  now; validating each once, 216 million before, 76 million now; the 48
  booleans that build and cut them, 6.62 billion before, 4.52 billion
  now.
* **The band of a NURBS course in a plane's chart split the whole curve
  into its Bézier pieces again for every chord**: 9.8% of the
  instructions of cutting the tilted countersunk holes one at a time.
  It now uses the pieces its trace split once, with the same arithmetic:
  1.89 billion instructions before, 1.70 billion now.
* **Binomials from a table.** The products of Bernstein polynomials made
  in the bands' innermost loops computed each binomial as a product; a
  table built by Pascal's rule at compile time, exact, replaces it (not
  measured on its own).
* **Solids with voids:** the check of which outer shell holds each void
  builds the shells' outlines once, and only when there is a void.

Meshing and mass properties are unchanged: the same instruction counts
on the tilted and flattened plates (745 and 746 million for meshing,
369 million for mass properties, before and after). Against main, the
median of 12 runs of each build interleaved, on the reference machine
(2026-10-15):

| Operation | main | now |
|---|---|---|
| Cut 16 countersunk holes leaning 10° into a plate in one boolean | 41.1 ms | 24.3 ms |
| Flatten 16 drafted pins in one boolean | 33.2 ms | 22.1 ms |
| Cut 16 countersunk holes, one boolean each | 262.3 ms | 101.6 ms |
| Dimple a plate with 16 balls, one boolean each | 122.2 ms | 89.4 ms |
| Validate the tilted countersunk plate | 14.3 ms | 4.2 ms |
| Validate the flattened pins' plate | 9.6 ms | 3.8 ms |
| Validate the countersunk plate | 11.4 ms | 4.0 ms |
| Validate the dimpled plate | 4.1 ms | 2.5 ms |
| Drill 64 round holes into a plate, one boolean each | 813.6 ms | 745.8 ms |

Cutting the tilted holes or flattening the pins one at a time was
refused on main (a NURBS edge on a curved face). Drilling the round
holes runs as many instructions as on main: 7.851 billion under
callgrind for the test that drills them and cuts a box across the
plate, against 7.833 billion (0.24% more); its medians swing with the
machine's noise (an earlier set of 12 runs had it 5% slower than main,
this one 8% faster). The guard of the solids' meshes and mass
properties likewise: 660.4 million instructions, against 660.5 million
on main, though some of its medians here differ by up to 8% (the
1024-sided prism's mass properties, 71.0 ms against 65.7 ms; 2% in the
earlier set).

## What changed, 2026-10-16

### Seventh pass: contacts that do not cross (design note 0012 §19)

The slice decides contacts in each face's split and checks every result
for a place within the tolerance of an exact face it is not on
(`assemble::check_apart`). Counted by callgrind against main, one run
each:

| Operation | main | now |
|---|---|---|
| Union of two 256-sided star prisms | 295.5 million | 311.9 million (+5.6%) |
| Cut 64 pins into a box, one boolean each | 2,168 million | 2,344 million (+8.1%) |

Of the difference, the check takes 11.3 million and 95.5 million, and
`build`, which now keeps only the origins that bound each entity, 1.5
million and 18.9 million more than on main.

* **The check was quadratic in a face's loops.** Its first version
  tested each point against a face by walking all the face's sides
  (parity across every loop), so a box across a corner of a plate whose
  square holes each hold a flush pin (islands on the top face) took
  38.1 ms, 257.8 ms and 1,119.8 ms at 100, 400 and 900 pins (found by
  review). Each face's outline is now built when a point first needs it:
  its sides in the plane, a tree of them, and a tree of its loops'
  boxes, a loop of more than 32 sides with a tree of its own; each place
  is tested once, and the faces with a corner at each point and the
  faces along each segment are grouped by a counting sort. From 100 to 900 pins the
  union's instructions now grow 9.7 times (179 million to 1,739
  million) for 9 times the pins, the check's 11.1 times (19.6 million
  to 216.4 million). The 400-pin union is guarded (above), with a
  scaling check.

### Eighth pass: a meeting beside a lying run (design note 0012 §20)

The strip check (§20) searches, for each chain of a face's boundary
lying within the tolerance of the other face's plane, the strip between
the chain and the line where the planes cross. Counted by callgrind
against main, one run each of plain release builds, each operation a
test of its own (building the operands included):

| Operation | main | now |
|---|---|---|
| Union of two 256-sided star prisms | 361.1 million | 365.5 million (+1.2%) |
| Cut 64 pins into a box, one boolean each | 2,345 million | 2,364 million (+0.8%) |

Of the pins' 19.0 million more, the strip checks take 13.7 million: a
face whose boundary has no mark within the tolerance of the other's
plane is skipped (the faces' relation knows it), and the outlines and
the tables of them are built only when a strip needs them. The star's
strip checks take 86 thousand; most of its 4.4 million more is in
classification (27.1 to 30.7 million), whose code differs from main's
only in taking its refusal's text from a constant, but whose walk of the
other body's tree (`Bvh::visit`) the compiler no longer inlines.
Profiling two bodies made for it found three places doing work that
grows with the square of their size, two of them on main:

* **Each strip's search walked every edge near the plane.** The time
  guards' zigzag is a face whose edge of n segments lies on the other
  body's top beside the line where the faces cross, a run with a reflex
  corner at every other vertex: n / 2 + 1 chains. Each chain's search
  took the edges of both faces within its distance of the plane, all n
  of the run's among them: the strip checks took 25.6 million
  instructions at n = 256 and 224 million at 1,024, 8.8 times for four
  times the segments, of which that walk (`Face::near_surface`) grew 15
  times. The walk of each face's tree now keeps to the strip's box as
  well (`Face::near_surface_in`), and each chain takes the pair's
  meetings over its own span from a list sorted once: 16.2 and 66.1
  million (4.1 times). The intersection itself (`Kernel::boolean`) at
  1,024 segments went from 555 to 396 million instructions; by the wall
  clock it takes 49–55 ms (best of three, three runs of it alone).
* **The end checks of a pair's stretches found the faces' lying runs
  again for each stretch** (main): every stretch walked both faces'
  edges near the other's plane, and every pair of the faces' intervals
  on the line was tested. On the time guards' comb, a face crossing a
  block's top 1,024 times along one line, the meeting of each pair, the
  faces' intervals on the line included, took 216 million instructions,
  206 million of them in finding the runs. The runs are now found once
  per pair, sorted along the line with their ends, and each interval of
  one face meets only those of the other within the tolerance of it:
  the meetings take 12.8 million, 2.5 million of them in finding the
  runs.
* **Points were placed in long cycles by walking every side**
  (`planar::regions`, main). Each part of a split face that does not
  touch its outer boundary is placed in the smallest region around its
  first point, tested against each candidate region's cycle by walking
  all its sides; on the comb, the faces' regions have cycles running
  along all the crossings. A cycle of more than 32 sides is now searched
  through the outline of `check_apart` (§19, seventh pass), built when a
  point first needs it: the regions took 174 and 1,373 million
  instructions at 1,024 and 4,096 crossings (7.9 times), and take 128
  and 564 million (4.4 times).

The comb's intersection (`Kernel::boolean`) went from 769 to 518
million instructions at 1,024 crossings and from 7,342 to 3,136 million
at 4,096. By the wall clock, against main in two interleaved runs of an
earlier build of these fixes, it went from 88–91 ms to 64–65 ms and
from 783–789 ms to 374–384 ms; the final build takes 64–65 ms and
391–439 ms (best of three, three runs of each alone). Both bodies are guarded by time
(above), the zigzag with a scaling check; the comb has none, for the
classification below.

The rest of each boolean dilutes these costs in its time, so a time
guard need not fail with one of the fixes undone: with the walk of every
edge near the plane restored, the zigzag's scaling grew 9.7 times
against a bound of 24, and with long cycles walked whole the comb took
71 ms against its limit of 700 (both measured in review). Each fix has a
test of the work it saves instead, in `crates/kernel/src/boolean/mod.rs`,
which fails with it undone: the zigzag's intersection is answered
within 150 strip steps a segment at 256 segments, against a block whose
top also has 256 edges crossing the line beyond the notch, and within
300 a segment at 1,024 against the plain notched block (with the walk
of either face's edges restored, the first takes about 58,000 steps
against 25,000, and with the walk of the edge's restored the second
about n² / 2 more: each is refused), and the
coedges looked at in finding the comb's lying runs, and the sides of
cycles looked at in placing its points, grow at most six times from 256
to 1,024 crossings (with the runs found for each stretch, the coedges grew 15
times; with long cycles walked whole, the sides grew 24 times).

A second review found two more places doing work that grows with the
square of their input. The pieces of a strip were cut by every meeting
beside them, each part by each meeting, and each part carried a side for
every split before it: a segment beside n meetings took about n² steps,
each test of a piece took time in proportion to n, and past about 2,050
meetings the boolean was refused as too large. And which of a face's
lying runs are whole intervals was found by scanning the intervals for
each run. The parts are now cut in order along the line, each meeting
cutting only those not ended before its start, a split replacing the
side it tightens; the intervals are sorted and searched. A unit test
cuts a segment beside 256 and 1,024 meetings within 8 steps a meeting,
each piece keeping at most five sides. The scan has no test of its
work: in review it cost 27 million instructions of a 5.2 billion
boolean built to show it.

## Known slow spots

Measured by wall clock, or counted by callgrind where a bullet gives
instructions, not yet fixed; each is a candidate for a later pass.

* **Delaunay tests on points of one circle.** Every triangulation here is
  a constrained Delaunay triangulation, whose incircle test (Shewchuk's,
  from the `robust` crate) is exact: when four points lie on one circle,
  as the samples of a round hole or the corners of a rectangle do, its
  fast stage cannot decide and it computes in exact arithmetic. Meshing
  a plate with 100 round holes takes about 150 ms, of which the two flat
  faces take 96 ms; the same flat faces with their points moved 1e-9 m
  off their circles triangulate in a third of the time. Meshing grows in
  proportion to the holes (1,600 holes: 2.4 s, 1.5 ms a hole). Flat faces
  need no Delaunay property (any triangulation on their boundary points
  is as good); a triangulation built from orientation tests alone, or a
  cheaper exact incircle stage, would remove most of it.
* **Exact incircle tests on curved faces.** Refinement splits a
  triangle's longest edge at its middle in the parameters, and on a
  regular grid of parameters (a ball's, a torus's) the new point and
  three others often lie on one circle: the Delaunay test then reaches
  its exact stage, which zeroes 26 kB of buffers each call. A unit
  ball's default mesh makes 38,282 such calls (main, adding centroids,
  364): 18% of its instructions, 12% of all of them in clearing memory
  alone. Adding only centroids instead took more triangles (2026-10-15:
  24% more on the S-spline at 10⁻³, 31% on the ring, 5% on the ball,
  about as many on the creased revolve). A cheaper exact stage would
  remove most of it, as for the flat faces above.
* **Sines and cosines** in evaluating a ball's surface take 14% of its
  default mesh's instructions (a torus's points take two `sin_cos` each;
  not profiled).
* **Curved segments carry many chords near their ends.** The halvings at
  a node (third pass) go on until the chords' bands are within a
  sixteenth of the tolerance, about eleven halvings at each end of each
  arc, so a face with 100 round holes is split on 5,600 chords where
  1,600 would follow the circles closely enough elsewhere: 9 ms a face,
  half of drilling a hole into a plate with 99. Two segments leaving a
  point at a wide angle cannot hide a contact near it once each stays
  within a cone about its first chord, which a bound on the segments'
  turning would show; that needs the turning of each course in its
  chart (design note 0012 §11).
* **Bands of NURBS courses on cones and cylinders** (design note 0012
  §16.2) are made afresh for every chord, from products of Bernstein
  polynomials of the curve's piece clipped to the chord's stretch: 31%
  of the instructions of cutting the tilted countersunk holes one at a
  time (the products alone 11%), whose cuts take about twice as long
  each as the straight countersinks' (197 ms for 16 against 102 ms,
  medians). The products could be subdivided with the
  chords (de Casteljau's algorithm on their coefficients) instead of
  made again.
* **Validation meshes the result** to check each shell's orientation by
  its volume, and to check that every face can be meshed: about a fifth
  of a 256-sided union (8 ms of 43) and a quarter of each of the 64 cuts,
  most of it the two faces with many holes, for the reason above.
* **Points inside pieces that are not convex** still come from a
  triangulation: 7.5 ms of a 256-sided union.
* **Every boolean reads both whole bodies.** Converting and indexing each
  operand, triangulating its faces for the parity test, splitting every
  face and validating the result take time in proportion to the bodies,
  not to where they meet: each of the 64 cuts takes 3.5–4 ms on a body of
  up to 262 faces, most of it outside the pin's neighbourhood.
* **Faces grouped by surface with a linear search.** Each boolean
  groups a body's faces by their surface (`boolean::shape::shape`),
  looking each one up among the groups found so far, so the work grows
  with the faces times the distinct surfaces: from 100 to 900 flush pins
  it grows 14.5 times for 9 times the pins (5.7 million to 82.8 million
  instructions, counted by callgrind), 4.8% of the 900-pin union (on
  main too). The candidate, not built or measured: an index of the
  groups by their surfaces' rounded parameters. To keep today's
  grouping it needs cells at least as wide as the matching tolerance,
  keys that ignore a normal's sign and where a point lies along an
  axis, and a search of the neighbouring cells for the earliest group
  that matches.
* **Classification casts a ray from each piece.** Each piece of a split
  face is classified by a ray from a point inside it to the other body
  (`classify::Boundary::classify`), counting the faces it crosses with
  a tree over them, so where a line crosses the other body's face many
  times, each of the pieces along it crosses the others: on the comb
  (eighth pass), 93.7 million instructions at 1,024 crossings and 1,073
  million at 4,096, 11.4 times for four times the crossings, 34% of the
  larger intersection (`Kernel::boolean`). Pieces that share an edge where the bodies do not
  meet lie on the same side of the other body, so a classification
  carried from piece to piece across such edges, with a ray only for the
  first piece of each group, would make it linear; it needs a design
  note first (which edges carry a side, and how a contact edge stops
  it). The comb's guard has no scaling check until then.
* **Meshing to a very fine tolerance can hit the vertex limit:** a full
  ring meshed to 1 µm needs more than 2,000,000 vertices on one face and is
  refused.

## Profiling

```sh
CARGO_PROFILE_RELEASE_DEBUG=line-tables-only cargo test --release -p opendatum-kernel --test timing --no-run
valgrind --tool=callgrind --callgrind-out-file=cg.out target/release/deps/timing-<hash> booleans --nocapture
callgrind_annotate --inclusive=yes cg.out | head -60
```

Callgrind counts instructions, so its numbers are stable from run to run
and compare well between two builds; wall-clock numbers come from the time
guards. An instruction is not a unit of time: zeroing memory runs many
instructions a cycle, and a mispredicted branch costs many cycles, so a
hotspot in an instruction count is confirmed by the wall clock (the best
of several interleaved runs of the two builds, or the stage timed in
place) before it is acted on or written down.
