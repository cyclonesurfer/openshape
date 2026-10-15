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
   loose.
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
and the median of 12 runs interleaved with main's. The limit is what the
test enforces.

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
| Dimple a plate with 16 balls, one boolean each | 82–89 ms | 900 ms |
| Mesh it (173,680 triangles) | 382–406 ms | 4,000 ms |
| Validate it | 2.4–2.6 ms | 40 ms |
| Its mass properties to 1e-9 | 14–17 ms | 150 ms |
| Cut 16 countersunk holes (a shank and a cone each) into a plate, one boolean each | 100–116 ms | 1,200 ms |
| Mesh it (9,612 triangles) | 31–37 ms | 300 ms |
| Validate it | 3.6–3.8 ms | 50 ms |
| Its mass properties to 1e-9 | 22–24 ms | 250 ms |
| Cut 16 countersunk holes leaning 10° into a plate in one boolean (ellipses) | 20–26 ms | 350 ms |
| The same holes one at a time (each cut's cones carry the ellipses of the cuts before) | ONE_SINK | 2,500 ms |
| Mesh it (12,956 triangles) | 44–59 ms | 600 ms |
| Validate it | 4.0–4.9 ms | 60 ms |
| Its mass properties to 1e-9 | 17–20 ms | 200 ms |
| Flatten 16 pins of 5° draft along their axes in one boolean (hyperbolas) | 17–24 ms | 350 ms |
| The same pins one at a time (each cut's cones carry the hyperbolas of the cuts before) | ONE_FLAT | 2,500 ms |
| Mesh it (9,868 triangles) | 26–35 ms | 350 ms |
| Validate it | 3.7–4.5 ms | 60 ms |
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

One more checks a cut against four times the holes before it (bound: 12,
three times the factor):

| Operation | Measured growth |
|---|---|
| Cut the 4th or the 16th countersunk hole leaning 10° into a plate | GROWTH |

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
| Cut 16 countersunk holes leaning 10° into a plate in one boolean | 41.7 ms | 25.7 ms |
| Flatten 16 drafted pins in one boolean | 33.4 ms | 23.7 ms |
| Cut 16 countersunk holes, one boolean each | 273.7 ms | 116.2 ms |
| Dimple a plate with 16 balls, one boolean each | 124.5 ms | 89.1 ms |
| Validate the tilted countersunk plate | 13.8 ms | 4.9 ms |
| Validate the flattened pins' plate | 9.6 ms | 4.5 ms |
| Validate the countersunk plate | 11.8 ms | 3.8 ms |
| Validate the dimpled plate | 4.2 ms | 2.6 ms |
| Drill 64 round holes into a plate, one boolean each | DRILL_BASE | DRILL_NEW |

Cutting the tilted holes or flattening the pins one at a time was
refused on main (a NURBS edge on a curved face). DRILL_NOTE

## Known slow spots

Measured by wall clock, not yet fixed; each is a candidate for a later
pass.

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
  time (the products alone 11%), whose cuts take TWICE as long each as
  the straight countersinks'. The products could be subdivided with the
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
