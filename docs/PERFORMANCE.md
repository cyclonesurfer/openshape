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
   and at four times that size; the ratio is bounded (8, where an
   n log n operation grows about 4.5 times and a quadratic one 16 times).
   A ratio holds on any machine, so these checks catch a change of
   complexity even where the budgets are loose.
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
machine. The limit is what the test enforces.

| Operation | Measured | Limit |
|---|---|---|
| Extrude a 1024-sided star prism (build and validate) | 12 ms | 120 ms |
| Mesh it | 9 ms | 100 ms |
| Validate it | 9.5 ms | 100 ms |
| Its mass properties to 1e-9 | 70 ms | 700 ms |
| Mesh a revolved ring (115,352 triangles) | 270–340 ms | 3,000 ms |
| Extrude a plate with 100 round holes | 12.5–23 ms | 150 ms |
| Mesh it | 135–220 ms | 1,500 ms |
| Union or subtract two 256-sided star prisms | 40–60 ms | 500 ms |
| Cut 64 holes into a box, one boolean each | 240–330 ms | 3,000 ms |
| Sketch regions of 122 curves (a 41 by 41 grid and 40 circles) | 105–115 ms | 1,200 ms |
| Regenerate the demo document | 36.5–41 ms | 400 ms |

Scaling, each size against eight times it (bound: 24, three times the
factor):

| Operation | Measured growth |
|---|---|
| Extrude a star prism, 256 to 2048 sides | 9.7–10.3 times |
| Union of two regular prisms, 128 to 1024 sides | 11.3–11.5 times |
| Mesh a plate, 16 to 121 holes (7.6 times) | 7.4–7.6 times |

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

## Known slow spots

Measured, not yet fixed; each is a candidate for the next pass.

* **The exact incircle test's fallback.** The triangulation library's
  Delaunay test (Shewchuk's incircle, from the `robust` crate) falls back
  to exact arithmetic whenever four points are on one circle, and the
  fallback zeroes a 9 KB buffer on every call. A cylinder face's points in
  its parameter plane form rows of rectangles, all such cases, so meshing
  a plate with 100 holes spends about 60% of its time there, mostly in the
  zeroing; booleans spend about a tenth. The fix is a smaller buffer in
  that fallback (in the library, upstream or by a patched copy).
* **One triangulation per boolean piece.** Each piece of a split face gets
  a point well inside it from a constrained triangulation, about 100,000
  instructions a piece, mostly fixed cost; for two 512-sided star prisms
  that is a quarter of the boolean.
* **Validation meshes the result** to check each shell's orientation by
  its volume, a quarter of a boolean's time.
* **Holes mesh one at a time:** a plate with 1,600 holes takes about 2.5 s
  to mesh, 1.5 ms a hole, most of it the incircle fallback above.
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
guards.
