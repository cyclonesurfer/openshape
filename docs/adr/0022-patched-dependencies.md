# ADR 0022 — Patched copies of dependencies

Status: accepted, 2026-10-13.

## Context

The kernel triangulates polygons with the `spade` crate (constrained
Delaunay triangulation, MIT OR Apache-2.0). Its bulk load, which the
kernel needs because inserting a polygon's points one at a time costs a
number of Delaunay tests growing with the square of their count when they
lie on a circle (docs/PERFORMANCE.md), keeps the vertices and edges it
defers in a `std::collections::HashSet` and inserts them in that set's
iteration order. The standard library seeds that order at random per
process, and points on a common circle (an annulus, a plate with holes, a
rectangle of samples) triangulate differently in different orders. So
the same model meshed differently from run to run, and since booleans
triangulate their pieces, a boolean could in principle decide differently
too. The project requires the same output for the same input on every run
and machine (CLAUDE.md).

Found 2026-10-13 by meshing one revolve twice in one process; spade 2.15.1
is the latest release.

## Options

* **A. Use only spade's unconstrained bulk load** (whose deferred points
  are kept in a list) and add the boundary edges as constraints one at a
  time. Deterministic, but on a boundary with exactly collinear runs of
  samples (a straight edge sampled finely, as in the random revolve test)
  spade's constraint insertion reported a crossing that is not there.
  Rejected after that test failed.
* **B. Insert points one at a time in a spatially sorted, randomized
  order** (BRIO with a Hilbert sort, as CGAL does), seeded deterministically.
  Deterministic and expected O(n log n), but a different, slower path than
  the bulk load the kernel has been tested on, and more code of our own.
* **C. A patched copy of spade** (`vendor/spade`, used through
  `[patch.crates-io]`) whose only change is ordered sets for the deferred
  vertices and edges. The algorithm and everything the kernel has been
  tested on stay as they are; the triangulation becomes a function of the
  input alone.
* **D. Wait for an upstream fix.** Leaves the defect in place meanwhile.

## Decision

**C**, with the change offered upstream when the project chooses to (that
is an outward step, for the maintainers of this project to take). Rules
for any patched copy:

* The copy is the published crate's source, unmodified except for the
  listed patches, with its licence files. `vendor/README.md` lists every
  patch, its reason and the test that shows it.
* The copy is excluded from the workspace (it is linted, formatted and
  tested as the upstream crate is, not by our rules), and cargo-deny still
  checks its licence.
* On each upgrade of the crate, the copy is replaced by the new release and
  the patches reapplied, or dropped if upstream has fixed the defect.
* A golden test fails without the patch:
  `meshes_and_booleans_are_the_same_every_time` in
  `crates/kernel/tests/golden.rs` meshes and combines models several times
  with fresh hash seeds and requires identical output.

## Consequences

The repository carries about 1.8 MB of third-party source (`spade`'s
`src` and the images its documentation includes). The same mechanism is
available for the next known defect of this kind: the exact incircle
predicate of the `robust` crate zeroes a 9 KB buffer on every call of its
fallback (docs/PERFORMANCE.md, known slow spots).

Update, 2026-10-13: measured by wall clock, a patched copy of `robust`
without that zeroing saved only 3–6% of meshing, so it is not carried
(docs/PERFORMANCE.md, second pass). A patch is carried for a defect or a
large measured gain, not for a few percent.
