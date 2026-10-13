# Design note 0012 — Boolean operations on B-rep solids (Phase 3)

Status: design, before implementation. Builds on the intersections of
design note 0011.

## 1. Problem

Union, intersection and difference of two valid closed solids (and the
features built on them: extrude and revolve that add to, remove from or
intersect existing parts, holes, patterns, mirror), producing a valid
solid or a clear error. The result must pass the kernel's validator
(closed loops, manifold edge use, edges on their surfaces within their
tolerance, Euler–Poincaré, outward orientation), and must be the right
solid: the measured volume identities of §7 hold.

Booleans are where most kernel failures happen in practice. The causes
are well understood: decisions about the same geometric fact (does this
vertex lie on that face?) taken twice with different results; missed
intersection branches; tangential and coincident configurations; and
tolerances that do not compose. The design is organised around avoiding
each of these, and around refusing loudly where it cannot.

## 2. Literature

* A. A. G. Requicha & H. B. Voelcker, "Boolean operations in solid
  modeling: boundary evaluation and merging algorithms", *Proc. IEEE*
  73(1), 1985: boundary evaluation as split, classify, select, merge.
* M. Mäntylä, "Boolean operations of 2-manifolds through vertex
  neighborhood classification", *ACM TOG* 5(1), 1986, and *An
  Introduction to Solid Modeling*, Computer Science Press, 1988:
  classification of the neighbourhoods of intersection vertices and
  edges, including the on/on (coincident) cases.
* D. H. Laidlaw, W. B. Trumbore & J. F. Hughes, "Constructive solid
  geometry for polyhedral objects", *SIGGRAPH* 1986: face subdivision
  and inside/outside/same/opposite classification with propagation to
  neighbouring faces.
* C. M. Hoffmann, *Geometric and Solid Modeling: An Introduction*,
  Morgan Kaufmann, 1989, ch. 4: robustness — numerical decisions must be
  made once and their consequences derived symbolically, or inconsistent
  topology follows.
* I. Stroud, *Boundary Representation Modelling Techniques*, Springer,
  2006: the stages of a B-rep boolean and the data each stage needs.
* D. J. Jackson, "Boundary representation modelling with local
  tolerances", *Proc. 3rd ACM Symposium on Solid Modeling*, 1995:
  tolerant vertices and edges, so that approximate intersection curves
  can be part of a valid model.
* J. R. Shewchuk, "Adaptive precision floating-point arithmetic and fast
  robust geometric predicates", *DCG* 18, 1997 (already used by the
  profile and tessellation code).
* A. van Oosterom & J. Strackee, "The solid angle of a plane triangle",
  *IEEE Trans. Biomedical Engineering* BME-30(2), 1983, and A. Jacobson,
  L. Kavan & O. Sorkine-Hornung, "Robust inside-outside segmentation
  using generalized winding numbers", *ACM TOG* 32(4), 2013: the winding
  number of a closed triangulated surface as a sum of solid angles, the
  first classification of §9 (replaced by exact ray parity, below).
* E. Haines, "Point in polygon strategies", *Graphics Gems IV*, 1994, and
  J. O'Rourke, *Computational Geometry in C*, 2nd ed., Cambridge, 1998,
  §7.5: inside/outside by the parity of a ray's crossings, with crossings
  decided by orientation predicates (§9).
* C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann, 2005:
  the closest point of a triangle (§5.1.5), the slab test of a segment
  against a box (§5.3.3) and bounding-volume hierarchies built top down
  by median splits (§6.1–6.2), used to index every search in §9.
* Our design notes 0006 (the tolerant 2D arrangement, reused in §4.2),
  0008 (periodic faces split at π) and 0011 (intersections).

## 3. Inputs, outputs and invariants

`Kernel::boolean(a, b, op) -> Result<Body, KernelError>` with
`op ∈ {Union, Intersect, Subtract}`; feature-level variants pick which
existing parts a new body combines with. Both inputs must validate.
The output validates, or the call fails with an error naming the stage
and the location (a point in model space), and nothing is produced.
Results that would be non-manifold (two blocks sharing exactly an edge
or a vertex) are refused in Phase 3 with that message; non-manifold
bodies are their own later decision (an ADR when it comes up).

Every face, edge and vertex of the result carries a generation record
(design note 0013): which input entity it came from, or which pair of
faces created it, so references downstream can follow it.

## 4. Algorithm

### 4.1 Intersect

1. Broad phase: a bounding-volume hierarchy over the faces of each body
   (boxes inflated by the tolerance); candidate face pairs whose boxes
   overlap.
2. For every candidate pair, design note 0011: intersection curves with
   both pcurves, piercing points of each body's edges with the other's
   faces, edge/edge contacts, and coincident regions.
3. **One decision per fact.** Every vertex/face, vertex/edge and
   edge/face incidence is computed once and stored in an incidence table
   keyed by the entity pair; every later stage reads the table instead of
   re-evaluating a predicate. A piercing point found from the edge side
   and the same point found as the end of an SSI branch are merged within
   the tolerance into one new vertex, by the same union–find clustering
   as the sketch arrangement (design note 0006 §4), so the two faces and
   the edge agree on it by construction.

### 4.2 Split (imprint)

Edges are split at their piercing vertices. Each face is split along the
intersection curves that cross it by building, in the face's parameter
space, the arrangement of its boundary loops and the pcurves of its
intersection curves: the problem design note 0006 solves for sketches,
with the tolerance mapped into parameter space from the surface's
metric (`|S_u|`, `|S_v|` bounds over the face, so a parameter-space
distance under the mapped tolerance is under the 3D tolerance). The
faces of that arrangement inside the original face are the new faces;
reusing it brings its tested handling of tangencies, overlaps,
near-coincident vertices and sliver absorption. Periodic faces are
already split at π (ADR 0008), so parameter spaces are simply
connected rectangles.

### 4.3 Classify

Each new face of A is classified against B as inside, outside, or on
(coincident, with the same or the opposite normal), and likewise each
face of B against A:

* Faces adjacent to an intersection edge are classified from the local
  geometry at that edge: the side of the other body's face, from the
  surface normals at a point of the edge, with the tangential cases of
  Mäntylä's vertex/edge neighbourhood tables.
* The classification propagates to every face connected to a classified
  face without crossing an intersection edge (flood fill), as in Laidlaw
  et al.; faces in components untouched by any intersection are
  classified by one point-in-solid test (a ray cast with robust
  predicates, repeated along a second random direction when the first
  passes within the tolerance of an edge).
* A face whose propagated class conflicts with a direct test is an
  internal inconsistency: an error, never a coin toss.

### 4.4 Select

| op | kept from A | kept from B |
|---|---|---|
| union | outside, on-same | outside |
| intersect | inside, on-same | inside |
| subtract | outside, on-opposite | inside (reversed) |

### 4.5 Merge and validate

Kept faces are sewn along the shared intersection edges (each such edge
appears once in A's split and once in B's; they are the same tolerant
edge by construction of §4.1) and vertices. Then:

* faces on the same surface that now meet along an edge whose two sides
  are the same surface and continue each other are merged (so a cut and
  re-joined plane does not leave a seam), and edges between two pieces of
  one curve meeting at a vertex of valence two are merged; this
  simplification is optional per call (patterns want it, debugging may
  not), on by default;
* the validator runs on the result; a failure is reported with the stage.

## 5. Tolerances

The model tolerance is the size below which points are the same.
Approximate intersection curves carry their measured deviation as a
local edge tolerance (Jackson), at most a tenth of the model tolerance;
vertices where tolerant edges meet carry the largest distance to their
edges' ends. The validator checks each edge against its two surfaces
within its own tolerance. Tolerances never grow across operations
without bound: an operation whose inputs already carry local tolerances
near the model tolerance is refused rather than accumulating (a later
healing step can recompute such edges exactly).

## 6. What is refused

Named with a location, never guessed: missing or inconsistent
intersection branches (design note 0011); classification conflicts;
non-manifold results; faces or edges thinner than the tolerance left by
near-coincident geometry (the same choice the sketch arrangement makes
for faces it cannot absorb as artefacts of the tolerance).

## 7. Tests

* **Volume identities** on every random pair: `|A ∪ B| = |A| + |B| −
  |A ∩ B|` and `|A − B| = |A| − |A ∩ B|`, with volumes from the exact
  mass properties (design note 0003) to 1e-9 relative, plus the same for
  surface area where the pair has no coincident faces.
* Golden cases: blocks overlapping in every combination of shared faces,
  edges and vertices (coplanar, touching, nested, disjoint); cylinders
  through blocks (holes), coaxial and crossing cylinders, a sphere
  through a block, revolves cut by extrudes; each with volume, area,
  topology counts and validity, at scales 1e-3 to 1e3.
* Fuzzing: random rigid motions of golden pairs, including motions by
  fractions and small multiples of the tolerance (where coincidences
  appear and vanish), asserting validity and the identities or a named
  refusal.
* **Failure-rate metric**: refusals per thousand random pairs per class,
  printed by CI on every run, so robustness is measured, not claimed
  (the Phase 3 roadmap item).

## 8. Order of work

1. Plane/plane (polyhedral solids from extruded polygons): the whole
   pipeline end to end with the parameter-space arrangement, the
   incidence table, classification and sewing, and the identity tests.
2. Planes with cylinders, cones and spheres in closed form (holes,
   bosses, revolves of lines and arcs).
3. Extrusions and revolutions of NURBS profiles via the reductions of
   design note 0011 §4.1, then the general intersection path.
4. Features: extrude and revolve with add, remove and intersect; the
   first part-combining UI; patterns and mirror once booleans are
   measured robust.

## 9. As built: step 1, polyhedral solids (2026-10-13)

`Kernel::boolean(a, b, op)` (ADR 0020) handles solids bounded by planes
with straight edges: boxes and extrudes of polygons, holes included.
Bodies with curved faces or edges fail with `KernelError::Unsupported`
until step 2. The code is `crates/kernel/src/boolean/`, one module per
stage.

**Points** (`points.rs`). Each vertex of either body gets one sign
(−1, 0 within the tolerance, +1) per face plane of the other body, once.
The candidate vertices of the result are the bodies' vertices, the
points where an edge crosses a plane of the other body (its ends'
signs strictly opposite) and the points where two non-parallel edges
come within the tolerance; points within the tolerance of each other
are joined into clusters (union–find on a grid), each placed at a
body's vertex if it holds one, else at the mean. A cluster holding two
vertices of one body, or spreading wider than twice the tolerance, is
refused. Whether a cluster lies on an edge is then decided once, by its
distance from the edge, and the edge is split at the clusters on it.

**Imprint** (`mod.rs`). Faces on two planes meet along the stretch of
the planes' common line inside both closed faces, found from the
vertices' signs by two half-open parity rules (counting a vertex on the
plane as on its positive and then its negative side; their union is the
closed face); no new geometric test is made. Faces on one plane (all
vertex signs zero both ways) imprint the pieces of each other's edges
that lie in the other closed face; a face on the other's plane only one
way round is refused as nearly coincident.

**Split** (`planar.rs`). Each face's edge pieces and imprints form a
planar graph on known points, walked with exact orientation predicates
(Shewchuk); counter-clockwise cycles are regions, clockwise ones the
holes of the smallest region of another component around them.

**Classify** (`classify.rs`). This replaces the propagation of §4.3. A
piece inside a coincident face of the other body is on it, with the same
or the opposite normal; any other piece is classified at the incentre of
its triangle with the largest inscribed circle in a constrained Delaunay
triangulation (a convex piece, at the better of its centroid, which is
at least two thirds as far from the boundary as the best point, and the
incentres of a few of its triangles; a face farther than the tolerance
from the other body's box is outside it without a test). Each piece is
decided by its own test, so no decision depends on another, and the
test's margin is measurable: a piece whose
inscribed circle is no wider than the tolerance is a sliver, and a point
within the tolerance of the other body's surface is on it; both are
refused. Any other point is decided exactly by ray parity: the other
body's faces are triangulated on its own vertices, so the triangles
close up edge to edge, and a segment from the point to beyond the body
counts the triangles it passes through, +1 leaving and −1 entering, each
decided by orientation predicates (Shewchuk). The count is the winding
number of the triangulated surface: 1 inside, 0 outside, anything else
an error. A segment through a triangle's edge or corner is not counted;
the test is repeated along the next of eight fixed directions. (The
first version summed solid angles, van Oosterom & Strackee, Jacobson et
al.: right, but approximate and linear in the other body's size for each
piece.)

**Indexes.** Every search is indexed (`crates/kernel/src/bvh.rs`, a
bounding-volume hierarchy built by median splits, deterministic): the
face pairs whose boxes meet, the edge pairs that may touch, the clusters
near an edge, the ray and distance tests of classification, and, for
faces with more than 24 segments, the split's crossing and closeness
checks and the hole-to-region search. The index only prunes; every
candidate is decided by the same exact test as before, so the results
are unchanged (the fuzzing below gave identical counts before and after).
Triangulations are bulk loaded (with the triangulation library patched so
the same input triangulates the same way on every run, ADR 0022), joins of
collinear edges are made incrementally, and the stage's work is capped
(docs/PERFORMANCE.md):
the union of two 1,024-sided prisms (about 4,000 result faces) went from
5.4 s to 0.6 s.

**Select, sew, merge** (`assemble.rs`). Pieces are kept by the table of
§4.4. Every segment must then be used exactly once in each direction (a
segment used more often is a non-manifold result, refused; fewer is a
kernel error), adjacent coplanar faces with the same orientation are
merged by re-running the split on their joint boundary, collinear edges
meeting at a point no third edge reaches are joined, and a vertex whose
faces form more than one fan is refused as non-manifold. Shells are the
connected components; outer shells have positive volume, and each void
goes to the solid whose outer shell winds around it. Origins follow ADR
0020.

**What is refused**, each with a location: an empty result; a
non-manifold result (bodies touching along an edge or at a point);
points of the bodies closer than the tolerance (two vertices of one body
in a cluster, a cluster spread wider than twice the tolerance, two
points on one edge or one line within the tolerance of each other, a
point within the tolerance of two edges of one body away from their
vertex, or of an edge it is not on); edges passing within the tolerance
of each other without meeting; faces that nearly coincide; slivers
thinner than the tolerance; a face running within the tolerance of the
other body's surface without lying on one of its faces. Internal
inconsistencies (a piece without a partner across an edge, a region
both inside and outside a face) are `InvalidResult` errors, which the
tests count as failures.

**Measured** (`crates/kernel/tests/booleans.rs`). Golden cases: two
overlapping cubes (volume, area, topology), stacked and flush boxes that
merge into one, a box against part of a side (a T-shaped prism), a slot
and a through hole (genus 1), nested bodies (a void), disjoint and
identical bodies, bodies touching along an edge or at a corner, the
cubes and the hole at 1 mm and 1 km, and the origin tables of a union.
Random pairs of boxes and extruded star-shaped polygons (some with
holes; half on a quarter-unit grid so that faces, edges and vertices
coincide; about two thirds on axis-aligned frames):
all four operations, checked by `|A ∪ B| + |A ∩ B| = |A| + |B|`,
`|A − B| + |A ∩ B| = |A|` and `|B − A| + |A ∩ B| = |B|` with exact
volumes to 1e-9 relative. Every result in every test is validated and
its origins checked against the operands' geometry: each face on its
operand face's plane and facing the right way, each edge along its
operand edge (a piece of it, or reaching it when joined from several)
or on both planes where the bodies meet, each vertex at its operand
vertex, all within three tolerances. On 6,000 pairs (20 times the default run)
every result is valid, every identity holds where the three results
exist, and 11 of 24,000 operations are refused (0.5 per thousand), each
because the result would be non-manifold (bodies on the grid touching
along an edge or at a point), which is the right answer. The tolerance
test moves a copy of a random body by 0.25 to 10 tolerances and checks
the same identities with a slack of the area times the motion; it is the
adversarial case (every face against its own copy). On 1,400 pairs no
result is invalid or inconsistent. Refusals per shift (of 800
operations each): 4 at a quarter of the tolerance and 24 at a half (for
points closer than the tolerance: these are coincidences the stage
should absorb, the first robustness target), 108 at 0.9, 564 at 1.1,
all 800 at 1.5, 524 at 3 and 156 at 10 tolerances, where faces nearly
parallel to the motion move by about a tolerance and edges nearly
parallel to it leave slivers thinner than the tolerance. CI prints
these rates on every run. The surface-area identity waits for step 2's
tests.

