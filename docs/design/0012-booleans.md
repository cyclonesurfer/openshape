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

**Points** (`points.rs`). Faces of one body built on one plane share
it: their normals parallel and their offsets equal up to rounding
(relative 1e-10, the sine of the angle between the normals; a test on
the cosine, `|n · m| ≥ 1 − 1e-10`, passed angles up to 1.4e-5 rad and
took two sides of a prism meeting at a slight crease for one plane,
leaving the far corners of the second microns off it). Each vertex of
either body gets one sign (−1, 0 within the tolerance, +1) per plane of
the other body, once.
The candidate vertices of the result are the bodies' vertices, the
points where an edge crosses a plane of the other body (its ends'
signs strictly opposite) and the points where two non-parallel edges
come within the tolerance; points within the tolerance of each other
are joined into clusters (union–find on a grid), each placed at a
body's vertex if it holds one, else at the mean. A cluster holding two
vertices of one body, or spreading wider than twice the tolerance, is
refused. Two clusters of derived points (crossings and contacts, no
vertex of either body) on one edge within the tolerance of each other
along it are then joined, and the joining repeated until none are left:
where an edge meets a plane at a shallow angle θ, the tolerance across
the plane is a stretch of ε / sin θ along the edge, so the crossing found
from each side (one body's edge through the other's plane, and the
other's edge through the first's) can land just over the tolerance
apart in space while within it along the edge. Vertices are the bodies'
own data and are never joined this way: joining copies of vertices moved
by just over the tolerance made an inconsistent result (caught by the
validator) where a refusal is right. Whether a cluster lies on an edge
is then decided once, by its distance from the edge, and the edge is
split at the clusters on it.

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
or the opposite normal. Every edge of that face inside the piece's face
bounds pieces, so a piece is wholly on it or wholly off it, and the
piece's test point, more than the tolerance inside it, is on the piece's
side of each of the face's edges: parity alone decides. (A band of the
tolerance about the face's edges took a strip 2.4 tolerances wide
between the holes of a body and of its copy moved by three tolerances,
off the face but within the tolerance of a corner of it, for a piece on
it, and made an inconsistent result; the strip is now classified
against the other body, which refuses it as within the tolerance of its
surface.) Any other piece is classified at the incentre of
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
kernel error), adjacent faces on one plane facing the same way are
merged by re-running the split on their joint boundary (faces of one
body on one of its planes, or of the two bodies found to coincide in
the imprint: the decision already made, not a new test on the normals), collinear edges
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
operations each): none at a quarter or a half of the tolerance (4 and
24 before crossings close along an edge were joined, above; the test
now requires none), 8 at 0.9, 548 at 1.1, all 800 at 1.5, 524 at 3 and
156 at 10 tolerances, where faces nearly parallel to the motion move by
about a tolerance and edges nearly parallel to it leave slivers thinner
than the tolerance. CI prints these rates on every run.

A run at a hundred times the default size (30,000 random pairs, 7,000
motions) found two faults the smaller runs had not, both now fixed and
kept as tests. Two sides of a prism at a slight crease were taken for
one plane (Points, above) and merged into one face with a corner 4 µm
off its plane; the validator let it through, because it held planes to
the hundred tolerances it allows curved surfaces for their iterative
projections, and planes are now held to ten (a result's vertices are
at most 1.06 tolerances off their faces' planes over the fuzz at twenty
times its default size).
And a strip between moved holes was classified on a face it only
touched (Classify, above), an inconsistent result the kernel caught.
At that size every result is now valid, 47 of 120,000 operations on
random pairs are refused, each a non-manifold result, and the motions'
refusals per shift are 0, 0, 20, 2,828, 4,000, 2,476 and 680 of
4,000. The surface-area identity waits for step 2's tests.


## 10. Step 2: cylinders, cones and spheres (design, 2026-10-13)

Step 2 adds faces on circular cylinders, cones and spheres (what extrudes
and revolves of lines and arcs make) and edges on circles and conics, with
closed forms wherever two such surfaces meet in curves of degree two or
less. It keeps step 1's architecture rather than the per-face
parameter-space arrangements of §4.2: every fact about the two bodies'
contact is decided once, in space, and every later stage reads it
(Hoffmann, §2). A per-face arrangement would decide again, in each face's
own coordinates, where curves cross, and two faces could disagree about
the same point; step 1's fuzzing found exactly such disagreements between
a face's plane and space (§9), which its 3D-first stages now rule out.

### 10.1 Surfaces and curves

* **Surfaces.** Each operand face's surface is recognised as a plane
  (point, unit normal), a cylinder (a point of its axis, unit axis,
  radius), a cone (apex, unit axis, half-angle) or a sphere (centre,
  radius), exactly from the kernel's surface records (an extrusion of a
  circle is a cylinder; a revolution of a line or of an arc centred on
  the axis is a cylinder, cone, plane or sphere), else the boolean is
  `Unsupported` (tori and NURBS: step 3). Faces of one body on one
  surface share it, recognised up to rounding as planes are in step 1
  (the sine of angles, offsets relative to 1e-10): the two halves of a
  cylinder (ADR 0008) are one surface. Each surface has an implicit
  form G (degree one for a plane, two for the others: for a cylinder
  `|p − a|² − ((p − a)·u)² − r²`, for a sphere `|p − c|² − r²`, for a
  cone `((p − a)·u)² sin²α − (|p − a|² − ((p − a)·u)²) cos²α` on its
  nappe) and a signed distance, positive outside (away from the axis or
  centre; the side the plane's normal points to).
* **Curves.** Edges and the curves where faces meet are lines, circles,
  and conics (ellipses from a plane across a cylinder or a cone,
  parabolas and hyperbolas from a plane across a cone), the conics as
  rational quadratic B-splines, exact (Piegl & Tiller, *The NURBS Book*,
  2nd ed., ch. 7). For the root finding below every curve is split into
  rational Bézier pieces of at most a quarter turn (a circle exactly,
  with weights cos(θ/2); Piegl & Tiller §7.5).

### 10.2 Points

The candidate vertices are found as in step 1, with the decisions
generalised from planes to surfaces:

* **Signs.** Each vertex gets one sign (−1, 0 within the tolerance, +1)
  per surface of the other body, from its signed distance, once.
* **Crossings of an edge with a surface.** On a straight edge and a
  plane the signed distance is monotone, so the ends' signs decide, as
  in step 1. Otherwise it is not (a line passes into and out of a
  cylinder; an arc crosses a plane twice with both ends on one side).
  Composing G with a rational Bézier piece gives a polynomial (degree at
  most four over the square of the weight), and the roots of its
  derivative split the edge into pieces on which the distance is
  monotone. Each piece end that is not a vertex (a turning point of the
  distance) gets a sign from its distance, decided once like a vertex's;
  a crossing is where a piece's end signs are strictly opposite, found
  by bisection of the monotone piece, and recorded by (edge, surface,
  piece). A turning point within the tolerance of the surface is a
  tangent contact: the edge touches the surface there. It is a
  candidate point too, and step 2 refuses tangent contacts between an
  edge and the interior of a face (a cylinder resting on a plane) with
  that message until tangency is handled (design note 0011 §4.4).
* **Roots** of these polynomials are isolated in Bernstein form by
  subdivision, the number of sign changes of the coefficients bounding
  the number of roots in an interval (the variation-diminishing
  property; J. M. Lane & R. F. Riesenfeld, "Bounds on a polynomial",
  *BIT* 21, 1981; Farin, ch. 5), to a fixed depth; an interval still
  undecided at that depth holds a near-double root, a near-tangency,
  which the distance at its middle decides like any turning point.
* **Contacts of two edges** (crossings within a face both bodies lie
  on, such as two coplanar arcs) in closed form for lines and circles,
  and by subdivision of rational Bézier pieces with the distance
  minimised on the exact curves for conics (design note 0011 §6).
* Clustering, the joining of derived points along an edge, and the
  points' incidences with edges are step 1's, with distances to circles
  and conics exact (closed form for circles, Newton on the exact curve
  for conics, from the Bézier piece's nearest chord point).

### 10.3 Imprint

Two faces on crossing surfaces meet along the surfaces' intersection,
which within step 2 is one of: a line (two planes), up to two lines (a
plane along a cylinder's axis, or along a cone's apex; two cylinders with
parallel axes), a circle (a plane across the axis of a cylinder, cone or
sphere; two coaxial surfaces of revolution; two spheres), or a conic (a
plane obliquely across a cylinder or cone). Other pairs (two cylinders
at an angle: a quartic) are `Unsupported` until step 3's general path.
Surfaces that coincide (one plane, one cylinder: equal axis and radius
within the tolerance at every vertex) are coincident faces as in step 1.
Two surfaces tangent along a curve or at a point (a cylinder on a plane,
two cylinders touching side by side) are refused as tangent contacts.

The stretch of an intersection curve inside a face is found exactly as
step 1 finds a line's: the curve enters and leaves the face only where
the face's boundary crosses the other surface, so the cuts along it are
the boundary's crossings with that surface (§10.2) and its vertices on
it, paired by the two half-open parity rules; a closed curve (a circle)
that no boundary crosses is wholly inside or outside the face, decided
at one of its points. Positions along a curve are its length from a
cut (a circle's angle times its radius), so the comparisons within the
tolerance are in space, as in step 1.

### 10.4 Split

A face's segments now include arcs. Segments are identified by an index
rather than by their two clusters (two arcs, or an arc and a line, can
join the same two points); two segments between the same clusters
whose midpoints are within the tolerance are the same segment, with both
sources. The face's planar graph is walked in a chart of its surface in
which lengths are lengths in space (a plane's frame, a cylinder or a
cone unrolled; the sphere, which has no such chart, in its angles with
the tolerance scaled by the metric's bounds over the face), so the
tolerance means the same everywhere. At each point the segments are
ordered by their exact tangent directions in the chart, then by their
curvature where tangents agree; two segments leaving a point in one
direction with one curvature are refused. Areas (orientation of
cycles), containment of holes and the checks that segments neither cross
nor pass within the tolerance of a point they do not end use the exact
curves where closed forms exist (lines and arcs in a plane, rulings and
circles on a cylinder or a cone, which unroll to lines) and otherwise
polylines within a sixteenth of the tolerance of the curve, so a
decision made on a polyline can be off only within that band (a refusal
for points within 1.06 rather than 1 tolerances, never a wrong region).

### 10.5 Classify

A piece's test point is the incentre of the largest triangle of a
triangulation of its boundary in the chart (polylines for curved
segments), with its clearance measured in the chart (a lower bound on
the distance in space up to the polyline band). Ray parity stays exact
in the sense that matters: the ray meets each face of the other body at
the roots of a polynomial of degree at most two (closed forms), each hit
is counted with the sign of the face's outward normal along the ray if
it lies inside the face (decided in the face's chart against its
boundary, with the hit farther than the tolerance from every boundary
curve, measured exactly), and a ray that passes within the tolerance of
an edge, or grazes a surface (a double root), is abandoned for the next
direction, as one through a triangle's edge is in step 1. Planar faces
with straight edges keep step 1's exact triangles. Pieces on coincident
faces are decided by parity in the coincident face's chart (§9).

### 10.6 Assemble

Pieces are kept by the table of §4.4 and sewn along segment indices;
edges carry their curves (a line, an arc of the circle, a piece of the
conic), and faces their operand's surface. Pieces of one face are merged
again; pieces of coincident faces of the two operands are merged on a
plane, but on a curved surface only when the merged face spans no more
than π about the axis (ADR 0008); otherwise the edge between them stays,
and displays hide it as a seam (both faces are pieces of one surface).
Collinear segments are joined as in step 1, and arcs of one circle or
pieces of one conic likewise. Shells' volumes, for orienting shells and
nesting voids, come from the kernel's exact mass properties.

### 10.7 Order and tests

1. The geometry layer, unit- and property-tested against sampling: the
   surfaces and their recognition, Bernstein root isolation, crossings of
   lines, arcs and conics with planes, cylinders, cones and spheres, and
   the closed-form intersection curves.
2. The pipeline for planes and cylinders (holes, bosses, slots, rounded
   blocks, counterbores, at any orientation), with golden cases (a plate
   with a round hole, a boss, a D-cut shaft, two parallel cylinders, a
   counterbore) checked for volume, area, topology and origins, and the
   fuzz of §7 on random extrudes of profiles with arcs and their motions.
3. Cones and spheres (revolves), with the same tests.

Volumes in the identities come from the exact mass properties, so the
identities hold to 1e-9 for curved results as for flat ones.

## 11. As built: step 2, planes and cylinders (2026-10-13)

§10.7 item 2. Solids bounded by planes and circular cylinders (extrudes
of profiles with lines, arcs and circles, holes included, at any
orientation) combine; cones and spheres (item 3), cylinders at an angle
to each other (a quartic), tangent contacts and free-form edges on a
cylinder (an ellipse a previous oblique cut left, stored as a NURBS
curve) are refused as `Unsupported` with the reason. Polyhedra take
step 1's path unchanged: the 100-times fuzz of §9 gives the same
results, operation for operation, before and after.

**Contacts** (`imprint.rs`, `meet.rs`). How two surfaces meet is found
once per pair of surfaces, not per pair of faces, over the region where
faces of both lie (the intersection of each surface's faces' boxes,
padded by the tolerance): apart, the same surface, touching (refused as
a tangent contact), unsupported, or tracks, the curves of §10.3 (lines,
circles, ellipses), each with an id that every later stage uses. A face
whose boundary lies strictly on one side of the other surface is apart
from it when the face lies in the convex hull of its boundary (a face on
a plane or a cylinder) and that side is one side of a plane or the
inside of a cylinder (convex); and when the face is on a cylinder and
the surfaces meet in curves found in closed form, which are lines along
its axis and circles about it, so they would cross its boundary to reach
it. (As first committed, any curved face with its boundary on one side
was taken to be apart, which misses a closed quartic curve inside a
face: a rod dipping into a pipe's wall between its seams, where no edge
of either crosses the other, was returned as two bodies that do not
meet. Found the same day, before cones and spheres, and fixed: such a
pair is refused as a quartic; `a_rod_dipping_into_a_pipe_wall_is_refused`
and the membership check below keep it so.)

**Points** (`points.rs`). Along each edge and surface of the other body,
the *stations* are the edge's ends (with their vertices' signs) and the
turning points of the surface's implicit function along the edge
(§10.2), each with a sign decided once. A crossing lies between two
stations of strictly opposite signs and is found by bisection; a
station within the tolerance of the surface, between stations off it,
is a point too (an arc touching a plane). Crossings are keyed by edge,
surface and the stations around them, so the two faces of an edge find
the same point.

**Imprint.** A track's cuts in a face are its boundary's crossings with
the other surface and its stations on it, by step 1's two half-open
rules. Each cut is assigned to the nearest track within 16 tolerances
(a point within the tolerance of both surfaces lies within the
tolerance over the sine of their angle of the curve where they meet);
a cut as near to a second track is refused. On a closed track a cut
enters the face where the boundary goes to the positive side of the
other surface and the track runs across the face one way
(`T′ · (n_F × ∇G)` > 0, `T` the track's tangent, `n_F` the face's
outward normal, `G` the other surface's implicit function), and leaves
it otherwise; the cuts at one point are taken together (their net is
±1, or 0 for the boundary touching the track there), and the enters and
leaves must alternate around the track. A closed track that no cut
reaches is wholly in the face or wholly out of it, decided by the face's
outline at four of its points. The stretch in both faces is the
intersection of the two faces' spans, along a line as in step 1 and
around a closed track as arcs. A track touching a face's boundary from
inside it is refused (a tangent contact).

**Segments.** A segment's course is straight, a piece of an operand's
edge, or a piece of a track; curved segments between the same two
clusters are one segment when their midpoints are within the tolerance.
Pieces of circles are reparametrised so that the arc's middle is at a
quarter turn from the parameter's origin, so no arc crosses it.

**Split** (`chart.rs`, `planar.rs`). A face is split in its chart (a
plane's frame, a cylinder unrolled about the face's middle, both
isometric), each curved segment a chain of chords with bands from a
bound on its course's second derivative (§10.4). Chords are halved
until no band can hide a contact; the chords at a node shared by two
segments are halved at once (they always meet there), and each later
round tests only its new chords. A chord's band is its course's sag
about it plus how far the chord's ends were moved to their clusters
(which take in points within the tolerance); halving narrows only the
sag, so a chord is halved while its sag, not its band, is over a
sixteenth of the tolerance. Ordering segments at a point by exact
tangents and curvature (§10.4) was not needed: the chords near a node
are within a sixteenth of the tolerance of their courses, besides the
move of their ends to the node's cluster.

**Classify** (`classify.rs`). A piece's point is found as in step 1,
with the chords' bands subtracted from its clearance and the clearance
measured in space; chords are halved while the clearance left is within
the tolerance. A piece with more than 64 boundary points is first probed:
from the middle of a few sides inward to the first side met, halfway,
each candidate checked by its winding number (exact orientation tests;
Hormann and Agathos, "The point in polygon problem for arbitrary
polygons", *Computational Geometry* 20, 2001) and measured against every
side; the triangulation serves when no probe is far enough inside.
Parity rays meet curved faces in closed form (§10.5); a ray grazing a
surface (the cosine of its angle with the surface's normal under 1e-3)
or meeting a face within the tolerance of its boundary is abandoned for
the next direction.

**Assemble** (`assemble.rs`). Pieces on a curved surface merge only with
pieces of the same operand face facing the same way (they span no more
than the face, half a turn); flat pieces merge as in step 1. A flat
face may be recorded as a surface of revolution (a revolve's cap, in
faces of half a turn each, ADR 0008): joined from pieces, it can go
around the axis, which that record's parameters cannot hold, so a joined
face on a plane is recorded as a plane (found 2026-10-13: before, a
turned rod cut square by a box came back with its cap's parameters
wrapped and its volume wrong, and validation now refuses any face of
revolution whose boundary leaves its frame's half turn). A point on the
axis of a face's surface of revolution stays a vertex: segments are not
joined there, since an edge running through a singular point of its
face's surface has no parameters on the face there (its angle jumps by
half a turn), and validation refuses such an edge. Curved
segments join end to end only when they are pieces of one operand edge
or one track; straight ones join as in step 1. A result whose faces are all flat with straight
edges gets step 1's exact shell volumes. Otherwise a lone shell is the
outer shell of the one solid (validation refuses it if it is inside
out); several shells are ordered by volumes from validation's coarse
mesh with an error bound (twice the area of its curved faces' mesh times
its distance from the surfaces; a flat face's mesh lies in its plane),
or by the exact mass properties where that bound could reach the sign,
and each void goes in the smallest solid holding one of its own corners.
Validation decides orientation by the same rule (§13).

**Tests.** Golden cases (`crates/kernel/tests/curved.rs`) check volume
and area to 1e-9 relative, topology and the origins of every face, edge
and vertex: a plate with a round hole, a boss through and on a plate, a
D-cut shaft, two parallel cylinders (all three operations), a
counterbore, a box cut from inside a cylinder and a cylinder from
inside a box (voids), and a slab cut across a cylinder (two solids).
Random pairs check the three volume identities of §7 with exact
volumes: boxes at any angle and prisms of circles, slots, rounded
rectangles, D shapes, polygons with round holes and discs with
polygonal holes, along one axis (the world's or a turned one), half of
them on a quarter-unit grid. On 10,000 pairs (40,000 operations) every
result is valid with its origins on its operands' surfaces and every
identity holds where the results exist; 525 operations are refused
(13 in 1,000): 516 as tangent contacts (on the grid, a box's face on a
cylinder's side is common), 5 as non-manifold results, and the 4 of one
pair whose box rests a corner exactly on the centre of a cylinder's end
face, from inside: the face's test point is its centre, within the
tolerance of the box, and the face is refused although the union
exists. Keeping test points clear of such contact points would decide
it. The run found one fault, now fixed and kept as a golden case: a
void inside a curved body was placed by a corner of the outer shell
instead of its own, so cutting a body from inside another was
refused (88 operations of the first run).

Each pair is also checked by point membership: 24 points sampled where
the operands' boxes meet, each decided on meshes within 1e-3 of the
operands and of every result (by the winding number, with points within
2e-3 of any mesh skipped), must be inside a result exactly when the
operation's rule on the operands says so. This shares no code with the
booleans and catches what the identities cannot: the missed contact of
the corrected rule above gave a union, a difference and an intersection
that satisfy all three identities. A second run takes prisms along two
different axes: on 2,000 pairs every result is correct and every
refusal (383 in 1,000 operations) is a quartic of two cylinders at an
angle.

**Turned solids** (2026-10-13). With rods, tubes and stepped rods
turned a whole turn or half of one about the pair's axis among the
random solids (a fifth of them), a run at 20 times the default size
found three faults, now fixed and kept as tests: an edge through the
centre of a cap (Assemble, above), whose mass properties did not
converge (102 results) or, for 7 pairs, converged to volumes that broke
the identities by up to 3e-7 of them; end chords halved without end
(Split, above; 8 operations); and a cap's centre that came out of its
chart at 8e-49, which the triangulation of a piece refuses as too small
(4 operations; coordinates under 1e-30 of a polygon's size are now
taken as zero). After the fixes, on 4,000 pairs along one axis and
2,000 at an angle, every result is valid and correct by the identities
and by point membership; along one axis 317 of 16,000 operations are
refused (20 in 1,000: 312 tangent contacts, 5 non-manifold results),
and at an angle 3,730 of 8,000, all quartics but one non-manifold
result and one elliptical edge on a cylinder (§10.7, later items).

## 12. Step 2, item 3: cones and spheres (design, 2026-10-13)

Item 3 adds faces on cones and spheres (revolves of lines and arcs) to
§11's pipeline. Their surfaces are already recognised, and the
geometry layer already finds crossings of edges with them and their
circles with coaxial surfaces (§10.7, item 1). What item 3 adds is a
chart of each such face, the curves where a plane crosses a cone off
its axis, the closed-form ray hits of classification, and care at the
apex and the poles, where a face's surface has no normal.

### 12.1 Charts

* **Cone.** A cone is developable: cut along a ruling and unrolled, it
  is a sector of the plane, an isometry (M. do Carmo, *Differential
  Geometry of Curves and Surfaces*, 1976, §4-2). A point at slant
  distance `s` from the apex and angle θ about the axis goes to polar
  coordinates `(s, θ sin α)`, α the half-angle: rulings go to lines
  through the apex's image, circles about the axis to circles about it.
  As for a cylinder the cut is opposite the face's middle (faces span at
  most half a turn, ADR 0008). The chart is anchored at the face: in
  coordinates relative to the point of the middle ruling at the cone's
  reference height, `(σ cos ψ − 2 s₀ sin²(ψ/2), (s₀ + σ) sin ψ)` with
  `ψ = θ sin α`, `σ = s − s₀`, `s₀` the slant distance of the anchor, so
  a slight taper (its apex kilometres away) costs no precision and the
  chart tends to the cylinder's as α goes to 0. Outer loops run
  counter-clockwise about the outward normal, as for the other charts.
  Lengths in the chart are lengths along the surface; a chord in space
  is shorter, by at most the factor `sin α / sin(π sin α / 2) ≥ 2/π`
  between two points of a face (their angles about the axis differ by
  at most π).
* **Sphere.** A sphere has no isometric chart (its curvature is not
  zero; Gauss's *theorema egregium*). Its chart is the stereographic
  projection from the antipode of a centre `c` of the face onto the
  tangent plane at `c`: a point at angle φ from `c` goes to distance
  `2R tan(φ/2)` from `c`'s image. It is conformal, and takes circles of
  the sphere to circles or lines of the chart (T. Needham, *Visual
  Complex Analysis*, 1997, ch. 3; D. Hilbert and S. Cohn-Vossen,
  *Geometry and the Imagination*, §36). Its scale at chart point `q` is
  `λ(q) = 1 + |q|²/(4R²) ≥ 1`, and the distance in space between the
  points at `q₁` and `q₂` is exactly `|q₁ − q₂| / √(λ(q₁) λ(q₂))`
  (the chordal distance; L. Ahlfors, *Complex Analysis*, 3rd ed., ch. 1,
  §2.4). The centre is the middle of the half of the sphere that holds
  the face: every face of a sphere lies within longitudes `[0, π]` of
  its surface's frame (ADR 0008), so within a quarter turn of the
  direction of longitude π/2 on its equator, where the scale is at most
  2; the antipode is a quarter turn from the face, and the projection is
  one to one on it. (The middle of the face's own boundary would do as
  well but must be found from samples, and a lune's boundary, two
  meridians, has no middle by averaging.) Operands are checked: the
  least cosine of the angle from the centre along each boundary edge is
  found in closed form, and a face reaching more than 120° from it
  (scale 4) is refused; within that the bounds below still hold, since
  they are local. The poles are ordinary points of the chart.
* **Lengths.** A chart compares lengths with space in two ways, each a
  bound: `reach(q, d)`, a chart distance from `q` beyond which every
  point is more than `d` from it in space, and `space(q, d)`, a distance
  in space within which no point `d` or more from `q` in the chart lies.
  A plane's are `d` and `d`. A cylinder's are the arc and the chord of
  `d` about its radius. A cone's are `d` over and times the factor above,
  with `reach` taken nearer to `d` away from the apex, where two points
  near in the chart differ little in angle. A sphere's follow from the
  chordal distance at `q`. Every test of the split and of point location
  that was "within the tolerance in the chart" becomes "within `reach(q,
  eps)`", so the tolerance means distance in space on every surface;
  bands (a course's distance from its chords) are measured in the chart
  and bound the distance in space, since no chart lengthens a chord.
* **Courses and bands.** On a cone (isometric) a course's band is the
  bound of §11 from its acceleration: for a circle about the axis, of
  radius `r`, exactly `r sin α` per unit angle squared; for an ellipse its
  greater semi-axis; for a parabola `2a`; for a hyperbola branch
  `√((a cosh T)² + (b sinh T)²)`, `T` the larger parameter of the stretch.
  On a sphere every course is a circle, whose image is a circle of
  curvature `|n · (C − S)| / (2 R r)` (`n` the circle's normal, `C` its
  centre, `r` its radius, `S` the projection's point): the band of a
  chord of length `L` is the sagitta `(L/2)² κ / (1 + √(1 − (κL/2)²))`
  of an arc turning less than a half turn, which the first division
  ensures from a bound on the turning (the curvature times the arc's
  length in the chart, at most `λ` times its length in space, with `λ`
  bounded over the stretch in closed form).

### 12.2 Curves

A plane across a cone meets its nappe in a conic. In the plane's frame
`(x, y)` with `x` along the axis's direction in the plane, `y` square to
the axis, and the origin at the foot of the cone's reference point, the
cone's implicit function is `A x² + y² − 2 B x + C` with
`A = 1 − g² / cos² α`, `g` the sine of the angle between the plane's
normal and the axis: an ellipse where `A > 0` (the plane steeper than
the rulings), a parabola where `A = 0`, one branch of a hyperbola where
`A < 0` (Apollonius, *Conics*, book I). Each is a track: the ellipse
closed, the parabola `vertex + a t² x + 2 a t y` and the branch
`centre + a cosh t x + b sinh t y` open and curved, imprinted along
their courses (not as chords) and stored as exact rational quadratic
B-splines (Piegl & Tiller §7.5). A configuration within a sixteenth of
the tolerance of a parabola over the region is taken as one, so no
ellipse or hyperbola is computed with a far centre. A plane through the
apex meets the cone in up to two rulings, half-lines from the apex: that
case is refused at first ("a plane through a cone's apex") and then
handled with the apex as an end of both half-lines, which a face holding
the apex enters at it. Every other pair of surfaces is as in §10.3:
circles about a common axis, and quartics refused.

### 12.3 Imprint, classify, operands

* **Faces apart.** The rule of §11 as corrected: a face on a plane, a
  cylinder or a cone lies in the convex hull of its boundary, so a
  boundary on one side of a plane, or inside a cylinder, cone or sphere
  (convex insides), keeps the face there; a face on a cylinder or a cone
  holds no closed curve of a meeting in closed form with a curved
  surface. A face on a sphere bulges out of its boundary's hull and can
  hold a whole circle of a meeting (a ball's face around a hole), which
  §11's test for a closed track that no cut reaches decides.
* **Singular points.** At the apex and the poles a surface has no
  normal. A track never passes through the apex: coaxial surfaces that
  meet on the axis touch there, which is refused, and a plane through
  the apex is refused until 3b. A face's facing is found where its
  parameterisation is most regular, at the corner or edge middle where
  `S_u × S_v` is largest, not at its first corner: at a pole of a
  revolved sphere the derivative along the parallel is rounding (about
  1e-17 of the radius) whose cross product points any way, and the
  surface's own singularity test, relative to the derivatives' sizes,
  does not see it (found on the first two-ball union: one face of a
  ball taken as facing inward). A vertex on the axis of a face's
  surface of revolution stays a vertex (§11, Assemble).
* **A whole circle in two faces.** A sphere's face can hold a whole
  circle of a meeting that no edge of either face reaches (a ball's face
  around a hole in a plate's face). The track is split at two new
  points, at its parameters 0 and π, into two halves with ends: no edge
  comes within the tolerance of the track (else a cut would reach it),
  so the new points are clear of every other, and the split of each face
  checks it. A track whole in both faces under one rule of zero signs
  and not under the other has an edge touching it, and is refused.
* **Classify.** A ray meets a cone's nappe and a sphere at the roots of
  a quadratic, anchored at the cone's reference point (a root on the
  other nappe is not a hit); a hit where the ray is within 1e-3 of
  tangent to the surface abandons the ray, as on a cylinder. Near the
  apex the cone's quadratic is flat and its roots lose precision: a
  root that rounding of the coefficients could move by more than a
  sixteenth of the tolerance (their errors over the quadratic's slope
  `2 √disc` there) abandons the ray too. The point of a cone or a sphere
  nearest a point, for the test that a piece's point is clear of the
  other body, is in closed form (on a cone, the foot on the ruling in
  the point's half-plane, clamped at the apex).
* **Boxes.** A face on a sphere is not in the hull of its boundary: its
  box also takes the sphere's extreme points along the axes that lie in
  the face's half of the sphere.

### 12.4 Order and tests

1. Charts, lengths and bands, unit-tested against sampling (map and
   unmap, orientation, `reach` and `space` against distances in space
   over random pairs of points of random faces, bands against dense
   samples of courses); ray roots and feet for cones and spheres.
2. Spheres and coaxial cones end to end: golden cases (a ball on a
   plate, a ball cut by a box, a spherical pocket, two balls, a
   countersink, a cone on a cylinder, a pointed cone cut square to its
   axis) checked as in §11, and a fuzz of revolves about one axis with
   boxes and prisms.
3. Planes across cones off their axis (the three conics), then planes
   through the apex, with golden cases for each conic and the fuzz with
   revolves and boxes at any angle.

Every fuzz checks, besides validity and the volume identities, that
points sampled where the operands' boxes meet are inside each result
exactly when the operation's rule says so, decided on fine meshes away
from every surface: the volume identities hold for a result that takes
two bodies as apart when they meet, and this check does not.

## 13. As built: item 3a, spheres and coaxial cones (2026-10-13)

Faces on spheres and cones now take part in booleans wherever every
curve where two faces meet is a circle: a plane with a sphere at any
angle, two spheres anywhere (their circle is about the line through
their centres), a plane square to a cone's axis, and coaxial pairs of
cones, cylinders and spheres centred on the axis. The charts, courses
and bands are those of §12.1 (a cone unrolled, a sphere projected
stereographically from the antipode of the middle of its half, with a
face reaching more than 120° from that middle refused), the imprint,
classification and operand rules those of §12.3. Everything else is
refused with its reason: a plane across a cone off its axis (a conic,
"not supported yet", item 3b; a test checks a hyperbola, an ellipse and
the two lines through the apex), a sphere or a cone with a cylinder or a
cone off its axis (quartics, as for two cylinders at an angle), and
tangent contacts.

Found while building it, each now handled and tested:

* **Pointed cones about a turned axis** were not recognised as cones:
  the radius at the cone's reference point, on the axis at the apex,
  came out of rounding slightly negative. A radius below zero by less
  than 1e-12 times one plus the profile point's distance from the
  frame's origin is taken as zero.
* **Tangent planes at the poles.** A sphere's normal at a pole is along
  the axis, taken from the side of the pole where the surface is
  regular; a cone has none at its apex, and a face's or an edge's
  origin check in the tests allows two points without one (its corners).
* **Thin shells.** Validation decided a shell's orientation from the
  sign of a coarse mesh's volume alone. The ring a ball keeps outside a
  cone whose side crosses its sphere twice close together (1 cm thick on
  a 1 m ball) meshes there to a volume of 1e-17, whose sign is
  rounding: 5 of 40,000 operations in two fuzz runs at 50 times
  their size were refused as inside out. Validation now decides as
  assembly does (§11): by the coarse volume where its bound cannot reach
  the sign, else by the shell's mass properties. The bound counts only
  curved faces' area: counting flat faces too made it exceed the volume
  of thin flat bodies, and validation then measured them exactly (ten
  times slower; docs/PERFORMANCE.md, fourth pass).

**Tests.** Golden cases (`crates/kernel/tests/curved.rs`), each with all
four operations, volumes exact to 1e-9, topology, validity and the
origins of every face, edge and vertex: a ball resting in a plate, a
ball cut by a box off its centre, a spherical pocket, two balls side by
side and two along one axis, a countersink, a cone on a rod, a pointed
cone cut square, and a ball crossed twice by a cone (the thin ring
above, its volume from the integral of the difference of the two
radii squared). Unit tests classify points near random balls and cones
in random frames against their exact signed distances, and measure the
thin ring's shell exactly. A fuzz of solids about one axis, the world's
or a turned one (balls; cones point up and down; frustums; rods with a
round end; tubes; stepped rods), pairs each with another such solid, a
ball with a box turned any way or another ball anywhere, or a slab
square to the axis, half of them on a grid. On 5,000 pairs (20,000
operations) every result is valid with its origins on its operands'
surfaces, every volume identity holds where the results exist, and
point membership agrees; 1,122 operations are refused (56 in 1,000):
1,096 tangent contacts (on the grid, equal balls touching and slabs at
a pole are common) and 26 non-manifold results where the bodies touch
along an edge.

The fuzz's membership check now finds a point's side by the parity of a
ray's crossings with the mesh, walking a tree of boxes over its
triangles, and falls back to the winding number where every ray passes
within rounding of an edge or along a triangle; a test checks that both
agree on thousands of points. With the meshes shared by the origins
check, the sphere and cone fuzz's own checks fell from about 8 to 2.4
times the time of its booleans (meshing, at a thousandth of the unit,
is most of what is left).

## 14. Item 3b: planes across cones (design, 2026-10-13)

§12.2 named the conics' tracks by their centres: the parabola
`vertex + a t² x + 2 a t y` and the hyperbola's branch `centre + a cosh t
x + b sinh t y`, with a configuration near a parabola snapped to one so
that no centre is far. Worked through before building, that leaves two
faults. Near the parabola the snap's boundary still allows centres up to
`R² / (eps/16)` away (`R` the region's radius): at a metre, a centre
1.6 × 10⁸ m off, whose rounding (about 2 × 10⁻⁸ m) is a fifth of the
tolerance. And the hyperbola's band grows with `cosh T` while the curve
straightens: the bound is loose exactly where the branch is followed
furthest. This revision replaces them; the rest of §12 stands.

**Literature.** The vertex equation of a conic, `η² = 2pξ − (1 − e²)ξ²`
(`ξ` along its axis from a vertex, `p` the semi-latus rectum, `e` the
eccentricity), one equation for the ellipse, the parabola and the
hyperbola, continuous through `e = 1` (G. Salmon, *A Treatise on Conic
Sections*, 6th ed., 1879; any text on conics). The rational
parametrisation of a conic by the lines through one of its points
(here the vertex), which is the rational quadratic Bézier form of
Piegl & Tiller, *The NURBS Book*, 2nd ed., §7.2–7.3, and E. T. Y. Lee,
"The rational Bézier representation for conics", 1987 (cited in design
note 0009). Apollonius for the section of the cone (§12.2).

**The conic in the plane.** In the plane's frame of §12.2 (`x` along
the axis's direction in the plane, `y` square to the axis, origin `O`
the foot of the cone's reference point `P`, which is kept near the
faces), with `γ = n · u`, `g = |u − γ n|`, `d = (P − Q) · n` (`Q` on the
plane), `r` the cone's radius at `P` and `k = tan α`, the nappe's
section is `A ξ² + η² − 2Bξ + C = 0` with

  A = γ² − g² k²,  B = g (r k − (1 + k²) γ d),
  C = (d g)² − (r − k γ d)²,

each a product or difference of numbers of the size of the faces, so
none loses precision to a far apex. Its discriminant is exact:
`B² − AC = σ²` with `σ = d k − γ r = k δ`, `δ` the apex's distance from
the plane. So `σ = 0` is the plane through the apex; else the semi-latus
rectum is `p = |σ|` and, for `A ≠ 0`, the centre is at `ξ = B/A` with
semi-axes `|σ|/|A|` along `x` and `|σ|/√|A|` along `y`. The plane
meets the nappe (not its mirror) when `A < 0` (the branch at larger
`ξ`, where the radius `r + k h` grows, `h = gξ − γd`), and when `A ≥ 0`
exactly when `σγ < 0` (the axis, from the apex along `u`, reaches the
plane).

**Tracks.** An ellipse whose centre lies within 16 radii of the region's
disk (the plane's section of the region) is the closed track of §10.3,
`centre + a cos t x + b sin t y`; there no coordinate is large and the
whole curve may lie in one face (a cone through a plate at an angle).
Every other section is an open *conic* track in the vertex form

  vertex + x · p t² / (1 + w) + y · p t,  w = √(1 − c t²),  c = A,

with `x` the way it opens from the vertex: the branch of a hyperbola or
a parabola (opening along `+x`), or the half of an ellipse at the vertex
nearer the region (an ellipse much larger than the region, or one near
a parabola). The vertex is the root of `A ξ² − 2Bξ + C` at that end,
taken as `C / (B + |σ|)` or `(B − |σ|)/A`, whichever has no
cancellation, so it is exact to rounding through `A = 0`, where it is
the parabola's vertex: no snap to a parabola is needed. `t = η / p`, so
the parameter resolves positions to the rounding of `η`; its first and
second derivatives are `p (t/w, 1)` and `p (1/w³, 0)`, so `|C''| = p/w³`:
at most `p` on a parabola or a hyperbola, the greatest at the stretch's
`|t|` nearest zero (where `w` is least for `c > 0` it is the greatest
`|t|`). That is the band's bound on a plane's chart and on a cone's
(an isometry, §12.1), now per stretch, tight where the branch straightens.
The first division of a stretch is by equal turning of its tangent,
`φ = atan2(t, w)`, at most an eighth of a half turn each: a conic turns
less than a half turn on a branch or a half ellipse, so at most eight.
A stretch's box is exact: each coordinate is convex or concave in `t`,
its one turning point where `t/w` is a ratio of the frame's components.

**Edges.** With `s = t / (1 + w)` the conic is
`vertex + 2p (s² x + s y) / (1 + c s²)`, a rational quadratic in `s`
(lines through the vertex), so a stretch is exactly a rational quadratic
B-spline: per piece `[s₀, s₁]` the weights `1 + c s₀², 1 + c s₀ s₁,
1 + c s₁²` and control points `vertex + 2p (X x + Y y)/W` with
`X = s₀², s₀ s₁, s₁²` and `Y = s₀, (s₀ + s₁)/2, s₁`. Pieces turn at
most a quarter turn, which keeps the middle weight positive (it is
`1 + c s₀ s₁`, positive on a hyperbola or a parabola always). Far along
a branch `1 + c s²` cancels (`s` nears `1/√|c|`; found by the edge test,
2000 times the rounding), so the weights are computed in `t`:
`1 + c s² = 2/(1 + w)` exactly, and the middle weight is
`N / ((1 + w₀)(1 + w₁))` with `N = 1 + w₀ + w₁ + (w₀w₁ + c t₀t₁)`, the
bracket taken as `(1 − c t₀² − c t₁²)/(w₀w₁ − c t₀t₁)` on a hyperbola
where its terms have opposite signs; the ends are the track's own
points.

**The apex.** A plane within the tolerance of the apex where the apex is
in the region is refused: "a plane through a cone's apex" where it
crosses the cone in two rulings (`A < 0`), and as a tangent contact
otherwise (only the apex, or along a ruling). A plane within `eps/16`
of the apex otherwise (the apex outside the region, as for a frustum
far from its apex) is taken as through it: two rulings (lines, exact on
the cone at the region's height, within `eps/16` of the plane), or a
tangent contact along the ruling `ξ`-ward when the plane is within the
tolerance of the tangent plane there over the region, or apart. Rulings
from an apex in a face, with the apex as an end of both, remain to do.

**Precision.** A conic track's points are found from its vertex, which
for a long branch (a slight taper cut nearly along its axis) can lie
far from the region; rounding is then about `ε` times the vertex's
distance. A track whose rounding could exceed a sixty-fourth of the
tolerance (a vertex some 10⁶ m away at the default tolerance) is
refused rather than computed loosely, with two thousandths to spare
since the fifth review, so that halving the chords of any stretch it
accepts ends under a sixteenth of the tolerance (§16.2, §17).

**Tests.** As §12.4 item 3: the meeting's curves on both surfaces to
rounding and none missed (random cones and planes, planes near a ruling's
direction, near the apex, nearly along the axis), parameters and edges
exact, golden cases for each conic with volumes in closed form (a solid
cone cut by a plane is `(1/3) Σ hᵢ Aᵢ` over its flat faces, `hᵢ` their
distances from the apex, since the conical face's tangent planes pass
through it: the areas of elliptical, parabolic (Archimedes) and
hyperbolic segments are closed forms), and a fuzz of turned solids with
boxes at any angle, with planes set along a ruling's direction.

## 15. As built: item 3b, planes across cones (2026-10-15)

Faces on cones now meet planes at any angle: the closed ellipses, the
open conic tracks of §14 in vertex form (halves of ellipses, parabolas,
branches of hyperbolas, continuous through the parabola with no snap),
and the rulings of a plane through the apex of a frustum whose apex lies
outside both bodies. The tracks, their exact stretch boxes, bands and
rational quadratic edges are §14's. What changed while building it, each
handled and tested:

* **Rays, not lines.** A plane through the apex within `eps/16` crosses
  the cone (one nappe) in two rulings that start at the apex. As lines,
  their halves on the mirror nappe cut a face the apex lies in (a slab
  through the apex of a frustum: the slab's face was cut along a line
  that does not lie on the cone, and the boolean failed as "meets the
  other body's surface inconsistently"). They are rays now,
  `origin + t dir` for `t ≥ start` with the apex at `start`, `dir` away
  from it. Their origin is on the cone at the height of the region's
  centre; where that centre lies below the apex (beside an oblique axis,
  the apex still outside the region), at the same distance beyond it
  instead, since below the apex is the other nappe (found by review: the
  rays ran down the mirror nappe and every operation was refused as
  meeting at too shallow an angle; a unit test now puts regions there).
* **Every crossing on a track.** A cut, where a face's boundary crosses
  the other surface, must lie within 16 tolerances of exactly one track,
  or the boolean refuses ("too shallow an angle", "nearly touch"). The
  region bounds where faces on *both* surfaces can be, but one face's
  boundary can cross the other surface well outside it, where that face's
  box reaches further. So both rays are returned when either comes near
  the region (a slab through a frustum's apex along one ruling: the
  frustum's base crossed the other), and an open half ellipse comes with
  its far half, from its other vertex (a tool 0.1 wide on a frustum 8
  high: the base and seams crossed the far half); both were refused
  until then (found by probing the fuzz's refusals, and by review).
* **The apex against the faces' box.** The region's ball holds the
  common box of the faces on both surfaces with room to spare, and the
  apex rules asked whether the apex was in the ball: a countersink's
  apex, in the hole a little below its cone's face, was, and every half
  section through the hole's axis was refused as a plane through a
  pointed cone's apex. The rules now ask whether the apex (or a small
  circle about the axis) is in the common box, where faces on both can
  be (`Region` keeps the box). The box is square to the world's axes,
  and beside a leaning axis it can hold the apex of a frustum whose
  faces stop well short of it (a frustum leaning along (1, 1, 1), its
  apex half a unit beyond its top, cut by a plane through the apex:
  refused, found by review). So the region also keeps the least height
  of the cone's faces above its apex along its axis (`Region::lift`),
  from the control points of their edges' Bézier pieces: a face on a
  cone lies in the convex hull of its boundary (each of its rulings ends
  on it), and that in the hull of those points (exact for the circles
  square to the axis that bound a frustum). The apex rules apply only
  where that height is within their margin too. A plane through an apex
  out of reach gives rays, from the apex where that is near the region.
* **Rays from the region, not the apex.** A ray now starts where it
  comes within `CLIP` = 256 tolerances of the region (`Track::Ray`'s
  `from_apex` is false), or at the apex where the apex is that near; a
  ray that stays farther starts at infinity. Cuts on a ray's line before
  its start are skipped: no point there is in faces on both surfaces,
  and near the apex they lie on both rays, which refused every section
  of a frustum by a box whose edge runs along the axis past the apex (a
  section view of a countersink, whose apex lies in the hole). The face's
  side at the start is found by locating the start in it, so a cut
  within `CLIP/2` tolerances of the start, where it could lie on either
  side, is refused, as is a start on a face's boundary (a pointed
  cone's seam running along the ray).
* **Open tracks that end inside a face.** A ray's apex, and a half
  ellipse's ends on its minor axis (`t = ±1/√c`), lie beyond the region
  where both faces can be, but can lie inside one face. The cuts of such
  a face no longer pair up: its spans gain an end at `±∞` (a sentinel
  cluster, no point of the result) where the track's end lies inside it,
  decided by locating the end in the face. A stretch with such an end
  inside both faces cannot occur and is refused as an invalid result; an
  end on a face's boundary within the tolerance is refused ("a curve
  where faces meet ends on an edge"). The same rule made a half ellipse
  wholly inside a slab 1,000 wide work (it has no cuts at all, and was
  refused). Where a face's ends then fail to alternate, the refusal
  names a cut's point, never the sentinel's (review found it indexed
  the points with the sentinel and panicked).
* **A plane near a pointed cone's apex.** A plane `δ` from the apex
  meets the cone no nearer its axis than `k |δ| / (|γ| + k g)` (`k |δ|`
  for a circle, `|δ|` for a plane along the axis; `γ` and `g` the cosine
  and sine of the axis's angle to the plane's normal); when that comes
  within the hundred tolerances by which validation keeps edges off a
  surface's singular points, the edge it leaves on a face reaching the
  apex is refused by validation, or the face's mass properties do not
  converge (planes 0.5–0.9 µm from the apex in the fuzz). Such a plane,
  its curve within 128 tolerances of the axis, where the apex is within
  `128 eps √(1 + k²) / k` of the region and the plane crosses the nappe
  (`A < 0` or `σγ < 0`), is refused: "a plane passes so near a cone's
  apex that their curve comes within the tolerance of its axis". The
  rule first held `|δ|` itself to that bound, the largest over all
  angles: planes along a slender cone's axis 0.2–0.35 mm from it (2°
  half-angle), their curves thousands of tolerances off the axis, were
  refused (found by review). With a plane within `eps` of the apex there
  (two rulings through it), these are the remaining apex refusals.
* **Circles about one axis within the tolerance of it.** Surfaces about
  one axis meeting in a circle of radius between `eps` (a point contact)
  and 128 tolerances, where faces on both can be (a box's face 16 µm
  below a pointed cone's apex), are refused: "surfaces about one axis
  meet in a circle within the tolerance of the axis" (found by review:
  the result's circle passed the cone's apex within the validator's
  hundred tolerances, and the kernel's own check of the result caught
  it as invalid).
* **Edges through an axis, refused by name.** Any result whose edge
  passes within the validator's hundred tolerances of the axis of its
  face's surface of revolution between its ends is refused before the
  result is checked, with a message naming a turned face's centre, a
  cone's apex or a sphere's pole (`assemble::check_axes`, with the
  validator's own test). A box whose face runs 0.5–5 µm off a turned
  rod's axis cuts its flat ends so; the kernel's check of its result
  caught it, but as an invalid result, a bug's error.
* **A plane nearly square to the axis** was taken as coaxial (a circle)
  when its normal's sine times the cone's radius at the axis was within
  the snap; near the apex that radius is small while the section, at a
  slant, is far from a circle (a plane through the apex crosses the cone
  in two rulings, not a point). The section of a plane tilted by `θ`
  lies within `2 ρ₀ tan θ (1 + k)` of the circle of radius `ρ₀` where it
  crosses the axis while `k tan θ ≤ 1/2`; that bound, not the radius, is
  held to the snap now.
* **Meshes of faces cut across their rulings** grew to tens of thousands
  of triangles (design note 0002 §3, 2026-10-15): a ruled face is now
  triangulated along its rulings.
* **Matching at cusps never settled** where three faces' edges leave one
  vertex together, as at the point where a box's edge pierces a cone:
  design note 0002 §5, "Groups of cusps".
* **Planes nearly tangent along a ruling, by the apex** (found by review
  of item 3b). Along the cone's ruling nearest a plane through its apex
  within the tolerance, the cone leaves the plane by `μ = cos α |A| /
  (|γ| + k g)` per unit of length (`α` the half-angle), so it stays within
  the tolerance of the plane out to `reach = (eps + |δ|) / μ` from the
  apex, the further the more nearly the plane is tangent. Where faces on
  both bodies can come within that reach of the apex (the region's box,
  and the faces' least height above the apex along the axis), the
  tolerance cannot decide how they meet: two rulings are refused, a
  plane on one side is a touch. The apex rules had been skipped wherever
  the cone's faces stop short of its apex by more than the tolerance (a
  frustum's top a thousandth of its height below it), and a plane that
  nearly touched the top rim made a union of two solids a nanometre
  apart, or a cut the other body's surface met inconsistently.
* **A rim touching a plane off the apex.** The nearest ruling runs from
  the apex along `r̂ = cos α u − sign(γ) sin α n⊥ / g` (`n⊥ = n − γu`),
  leaving the plane at `m = sign(γ) cos α (|γ| − k g)` per unit of
  length: `f(s) = δ + s m`. Where `γδ > 0`, every other ruling leaves
  toward the apex's side faster, or comes toward the plane more slowly,
  so the cone lies on the apex's side of the plane wherever this ruling
  does, and is nearest the plane on it. Along it, from the faces' least
  height above the apex, through the region's box: where it stays on
  the apex's side and comes within the tolerance of the plane, the
  surfaces touch (a tangent contact, refused). A frustum's rim within
  the tolerance of a plane nearly tangent along this ruling had given
  two solids a few hundredths of a micrometre apart.
* **Touches read along whole loops.** A face's boundary that runs within
  the tolerance of the other surface without crossing it (a run of
  stations on the surface between stations on one side) is a tangent
  contact and refused; the stations are read along each loop as a whole,
  across its vertices, where edge by edge a run spanning a vertex was
  missed. A run that holds a whole edge is that edge lying on the
  surface (a rim on a plate), not a touch. Nor is a run through a vertex
  that another edge leaves to the other side, by its first station past
  the vertex: there the body crosses the surface and the curve where they
  meet passes through the vertex (a box's corner on a cylinder with two
  edges going in and one out; a cylinder crossing a turned rod along its
  seam, touching one half of the rod's cap at the seam's end and entering
  the other). Read face by face, on an intermediate version, these had
  been refused as touches: 202 operations of 40,000 in the
  planes-and-cylinders fuzz at 50 times, every one of which the code
  before the rule did correctly (found by diffing the fuzz's refusals
  against it). Read along whole loops, the final code refuses 52
  operations there that the code before the rule did, and 11 of 20,000
  in the fuzz of solids at an angle (counted operation by operation
  against the code before item 3b): edges that touch a face without
  crossing it, genuine tangent contacts, which the code before the rule
  got right but which the rule refuses with the near-tangent planes that
  had given two solids a nanometre apart. They wait for tangent contacts to
  be handled, which comes after operands with conic edges (STATUS.md).
* **Closed tracks with no cuts.** Whether such a track lies in both faces
  was decided at one point of it, which could lie where only one of the
  faces is (an ellipse grazing a slab's edge metres from the frustum, an
  invalid result found by review). It is decided now at its quarter
  points within the box where faces of both bodies can be and at its
  point nearest the box's centre, and at all four quarter points when
  none of those decides.
* **Meshes of the results** (design note 0002 §4): slivers too small to
  split that face against the surface by an apex are unfolded by edge
  flips, or else split at their centroid, and needles at a pointed
  cone's apex whose corners lie in a line within the model tolerance are
  not tested for a side; both made the kernel's own check refuse valid
  results as "too small to split" (23 of 300 unions of a frustum with a
  slab nearly tangent to it by its apex, before the centroid split).
  Review measured meshes against the surfaces' closest points and found
  them far past the chord tolerance by apexes and poles (16 tolerances
  by a pointed cone's apex), and a sliver standing on edge across a cone
  0.45 mm high in a validation mesh; meshing was rebuilt to hold the
  tolerance, as estimated between the points it measures, and judges
  facing at the corners too (0002 §4).

**Tests.** Golden cases (`crates/kernel/tests/curved.rs`), each with all
four operations, volumes exact to 1e-9 against derivations from the
rulings and chords (the area of each flat face cut from the cone in
closed form: elliptical, parabolic and hyperbolic segments), topology,
validity, origins and one solid of genus 0: a pointed cone under a
roof tilted 20° (a whole ellipse), a countersink leaning 15° into a
plate, a cone cut along a ruling (a parabola) and a thousandth of a
radian either side (a long ellipse, also by a slab 1,000 wide, and a
hyperbola close to a parabola), a cone cut parallel to its axis (a
hyperbola, its area also integrated), a slight taper (1°) cut along its
axis (also by slicing), and a frustum cut by a plane through its apex
(two rays). Then, for the rules above: sections of a frustum by boxes
whose edges run along its axis past the apex, a micron either way;
halves and a quarter of a countersunk plate through the hole's axis;
frustums leaning along (1, 1, 1) cut through their apexes, at a slant
and along the axis; a box beside a leaning frustum in a plane through
its apex; small tools across large cones (both halves of an ellipse); a
slender cone (2°) cut along its axis 0.2–0.35 mm off it, exact, and a
box 13 µm off refused by name; a box across a cone just below its apex,
square and tilted, refused by name; and boxes a hair off a turned rod's
axis, exact or refused by name. Unit tests in `meet.rs` put every track on both surfaces to
rounding and miss none, near the parabola, the apex and the axis. A
fuzz of tapered solids (cones pointing up and down, frustums, tubes with
tapered bores, double cones) about turned axes, each with a box turned
any way, a slab along a ruling's direction (tilted by 0 to 0.1 rad, often
by 10⁻³ to 10⁻⁹), a slab through or near the apex, or a slab within
10⁻⁴ to 10⁻⁷ rad of tangent along a ruling, through the apex or within
a micrometre of it. At 50 times its default size (5,000 pairs, 20,000
operations; 2026-10-15, on the final code) every result is valid with
its origins on its operands' surfaces, every volume identity holds
where the results exist, and point membership agrees. 5,672 operations
are refused (284 in 1,000), none with a box or a slab across the
rulings. Through or near the apex, 2,856, all on pointed cones: 2,256
planes through the apex where both bodies reach it (two rulings), 596
whose curve comes within the tolerance of the axis, 4 faces nearly
touching. Nearly tangent along a ruling, 2,596 on pointed cones (724
two rulings, 1,088 near the axis, 784 faces touching) and 220 on
frustums, where the cone stays within the tolerance of the plane, or
between two rulings closer than it, out to the frustum's faces, so the
tolerance cannot decide how they meet (the near-tangent bullets above):
96 refused as two rulings, 100 as tangent contacts, 24 as slivers,
shallow angles or a face running within the tolerance of the other
surface. The other fuzzes at scale on the final code, against the code
before item 3b: polyhedra at 100 times refuse the same 47 of 120,000
operations, and the same in the tolerance motions; solids about one
axis the same 1,122 of 20,000; planes and cylinders at 50 times 879 of
40,000 (827 before) and at an angle 9,609 of 20,000 (9,598, quartics),
the additions (52 and 11) all edges touching a face (above), and none
dropped.

**Left refused.** A plane through a pointed cone's apex, or whose curve
comes within 128 tolerances of its axis, where the apex is near where
faces of both bodies can be: the result's edges would meet at the
cone's singular point, which the validator keeps edges away from
(nearly every section of a pointed cone through its axis). A second
boolean on a body whose conic edge lies on a curved face: the edge is
stored as a rational quadratic NURBS curve, which the stages do not yet
recognise as a conic, so they refuse it even where the bodies are apart
(since done: §16, §17). Cones with cylinders, cones or spheres off their
axes (quartics, §10.3).

## 16. Operands with conic edges on curved faces (design, 2026-10-15)

A boolean whose result keeps part of a conic (§14, §15) stores that
edge as a rational quadratic NURBS curve: a half ellipse in pieces of at
most a quarter turn with weights `1, cos(δ/2), 1`, a parabola or a
branch of a hyperbola in at most two pieces. A second boolean on such a
body is refused today, even where the bodies are far apart, by one
rule: tracing a NURBS edge in the chart of a curved face (§12.1) is
"not supported". It is reached from four places: the outline of every
curved face of both operands (before any pair of faces is looked at),
the split of each face, the merge of faces and the check of solids with
voids (§10.6). No other stage asks what kind an edge's curve is: the
splits of an edge by a surface (§10.2), the crossings, the clustering of
points on edges and the validation work on any rational curve of degree
at most five.

### 16.1 Choice: follow the stored curve

Three ways were weighed (ADR 0023). Recognising the conic a NURBS edge
came from, by its control points and weights (E. T. Y. Lee, "The
rational Bézier representation for conics", 1987; Farin, *Curves and
Surfaces for CAGD*, ch. 12), is ill-conditioned where it matters: on a
short or flat piece the conic's shape factor differs from that of a
parabola by the square of its turning, the axis of a nearly circular
ellipse is undefined, and a recognition that is wrong by more than a
band lets the split misplace a crossing while the result's edge, still
the right NURBS, passes validation. Storing exact conics on edges (new
curve kinds) touches every place curves are evaluated, sampled, boxed
and meshed, and is better done with exchange formats, which need it
anyway. The stored NURBS is already exact, so the stages follow it
directly: its course in a cylinder's or a cone's chart, with a band
(its distance from its chords) proven from the curve's own
derivatives. No conic is recognised and the NURBS stays the one record
of the edge. The same holds for any rational curve of degree at most
five lying on a cylinder or a cone.

### 16.2 The band of a rational curve in a developable chart

A cylinder's and a cone's charts are isometries of the surface onto the
plane (§12.1). The image `Y = φ ∘ X` of a curve `X` on the surface has
for its second derivative the part of `X''` along the surface (the
covariant derivative, which an isometry keeps; do Carmo, §4-4), so
`|Y''| ≤ |X''|`. Over a rational Bézier piece `X = P / W` of degree `p`,
clipped to the stretch of a chord and parameterised by `s ∈ [0, 1]`,

`X' = R / W²`, `R = P'W − PW'` (degree `2p − 2`),
`X'' = Q / W³`, `Q = R'W − 2RW'` (degree `3p − 3`),

computed as products of degrees `2p − 1` and `3p − 2` (for a conic
piece `Q` has five coefficients, well within the polynomials of
§10.7), with `P` taken relative to the piece's first point: `R` and `Q`
do not change, and the numbers they are computed from become small. A
Bernstein polynomial lies in the convex hull of its coefficients, so
where every weight `wᵢ` is positive (`W > 0` alone is not enough) `|X'|
≤ max |Rᵢ| / (min wᵢ)²` and `|X''| ≤ max |Qᵢ| / (min wᵢ)³`, the `wᵢ`
those of the clipped piece. The kernel's NURBS curves have positive
weights and subdivision keeps them so. Linear interpolation of `Y`
between `s = 0` and `s = 1` errs by at most `sup |Y''| s (1 − s) / 2`,
so the chord is within `sup |Y''| / 8` of the course. The bound is
loose where the curve's speed changes along it (`X''` then has a large
part along the curve, which does not bend the image). Measured on
random plane sections of cylinders and cones written as NURBS (2026-10-15,
the unit test of §16.5): the band is 4.7 times the farthest a dense
sample strays from its chord at the median, 140 at the 90th percentile,
`1.1 × 10⁴` at the 99th and `7.4 × 10⁵` at worst, far along a
hyperbola's branch; review's simulation of 7,500 sections (half-angles
2° to 89°, planes down to 10⁻⁸ from the apex) found up to `1.5 × 10¹¹`
on final chords. Looseness costs chords, never a wrong answer. The
division (below) halves a chord while its bend is more than a quarter of
its length; a bend shrinks with the square of a chord's length and the
length only linearly, so a band `L` times too loose takes about `log₂ L`
more halvings, where the curve is nearly straight. Review's sections
needed at most 116 chords a course, of the 4,096 allowed. After the
division the planar graph's refinement (§11) quarters a band a round
where the band comes from derivatives, so it needs about `log₄ (chord /
(ε / 16))` rounds, about 13 for a chord of a metre at the default
tolerance, of its 60. A band from control points (a NURBS course in a
plane's chart) can narrow more slowly: review measured a halving only
about halving it along a hyperbola's arm with very uneven end weights.

**Off the surface.** A chart maps a point near its surface to the image
of its foot: on a cylinder of radius `r`, `(r θ, h)`; on a cone of
half-angle α, slant position `ρ sin α + z cos α` and angle `θ sin α`,
in cylindrical coordinates `(ρ, θ, h)` about the axis, `z` the height
above the apex. An edge lies on its face's surface only within
rounding, so its image's second derivative has terms in its distance
`δ` from the surface: `δ = ρ − r` on a cylinder, `δ = (ρ − R(h)) cos α`
on a cone, `R(h) = r₀ + h tan α` its radius at height `h` along the
axis from the chart's point (`= ρ cos α − z sin α`). Differentiating
the chart along `X`:

* cylinder: `(r θ)'' = (r/ρ) (X''·e_θ − 2 θ' ρ')`, `h'' = X''·a`, with
  `ρ' = δ'`;
* cone: radial part `s'' − s ψ'² = X''·e_s + δ sin α cos α θ'²`,
  tangential part `s ψ'' + 2 s' ψ' = (1 − δ cos α / ρ) X''·e_θ − 2 cos α
  θ' (δ' − δ ρ' / ρ)`, with `e_s` the unit ruling direction and `s sin α
  = ρ − δ cos α`.

With `|θ'| ≤ |X'| / ρ` and `|ρ'| ≤ |X'|`, both give

`|Y''| ≤ (1 + |δ|/ρ) |X''| + 3 |δ| |X'|² / ρ² + 2 (1 + |δ|/ρ) |δ'| |X'| / ρ`,

which is `|X''|` on the surface. (Review found the constants generous:
the cylinder needs no `3 |δ| |X'|² / ρ²` term, the cone 2.5 in place
of 3; they are kept for one formula.) `δ` and `δ'` come from the
surface's implicit function along the piece (§10.7: `G = n / W²`, `n`
a Bernstein polynomial of degree `2p`): on a cylinder `G = ρ² − r² = δ
(ρ + r)`, on a cone `G = ρ² − R(h)² = (ρ − R)(ρ + R)` and `δ = (ρ − R)
cos α`; `G' = (n'W − 2nW') / W³`. With `g = max |nᵢ| / (min wᵢ)² ≥
|G|`, `|δ| ≤ g / (ρ_min + r)` (cone: `g cos α / (ρ_min + R_min)`), and
`|δ'| ≤ (|G'| + |δ| |X'|) / (ρ_min + r)` (cone: `(|G'| cos α + |δ| sec
α |X'|) / (ρ_min + R_min)`, as `|(ρ + R)'| ≤ |ρ'| + tan α |h'| ≤ sec α
|X'|`, `ρ'` and `h'` being orthogonal parts of `X'`). `R_min` is the
least of `R` at the piece's control points, whose hull holds it.
`ρ_min` is the larger of two bounds: `ρ` is 1-Lipschitz in space, so
`ρ ≥ ρ(X(½)) − ½ max |X'|`; and `|ρ − R| = |G| / (ρ + R) ≤ |G| / R`,
so `ρ ≥ R_min − g / R_min` (`R − g / R` grows with `R`). Review showed
the first alone too weak: along a long ellipse on a thin rod (a plane
nearly along the axis) it was not positive on most pieces, and a whole
ellipse 1,400 times longer than the rod's radius would have needed more
chords than a course may have. Where `ρ_min` or `R_min` is not positive
(a stretch reaching toward the axis or the apex) the band is infinite
and the chord is halved; a chord whose band is still infinite after the
halvings is refused with a message naming the axis or apex. The bounds
hold however far the piece lies from the surface, so an edge off its
face needs no rule of its own (the kernel's operands are its own
results, whose edges it has checked lie on their faces).

**Rounding.** Every bound is computed with its rounding allowed for, so
that it is a bound on the curve as stored, not only in exact
arithmetic:

* The control points' offsets from the reference point (on the axis
  near the piece, §10.7) round at a few ε of the coordinates they are
  computed from, not of the offsets' size, and the reference itself (and
  a cone's radius there) rounds off the stored surface. Each
  coefficient of `n` carries `2m (1 + tan α)² (reach + |R| w_max) + 2m²
  (1 + tan α)²`, how far a quadratic form of the offsets can move when
  each moves by `m = 8ε (max |Pᵢ| + (|X₀| + |c| + r₀ + |c'| + |R'|)
  w_max)` (`c`, `r₀` the stored point and radius, `c'`, `R'` the
  reference and its radius), and `10⁻¹²` of `n`'s size for its own
  products; `n'W − 2nW'` carries `4 (deg n + 1)(max |wᵢ| + max |w'ᵢ|)`
  times that.
* `R` and `Q` carry the rounding of `P`'s coordinates (4ε of the numbers
  they come from) and `10⁻¹²` of `P`'s coefficients, through the same
  products by the magnitudes of their terms (Higham, *Accuracy and
  Stability of Numerical Algorithms*, §3.1).
* Lower bounds (`ρ_min`, `R_min`) lose 8ε of the numbers they are
  computed from, and the final bend is raised by `10⁻¹²` relative for
  the formula's own few dozen roundings of positive numbers.

Review measured, in exact arithmetic on the stored curves, the distance
from the surface up to 31 times the bound computed without the first
allowance (on a near-apex piece), and the acceleration `6 × 10⁻⁶` over
its bound without the second; both far inside the slack below, but not
bounds. A test checks the bounds against references in double-double
arithmetic at samples with no tolerance, and fails without them. What
is left is positional: the chart's points at the chords' ends, and the
Bézier pieces' control points, which knot insertion and de Casteljau's
subdivision make by convex combinations (the weights being positive),
each within a few `pε` of the coordinates. Moving the curve by `η`
moves its course by about `η` (the chart's Lipschitz constant, about
`1 + |δ| / ρ`, is close to one for an edge on its face), so to each band
is added `10⁻¹³ (1 + size + the chart's size)`, as for tracks (§12.1):
the curve's size, and the size of the numbers the chart subtracts and
scales by (a plane's origin; a cylinder's or a cone's stored point plus
its radius; a sphere's centre plus four radii). The curve's size alone,
as first built, missed a chart anchored far from the curve: measured
against double-double references, a conic's chart points strayed up to
26 times that allowance with the chart's anchor 10⁴ m from the conic's
vertex (§16.5).

**Division.** The first chords of such a course are, as on a plane,
two per Bézier piece, so that no chord spans a knot: the conic edges
are only continuous in their parameter at their double knots (the
speed jumps there by a factor of 0.01 to 30, review measured), the
ellipses' once differentiable; a bound from `sup |Y''|` needs `Y` twice
differentiable between the chord's ends. Then every chord whose bend
(its band less the positional allowance below) is more than a quarter
of its length is halved, a round at a time, up to the 4,096 chords a
course may have or until the bend is down to four times its floor
(below), so that the outline's widest band, which bounds the
searches of point location, stays small. Rounds spend the chords evenly
along the course when the cap is reached; halving depth first, as first
built, could spend them all on an early stretch and leave a later one
too wide, or unbounded and so refused (review; not seen). A band asked
over a range that spans knots is the largest of its pieces' plus the
farthest the course is, at an interior knot, from the whole chord: each
piece's chord differs from the whole chord by a function affine in the
parameter, so most at the piece's ends.

**Halving and rounding.** Every band is a sag, which halving the chord
narrows, plus the positional allowance above, which nothing narrows. The planar graph's refinement (§11), the search for a
point inside a piece and point location halve a chord while its sag, not
its band, is over a sixteenth of the tolerance. Halving by the band never
ended where the allowance was over that, more than 62.5 km from the
origin at the default tolerance (6,250 km for a NURBS course on a plane
before this item raised its allowance, §16.4), and the search for a
point inside a sliver then doubled its chords every round until memory
ran out (review).

The sag has a floor of its own, the part that its terms' allowances for
their own rounding make, which halving does not narrow either:

* A NURBS course on a cylinder or a cone: an eighth of the acceleration
  bound's allowance for the control points' rounding (the longest error
  vector of `Q`'s coefficients over the least weight cubed, times the
  off-surface factor `1 + δ/ρ`), carried with the bound. Review found it
  about `4.3 × 10⁻¹⁴` of the coordinates for evenly weighted conics, 43
  times what a draft of this note estimated (`10⁻¹⁵`); the unit test
  below measured up to about `2.7 × 10⁻¹³` with weights uneven by up to
  a factor of 11.
* A NURBS course on a plane: the allowance its bound carries for the
  rounding of the distances it is computed from, `ε ((7p + 2) R + (6p +
  3.5p (1 + ρ)) D + 6.4 d)` for degree `p`, `R` the control points'
  largest distance from the world's origin, `D` their diameter (the
  diagonal of their box), `d` their largest distance from the chart's
  origin and `ρ` the ratio of the curve's largest weight to its least.
  The weights cancel from the arithmetic of de Casteljau's steps: with
  positive weights and coefficients every term a step sums is at most
  `R` times its weight, and the division by the weight undoes it,
  `(4p + ½) R` for the two runs of `p` steps a clip takes and `p D` for
  rounding `1 − s`. They enter through the clip's parameters, rounded
  before `1 − s` is taken from them: moving a split by `δ` moves a point
  by up to `(1 + ρ) D δ`, at most `2p (1 + ρ) D` for the clip. Splitting
  the curve into its pieces adds `2p R + 2p D`, and `1.5p (1 + ρ) D` for
  its parameters (sketched, not derived; measured a hundredth of it);
  the chord's ends, evaluated, `(p + 3/2) R + 3p D`; the chart and the
  distances `4.5 √2 d`. Against double-double references (§16.5), over
  ten million stretches of rational curves of degrees two to five with
  weights' ratios up to 10⁴, 1 m to 10¹⁰ m from the origin and down to
  10⁻¹² of a piece, the rounding was at most 0.47 of the allowance.
  First built as `32 ε ρ (1 + size)`, the weights' ratio times the whole
  size: it held the sags of the hyperbolas of cones cut a millimetre off
  their axes (weights' ratios in the hundreds) over a sixteenth of the
  tolerance from 1 km from the origin, refusing their second cuts there
  (third review), though their rounding grows with their own size, not
  their distance; and it was not carried in the bound, whose rounding
  then only the positional allowance covered, and too little of it where
  `ρ D` is large (measured up to `1.2 × 10⁴ ε (1 + size)`).
* A conic track: the allowance `32 ε` of the stretch's reach from the
  vertex that its sag carries. To first order the rounding of a conic's
  chart points across its course is at most `7.5 ε` of the reach in a
  plane's chart and `27.2 ε` in a cone's, plus a part that grows with
  the vertex's coordinates, its distance from the chart's anchor and a
  cone's radius there, which the positional allowance covers; measured
  against double-double references, at most `1.9 ε` of the reach over
  1.4 million points of conics cut from cones by planes (§16.5). First
  built as `10⁻¹³` of the reach (450 ε), which refused the first cut of
  a cone 20 km tall by a plane parallel to its axis (third review).
* Lines, circles and ellipses, whose bounds narrow to nothing: none.

As its chord shrinks a sag comes down to its floor: on stretches a
billionth of a course's parameter, every kind measured at most 1.21
times its floor. So every halving
stops at four times the floor as well as at a sixteenth of the
tolerance: a chord is halved while its sag is over both. A conic
track's bound settles on its floor itself (its other term is a product
of the stretch's parameters, rounded only in proportion to it, and
narrows to nothing), so its chord is halved while its sag is over the
sixteenth and over its floor by more than a thousandth of it, about six
halvings past four floors; stopped at four floors as the others are,
chords from 220 km of reach could stop over the sixteenth, and from
500 km always did, refusing booleans that more halving resolves
(fourth review, §17). A chord whose
sag is held over the sixteenth that way, where it could hide a contact
(the refinement) or where it leaves a piece's point within the
tolerance of the piece's boundary (the search), refuses the boolean,
naming the place: the bodies lie too far from the origin for the
tolerance. Point location takes a point within such a chord's band as
on the boundary. Without the floors, past the distance where a floor
passed a sixteenth of the tolerance (about 145 km for an evenly
weighted conic on a cylinder, review found), point location halved
every chord within its band of a point down to the parameter's last
bit, the chords doubling at each depth (some billion bends for one
point), and the refinement doubled the chords at a node every round;
the search refused a round that did not halve its widest sag, a rule
that could also refuse a refinement a round from done (a conic track's
allowance between `ε/24` and `ε/16`; review), and is gone. At the
default tolerance four floors first pass the sixteenth about 36 km from
the origin for an evenly weighted conic on a cylinder (a rod is cut
twice 30 km from the origin and refused at 300 km; §16.5); a conic
track's floor would pass it at about 880 km of reach, but meet refuses
the conic's vertex as too far away to compute from about 878 km of its
reach, the vertex's coordinates, its distance from where the faces meet
and that region's size (§14), which bounds every stretch's: two
thousandths short of the floor's limit, so that the halving, settling
within a thousandth of the floor, always ends under the sixteenth
(fifth review, §17; a cone at the origin 400 km tall is cut by a plane
parallel to its axis, one 600 km tall refused, a first cut); and
for a small NURBS curve on a plane, by its allowance's formula, from
about 160 km (degree five) to 310 km (a conic) from the origin where its
chart's origin is near the world's, farther where it is near the
curve.

**Spheres.** No conic lies on a sphere (a plane meets it in a circle,
stored as a circle), so a NURBS edge on a sphere's face stays refused,
with a message naming the sphere.

### 16.3 Coincident curved segments

Two curved segments between the same two clusters were taken as one
when their middles were within the tolerance, the middle of each being
its course's point at the mean of its end parameters. For circles every
course's parameter is the angle, so the middles agree; a conic's NURBS
parameter and its track's parameter (§14) do not: review measured the
two middles of one stretch 0.9% to 4.4% of its chord apart at the
median, by the kind of stretch, and up to 49%. A second cut in the
plane of an earlier one (the same conic from both bodies) would then
make two segments with one course, and the split would refuse the
result.

Nor is one point enough, by any definition of the middle. Conics
through the same two points have three more degrees of freedom; review
built a half ellipse and a parabola in one plane between the same ends
(5 × 10⁻⁸ apart), symmetric about one axis, so that their middles by
parameter, by distance and by turning all fall at their vertices,
within the tolerance, while they are 2,548 tolerances apart elsewhere.
(A draft of this note gave another reason, that the point of a course
equidistant from its ends need not be unique on an ellipse's arc.
Review showed it is unique: the difference of the squared distances to
the ends is affine in the point, so a sinusoid along an ellipse, with
one zero on each arc between the ends; the same holds on hyperbolas and
parabolas.)

So two curved segments between the same clusters are one when three
points of each course, a quarter, half and three quarters of the way
along it, are within the tolerance of the other's stretch (their feet
within it), and so are three more, at the same fractions of its
tangent's turning, on an ellipse, a conic or a NURBS course (a circle's
parameter is its angle already). A conic track's turning has a closed
form (§14). On an ellipse track or a NURBS course the stretch is broken
into pieces that each turn less than a half turn, an ellipse's at most a
quarter of its parameter's period (its tangents half a period apart are
opposite), a NURBS course's at its knots, and a knot span's piece
halved until each part is shown to: the tangent of a rational Bézier
piece with positive weights is a combination with nonnegative
coefficients of the differences of its control points (by the quotient
rule), so a piece all of whose differences lie on one side of a plane
through the origin turns less than a half turn. A rational quadratic
piece of a conic passes at once (its legs are not opposite); a cubic
can turn further, and the third review wrote an arc of 238° of an
ellipse as one cubic span with positive weights, whose probes crowded
as before. The angles between each piece's end tangents are summed,
and the fraction found by bisection in its piece. The angle from the
stretch's first tangent alone falls back after a half turn, which an
ellipse's stretch can pass: the probes crowded into its first part
(review). Each set alone has a blind spot: on a
hyperbola's stretch across its vertex the probes by parameter (linear
across the axis on a track) can all fall on the arms, and the probes by
turning crowd the vertex. Review built pairs of different coplanar
conics passing within the tolerance of one set's probes and measured
the farthest the stretches are apart over the farthest any probe is
from the other course: by parameter alone up to `4.4 × 10⁵` for two
tracks, 10.6 for two NURBS edges, 9.7 for an edge and a track; by
turning alone up to 143; by both at most 2.3. For two plane sections of
one cone, 1.4 by parameter alone. The probes are taken both ways; from
one course only, up to `3.1 × 10⁵`.

A face lying on the other body's face (coplanar) takes the edges of the
other face that lie in it. Whether a curved one does is asked at the
point of its course as far from one end as from the other (found by
bisection on that affine difference), at least half the chord from
each end; the parameter middle can lie within `1.5 × 10⁻⁴` of the
stretch from an end (review), where a shallow crossing of the face's
boundary leaves point location undecided and the edge was kept.

### 16.4 Smaller changes

* The search for a point's foot on a NURBS edge keeps to the edge's
  stretch of the curve (it searched the whole curve, then dropped feet
  outside the stretch, which loses a point on a stretch whose curve
  comes back near it).
* The check of solids with voids builds the outer shells' outlines once,
  and only when a void is found, not once per void.
* A NURBS course in a plane's chart gets the tracks' positional
  allowance, `10⁻¹³ (1 + size + the chart's size)` (§16.2), instead of
  `10⁻¹⁵`: it has not run in a boolean before (every plane face bounded
  by a conic also borders a curved face, which refused first). Every
  course's allowance now counts the chart's size as well as its own.
* The foot of a point on half an ellipse is also found by the ellipse's
  eccentric angle, and the nearer of the two taken: near the ends of the
  half, on its minor axis, the conic's parameter runs ever faster (`ξ' =
  p t / w`, `w → 0`), and Newton's method in it stopped short of the
  foot, 360 times as far from a point 10⁻⁸ from an end (review). It
  decides which curve a crossing lies on and which clusters lie on a
  curve.
* `Kernel::edge_kinds` gives the kind of each edge's curve (a line, a
  circle, or a NURBS curve: free-form or a conic other than a circle),
  so that tests can count the conic edges they exercise.

### 16.5 Tests

* Unit: bands of rational edges hold on random plane–cone and
  plane–cylinder tracks written as NURBS, over random sub-ranges either
  way, against dense samples, and with control points pushed off the
  surface; the derivative and off-surface bounds hold at dense samples
  of random pieces of degrees two to five with uneven weights; the
  bounds hold against references in double-double arithmetic (Dekker,
  "A floating-point technique for extending the available precision",
  1971) on thin cylinders and cones far from the origin, with no
  tolerance, and fail without the rounding allowances; a NURBS edge on
  a sphere is refused, naming the sphere; the foot of a point near the
  ends of half an ellipse; two hyperbolas apart only at their vertices
  make two segments; a sag that halving does not narrow is refused;
  halving every kind of course stops before the parameter's last bit,
  its sag below four times its floor (a conic track's within a
  thousandth of it): NURBS ellipses on planes and cylinders, circles on
  cones, conic tracks (some stretches 10⁶ m and more from their
  vertices) and rational Bézier pieces as large as their distance from
  the origin with weights' ratios up to 10⁴, from the origin to 10⁸ m
  from it, and the hyperbolas of cones cut a millimetre off their axes,
  to 5·10⁵ m (meet refuses conics farther out); where it stops the sags
  are under a sixteenth of the tolerance near the origin (within 10² m;
  the hyperbolas to 10⁴ m, the Bézier pieces at the origin, up to 4 m
  across), and a conic track's wherever its floor is under it by a
  thousandth (parabolas 260 km to 880 km from their vertices, one just
  past that held over it), and so for every reach meet accepts, swept
  across its limit; the exact courses of conic tracks in planes'
  and cones' charts (the anchor up to 10⁴ m from the vertex, stretches
  1 m to 10⁷ m from it, half ellipses to their ends, where `w` is under
  10⁻⁸) lie within their bands, as do those of NURBS arcs in cylinders'
  and cones' charts whose point lies up to 10⁶ m from them along the
  axis, and the hull distances of NURBS curves in planes' charts
  (degrees two to five, weights' ratios to 10⁴, 1 m to 10¹⁰ m from the
  origin, stretches down to 10⁻¹² of a piece) round within their
  allowance, against double-double references (`chart_rounding.rs`);
  each of these fails without its allowance, its floor, its weights'
  term, the chart's part of the slack of a track or of a NURBS course,
  or a conic's halving to its floor (ten mutations checked); points on
  the ellipse a cut leaves on a rod 1 m to
  3,000 km from the origin are located at once; the probes by turning,
  through the stages' own entry, sit at their fractions along stretches
  past a half turn of an ellipse, of a NURBS ellipse in quarters, and of
  one cubic piece holding an arc of 238° of an ellipse.
* Golden second cuts, all four operations with exact volumes, counts
  and origins: a conic-edged body against a box far away; a rod cut by
  a tilted plane and then by a second plane crossing the first ellipse
  (three surfaces at a point); a pointed cone under a ridge; the tilted
  countersink plate halved by its mirror plane and cut by a second
  countersink beside the first; a flattened drafted pin halved through
  the hyperbola's vertex and flattened again on its other side; re-cuts
  in the same plane and with the same cone (a stored edge and the track
  where faces meet coincide, §16.3); a cut that leaves a void beside
  an ellipse; a rod cut twice at the origin and from 1 km to 3,000 km
  from it, its volume within a millionth where it is cut (to 30 km),
  refused as too far beyond, and soon either way; a cone's hyperbola a
  millimetre or a centimetre off its axis cut again from the origin to
  100 km from it, within 2·10⁻⁹ in volume of the cut at the origin to
  30 km; and, first cuts, cones 1 m to 400 km tall cut by a plane
  parallel to their axis, within 10⁻¹² in volume over size cubed, and
  one 600 km tall refused (its hyperbola's vertex too far away, §14);
  needle cones, their hyperbolas' vertices at the origin, cut to
  877.8 km long and refused by meet from 879 km, never as too far from
  the origin.
* Fuzz: second cuts of random cut cones and cut rods, and chains of
  three or four cuts, with the volume identities and point membership,
  at 50 times their size, asserting that no refusal names a free-form
  edge and counting the conic edges exercised.
* Time guards: sixteen tilted countersinks and sixteen flattened pins
  cut one at a time, each against the same cuts in one boolean, with a
  scaling check.


## 17. As built: operands with conic edges on curved faces (2026-10-15)

Bodies whose conic edges lie on faces of cylinders and cones are
operands of every boolean now: the stages follow the stored rational
curve in the face's chart, with the band of §16.2 (ADR 0023). A
countersunk hole leaning 10° can be cut beside another, a drafted pin
flattened on both sides, a cut repeated in the plane or with the cone of
an earlier one. §16 describes the rules as built, after review corrected
them (the rounding allowances, the probes by turning, the foot near the
end of half an ellipse, halving by the sag and the refusal where
rounding holds it up); each correction has a test that fails without
it. Besides them:

* **A line that misses where faces of both bodies can be.** The fuzz of
  second cuts found a fault in a first cut: a slab whose face passes
  through a frustum's apex within 10⁻⁷ rad of tangent along the ruling
  of its seam. The plane meets the cone in two rulings 0.4 mm apart and
  stays within the tolerance of it between them out to about a metre
  from the apex, just short of the frustum's top. There the slab's face
  ends, its boundary grazing the cone, and the cut where one ruling
  leaves that face was taken for a touch; the ruling's cuts in the face
  did not pair up, and the boolean failed as meeting the other body's
  surface inconsistently, an internal error where a refusal belonged.
  Neither ruling reaches the box where faces of both bodies can be
  (§15), so no part of a result lies on them: a line or a ray that
  misses the box now has no stretches in a face, and all four
  operations on the pair are refused as faces meeting at too shallow an
  angle (`a_plane_by_a_frustums_apex_along_its_seam`).
* **Two edges on one conic.** Where a stored edge's course and a
  track's course of the same conic (a cut repeated in one plane) meet
  end to end, the merge of edges does not join them, as it joins pieces
  of one edge or of one track: the result keeps a vertex between two
  edges on one curve (review). It is valid and both edges are exact;
  joining them would need a curve made from pieces of two, which nothing
  needs yet.
* **What the validator samples.** It checks each edge against its faces
  at five points evenly spaced in its parameter, the blind spot that
  matching by parameter had (§16.3), and was left so: on a curved face
  it allows 100 tolerances (`check::on_tolerance`), while two courses
  the matching rule takes as one stay within 2.3 times the farthest
  probe of each other in review's adversarial cases, so samples by
  turning would find nothing over its limit there. What keeps different
  courses apart is the matching rule, whose test fails without the
  probes by turning.
* **Fuzz**, at 50 times, on the final code (2026-10-15). Rods and cones
  cut by a random slab, each cut again by a random box, slab or rod,
  with the volume identities and point membership: of 6,255 cut rods
  (3,451 with conic edges, 11,215 of them), all 25,020 operations done,
  valid and correct; of 4,671 cut cones (3,787 with conic edges,
  14,790 of them), 2,564 of 18,684 operations refused (137 in 1,000):
  2,092 where a plane passes through or near a pointed cone's apex, 400
  tangent contacts, and 72 the tolerance cannot decide (25 slivers, 23
  faces within the tolerance of the other surface, 12 shallow angles,
  12 points within the tolerance of an edge). Chains of three or four
  cuts (3,179 of 4,415 first cuts leaving conic edges, 13,763 of them):
  540 of 14,356 operations refused (38 in 1,000), 456 at an apex, 60
  tangent contacts, 24 others. No refusal names a free-form edge or the
  distance from the origin. The fuzzes main ran refuse what they
  refused on main, operation by operation, but for two unions of
  solids at an angle that main refused as having a free-form edge on a
  curved face, done now and correct by the identities and point
  membership; the flat fuzz at 100 times refuses exactly what it did.
* **Speed** (docs/PERFORMANCE.md, sixth pass). The check that no edge
  of a turned face runs through its axis sampled every edge at 65
  points, with a search at each local least, many of them along a line
  parallel to the axis; lines and circles now take closed forms and
  NURBS curves the convex hulls of their pieces first. Validating the
  tilted countersunk plate fell from 14.3 ms to 4.2 ms, the median of
  12 runs interleaved with main's, and cutting its holes in one boolean
  from 41.1 ms to 24.3 ms. The band of a NURBS course in a plane's chart
  uses its trace's pieces instead of splitting the curve again for every
  chord. A cut one at a time into conic-edged cones takes about twice as
  long as one into the straight countersinks (16 take 197 ms against
  102 ms), most of it in the bands of their courses (a known slow spot
  there).
* **The second review** (of the fixes and the speed work above) had six
  findings confirmed by at least two of three skeptics: five fixed, each
  with a test that fails without the fix, and a false claim in the
  code's documentation corrected:
  * the validator's quick test of a circle against its face's axis took
    the least singular value of its frame as the square root of a
    difference of eigenvalues, which loses about `√ε` where the circle's
    plane holds the axis: a ball's meridian of 700 m, leaning, was
    cleared at the default tolerance though it runs through both poles,
    and of circles of 10 m passing the axis a hair inside the tolerance,
    808 of 3,000. It is now `|X × Y| / σ₊`
    (`circles_along_a_leaning_axis_are_not_cleared_through_it`);
  * three on the floors of sags (§16.2): point location and the
    refinement halving without end beyond about 145 km, and the
    search's rule that could refuse a refinement a round from done
    (`sags_narrow_to_their_floors`,
    `points_on_a_far_conic_edge_are_located_at_once`,
    `second_cuts_far_from_the_origin_finish`);
  * the probes by turning crowding into the first part of a stretch
    that turns past a half turn (§16.3,
    `probes_by_turning_follow_stretches_past_a_half_turn`);
  * a claim in the code's documentation that every halving quarters a
    sag, false for a NURBS course on a plane (§16.2; a second copy of
    it, on the refinement's rounds, survived until the third review).

  Three more were not confirmed, and were acted on all the same: a
  broken link in the code's documentation; a conic track's rounding
  allowance carried in its sag, now that sag's floor; and an edge whose
  ends both lie within the tolerance of the axis, reported through it
  though sampling finds its least at an end. Such an edge runs along
  the axis within the tolerance all its length, and the closed form for
  lines now says so (`line_through_axis`).
* **The third review** (of everything above) had six findings confirmed
  by at least two of three skeptics:
  * the probes by turning broke a NURBS course only at its knots, but a
    piece of degree three or more can turn past a half turn within a
    span (an arc of 238° of an ellipse as one cubic piece), and the
    probes crowded as before; a piece is now halved until its control
    points' differences show it turns less (§16.3);
  * the test of the probes by turning called the helpers itself for a
    NURBS course, so it could not fail when the stages' own entry
    changed; it now goes through that entry
    (`probes_by_turning_follow_stretches_past_a_half_turn`);
  * a conic track's allowance for rounding, `10⁻¹³` of its reach,
    refused first cuts main did: a cone 20 km tall cut by a plane
    parallel to its axis, though both operands sat at the origin. The
    allowance is now 32 ε of the reach, set by measurement against
    double-double references (§16.2), and STATUS no longer says such
    refusals need a conic edge on an operand
    (`first_cuts_of_large_cones_finish`);
  * the halving test could not fail without the conic track's floor or
    the NURBS courses' on planes, only the floor of NURBS courses on
    cylinders and cones; it now covers far conic stretches, the
    hyperbolas of cones cut off their axes and large unevenly weighted
    pieces, and fails without either floor
    (`sags_narrow_to_their_floors`);
  * STATUS said a rod was cut twice up to 3,000 km from the origin; it is
    cut to 30 km and refused at 300 km, as its test allows;
  * the documentation of the refinement's rounds still said each round
    quarters the bands it halves, and the second review's entry above
    said all six of its findings were fixed with tests; both corrected.

  One finding fell short of a majority of the skeptics (the refusal it
  found was new, not a regression from main) and was acted on: the floor of a
  NURBS course on a plane, in proportion to its weights' ratio times its
  whole size, refused second cuts across a cone's hyperbola a
  millimetre off its axis from 1 km from the origin, where halving
  further resolves them exactly. Measured against double-double
  references, the weights' ratio enters only through the rounding of the
  clip's parameters, in proportion to the curve's own size; the floor
  is now the derived bound of §16.2, carried in the sag's bound as its
  allowance, and such cuts are done to 30 km
  (`second_cuts_of_uneven_hyperbolas_far_from_the_origin`). The same
  measurements found the positional allowance missing the chart's own
  size: a conic's chart points strayed up to 26 times it with the
  chart's anchor 10⁴ m from the conic's vertex; it now includes it
  (§16.2). The tests of both measurements stay
  (`conic_courses_lie_within_their_bands`,
  `hull_distances_round_within_their_allowance`), and each of the seven
  ways of dropping an allowance, a floor or a term of either fails one
  of them or the halving test.
* **The fourth review** (of the retuned allowances and their tests) had
  six findings confirmed by at least two of three skeptics:
  * halving stopped at four floors, though a conic track's bound settles
    on its floor itself: with the floor between a quarter of the
    sixteenth and the sixteenth (220 km to 880 km of reach) a chord
    could stop over the sixteenth and refuse a boolean that more halving
    resolves, a cone 400 km tall cut by a plane parallel to its axis
    among them. A sag now says whether it settles on its floor, and
    such a chord is halved to within a thousandth of the floor (§16.2;
    `conic_sags_settle_on_their_floors`, and the cone 400 km tall is
    now cut in `first_cuts_of_large_cones_finish`);
  * that test required the cone 300 km tall to be cut, which it was by
    0.5% (its end chord stopped at 0.995 of the sixteenth), and its
    documentation put the hyperbola's reach at 200 km where it is 281 km;
    it now holds a margin, and asserts the refusal of a cone 600 km
    tall by meet (two findings);
  * no test failed without the chart's part of a NURBS course's slack:
    the course in a cylinder's or a cone's chart whose point lies far
    along the axis is now checked against double-double references
    (`developed_courses_lie_within_their_bands`);
  * the conic band test's points near a half ellipse's end were never
    checked below `w` of about 3·10⁻⁷, its stretch there being under the
    parameter's last bit; they are now taken by the last bits, `w`
    under 10⁻⁸ (3·10⁻⁹ in the default debug run, 9·10⁻¹⁰ in the
    release one);
  * §16.5 and the halving test's documentation claimed the off-axis
    hyperbolas to 10⁸ m from the origin (they are tested to 5·10⁵ m) and
    large Bézier pieces down to the sixteenth (none was checked); the
    claims are corrected and the pieces at the origin now checked.

  Each of ten ways of dropping an allowance, a floor, a term of either,
  or a conic's halving to its floor now fails one of the rounding
  tests or the halving tests, and without that halving the golden test
  of large cones fails, its cone 400 km tall refused.

* **The fifth review** (2026-10-15), of the settling stop, its tests and
  the corrected documentation, had six findings, every one checked
  against the code:
  * meet accepted a conic whose stretches' floors came within a
    thousandth of the sixteenth, where the settling stop can end over
    it: needle cones 879 km to 879.6 km long, their hyperbolas'
    vertices at the origin, were refused as too far from the origin,
    though §16.2 and the ADR said meet refused first. meet now refuses
    with two thousandths to spare (§14; `meet::too_far`); the halving
    test sweeps every reach meet accepts across its limit, and
    `needle_cones_at_the_far_limit` cuts them to 877.8 km and has meet
    refuse them from 879 km, both failing without the margin;
  * four were false claims in documentation: the halving test's Bézier
    pieces are up to 4 m across, not 2 m; the conic band test's least
    `w` was claimed to 10⁻⁹ where it asserted 10⁻⁷ (it now asserts
    10⁻⁸, the default runs reaching 3·10⁻⁹ and 9·10⁻¹⁰); a scaling
    check compares four times the holes counting the one cut, not the
    holes before it; and a test helper's cone was described at twice its
    size. The sixth, the record's unfilled test count, was filled in
    with the commit.
