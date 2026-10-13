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
mesh with an error bound (twice the mesh's area times its distance from
the surfaces), or by the exact mass properties where that bound could
reach the sign, and each void goes in the smallest solid holding one of
its own corners.

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
  §2.4). The centre is the middle of the face's boundary in its
  surface's own parameters (longitude within `[0, π]` of the surface's
  frame, ADR 0008, and latitude), so the antipode is outside the face
  and the projection is one to one on it; the poles are ordinary points
  of the chart.
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
  normal. The way a track runs across a face, and the side of a surface
  a face's boundary goes to, are decided at cuts where both surfaces
  have normals; a face's facing is found at a regular point of it, not
  at its first corner.
* **Classify.** A ray meets a cone's nappe and a sphere at the roots of
  a quadratic, anchored at the cone's reference point (a root on the
  other nappe is not a hit); a hit where the ray is within 1e-3 of
  tangent to the surface abandons the ray, as on a cylinder. The point
  of a cone or a sphere nearest a point, for the test that a piece's
  point is clear of the other body, is in closed form.
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
