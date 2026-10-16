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
bodies are their own later decision (an ADR when it comes up). Bodies
that touch without crossing are decided by §18 and ADR 0024: no trace
of the contact, or a refusal as non-manifold. As built (§19) this covers
edges and corners resting on flat faces, the faces there all flat;
contacts where a face is curved are still refused as tangent contacts.

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
non-manifold results (§18 for bodies that touch); faces or edges
thinner than the tolerance left by near-coincident geometry (the same
choice the sketch arrangement makes for faces it cannot absorb as
artefacts of the tolerance).

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

## 18. Contacts that do not cross (design, 2026-10-16)

Two bodies can meet without crossing: a corner or an edge of one resting
on a face of the other, a rim on a plate, the root of a groove on a face,
a ball on a plate. The booleans have had no rule for this, and they
answer in three ways depending on the faces involved:

* **Curved faces refuse.** An edge touching a face's surface at a point
  of the face is refused by `refuse_touches`, and surfaces that touch
  are refused by `relation`, both as "tangent contacts". At fifty times
  its default size the curved fuzz counts 248 operations refused as
  edges touching faces and 3,152 as faces touching faces.
* **Flat faces return bodies that touch themselves.** `plane_meeting`
  skips `refuse_touches` when both faces are flat, and drops the
  interval of one point that a vertex resting on a plane leaves. The
  planar split then drops an edge lying inside a face with a free end
  as an antenna (`prune`). A plate's union with a cube standing on one
  corner, or with a wedge resting on its edge, inside the plate's top
  face, comes back as two solids that touch at the corner or along the
  edge, and the result validates. If the plate also has a round hole
  elsewhere, its top face is no longer `flat` and the same contact is
  refused as unsupported: the answer depends on an unrelated hole, and
  `Kernel::boolean` says such results are refused (found by review,
  confirmed on main on 2026-10-16).
* **A slit.** An edge lying in a face from its outer loop to a hole's
  loop is not an antenna, and `prune` keeps it: the face's loop runs
  both ways along it, joining the hole to the outer loop. A plate with a
  square hole, minus a knife lying from the plate's edge to the hole,
  comes back with (10, 27, 18) faces, edges and vertices instead of the
  holed plate's own (10, 24, 16), and validates (found by review,
  confirmed on main).

This section gives one rule for every such contact, on flat and curved
faces alike. It was chosen from three candidate designs, each reviewed
against the code (ADR 0024): detection at the imprint, from runs of a
loop's marks on the other surface; detection in the face's planar
split; and a local model of the contact's neighbourhood, decided before
the split. The second is the base, with ideas from the other two
grafted on where review showed them needed. The work comes in three
slices: an edge or a vertex of one body against a face of the other
(slice 1, built first); surfaces that touch (slice 2, designed here);
and the tangencies refused today for other reasons (slice 3, outlined).

### 18.1 What the result is

The booleans compute regularised sets, the closure of the interior of
the set operation's result (Requicha 1980), and the regularised
operations are closed on closed regular sets (Tilove & Requicha 1980),
so every boolean of touching operands has one exact result.
Regularisation drops the lower-dimensional parts: the line or point
where two bodies touch from outside. It does not remove a **pinch**,
where two lobes of the result meet at a line or a point only; the
result is then a regular set that is not a 2-manifold. Phase 3 refuses
non-manifold results (§3, ADR 0020), and this design keeps that.

**The four cases.** Let an edge or a vertex of body T touch face F of
body P at a contact x, and let

* σ = +1 when T's boundary next to x lies on the outward side of F, −1
  when on its inward side;
* τ = +1 when F next to x lies outside T (T is convex at x), −1 when it
  lies inside T (x is the root of a groove or a notch of T).

Every point of F near x has one class against T, and every piece of T at
x one class against P; these decide the selection (`kept`). "P" or "T"
below means that operand's surface alone with no trace of the contact;
"none" means no surface near x (if nothing else overlaps, "the result
would be empty").

| case | σ | τ | P ∪ T | P ∩ T | P − T | T − P |
|---|---|---|---|---|---|---|
| (i) external | + | + | pinch | none | P | T |
| (ii) nested | − | + | P | T | pinch | none |
| (iii) on a groove of T | + | − | T | P | none | pinch |
| (iv) a notch of T on F | − | − | none | pinch | T reversed | P reversed |

Each case pinches in exactly one operation: the union of a cube on its
corner on a plate; a void touching the outer wall; a cavity touching a
groove's root; two lobes of an intersection meeting at a notch's root. A
pinch fails with `InvalidInput` naming the contact. Volumes still obey
the identities of §7: the contact has no volume.

**No trace.** A result that leaves no trace of a contact equals, entity
by entity, the same boolean with the touching body moved clear along the
contact's normal, away from F, far enough to remove the contact and no
other relation: the same face, edge and vertex counts, the same face
origins in order, each face's loops as cyclic sequences of edge origins
up to rotation, the same vertex origins and merged origins, the same
volume and area. The reference is the moved boolean, not the operand:
the stages normalise some operands even when the bodies are apart (a
turned knife rim of (6, 10, 6) comes back as (5, 8, 5) from both the
touching and the separated difference, measured on main). For operands
bounded by planes and straight edges the two references agree, and the
tests assert both.

**Where the contact sits.**

* **Inside F.** The other surface stays on one side, so the contact
  separates nothing in F. A lying edge is a *bridge* of F's planar
  graph (the same region on both sides of it), a touching point an
  isolated node; neither ever enters a loop of a result face. The rule
  of §18.2 decides them.
* **On F's boundary.** The contact is a node of F's graph, or a segment
  F shares with the touching edge, and the existing tests decide it:
  `check_closed` refuses an edge used more than twice, and
  `check_vertices` a vertex whose faces form more than one fan.
  Together they are the neighbourhood test of a contact edge or vertex
  (Mäntylä 1986) in counting form. At a vertex of F the two bodies'
  vertices share a cluster and `check_vertices` decides.
* **Closed.** A rim lying on F with its faces on one side is a closed
  chain that separates F. `check_closed` decides it (four uses is a
  pinch), and in a result with no trace `merge_faces` joins F's two
  pieces again: main already returns the plate minus a knife rim
  resting on it unchanged.

**What is not such a contact.** A *crossing tangency* is a point where
the surfaces touch and also cross: their relative second fundamental
form is indefinite there, and the point is singular on the curve where
they meet (Patrikalakis & Maekawa 2002, §6.4.1). A curve where faces
meet that is tangent to an edge at the contact is one. The result can
be manifold, with a vertex of four branches; slice 3 handles these, and
until then they are `Unsupported`. Coincident faces are not contacts:
the coplanar path of §9 and §16.3 handles them.

**Touching bodies are never returned.** A union whose solids would touch
is refused, not returned as a body of several solids that touch: such a
body touches itself, and later booleans assume their operands do not.
Disjoint operands still unite into bodies of several solids.

**Error kinds** (ADR 0020): `InvalidInput` where the tolerance makes the
result non-manifold or undecidable (pinches, near-tangencies over a
stretch, nearly coincident surfaces, cusps, points within the tolerance
of edges they are not on, contacts too shallow to classify);
`Unsupported` for what is not built yet (crossing tangencies, and
surfaces that touch until slice 2); `InvalidResult` where stages
disagree. No flat input can reach the crossing-tangency refusal: a
straight edge meets a plane at one point, lies in it or misses it, so on
flat faces a contact that does not cross is always at a vertex, and the
flat fuzz keeps counting any `Unsupported` as a bug.

### 18.2 Slice 1: an edge or a vertex against a face

A contact reaches the split of the face it touches in one of three
forms: an edge lying on the face, imprinted as today; a point touching
it, recorded by the new `touches_in` as a node with no segments; or the
end of a lying edge. `planar::regions` reports what separates nothing,
bridges and isolated nodes, and each becomes a *place*, held by the
piece whose region contains it. After selection one rule decides every
place: **it pinches if and only if kept pieces of both operands are at
it** (`assemble::check_contacts`). Contacts on a face's boundary are
shared nodes or segments, which `check_closed` and `check_vertices`
decide.

**Lying edges stay imprinted.** An edge of T lying on F is imprinted by
`plane_meeting` or `stretches` as the meeting of F with T's faces at the
edge, and `Segments::add` and `add_curved` merge it with the edge's own
segment. It stays so: detection rests on what is imprinted, not on which
runs along a loop a test confirms as touches. A design that filtered a
confirmed lying edge out of the imprint was found by review to return a
result main refuses: a block crossing a plate's top face F, with a
through V-groove whose root lies inside F, the block's front and back
faces crossing F through the root's two ends (reflex corners of those
faces). With the root filtered out, the root's ends were left to
topology, which saw one fan at each, and A − B came back as two prisms
touching along the root. Kept in the imprint, the root is a chord
separating F's footprint, and `check_closed` refuses it with four uses.

**Whole lying loops.** A face's loop can lie wholly on the other surface:
a plate's hole rim resting on a ball, the plate's bottom face holding the
rim as an inner loop. The curve where the bottom face meets the sphere
is that circle, so the track meets no edge: `spans_on`'s branch for a
track with no ends finds every probe on the face's boundary and returns
the whole track, the whole-track path cuts it at two new seam points,
and the seams lie on the rim edge but are none of its vertices, so the
split refuses them as points within the tolerance of an edge (found by
review). The branch changes: when every probe is on the boundary, it
looks for a loop of the face whose edges all lie on the other surface
(all their stations zero) and returns that loop, whose own edge segments
the caller imprints in the other face, each gaining the `Meet` source,
with no seams. With no such loop the track runs within the tolerance of
the boundary without lying on it, refused with the existing text "a
point where the bodies meet lies within the tolerance of an edge it is
not on". Imprinting nothing for such a pair, the other choice reviewed,
would lose the contact from this pair's faces and leave it to other
pairs, of which none need exist.

**Points.** A vertex or a station of T's boundary touching F inside F
leaves no segment, so it is recorded separately:

* *Zero flags.* `relation` already computes, through `sides`, whether
  each face's boundary has a mark within the tolerance of the other's
  surface. `Relation::Crossing` carries the two flags, and only the
  flagged sides look for touches, so bodies with no contact pay nothing.
* *Runs.* `touches` reads each loop's marks on the other surface across
  vertices and reports a run of zero marks between neighbours strictly
  on one side. It skips runs holding an edge that lies on the surface,
  runs through a vertex that another edge leaves to the other side
  (`leaves_to`), and runs through a vertex inside a smooth piece of
  surface (`smooth_at`). It now also reports each run's side: the sign
  of the neighbour before it, times +1 or −1 as the touched face's normal
  is outward or not. Every mark of a run that is not inside a lying edge
  must fall in one cluster, or, in a run with lying edges, in the
  clusters of those edges' ends. Anything else is an edge running within
  the tolerance of the surface over a stretch, which the tolerance
  cannot place (an arc leaving a lying edge's end tangentially, its
  first stations within the tolerance beyond that end's cluster),
  refused with `InvalidInput`: "an edge of one body runs within the
  tolerance of a face of the other over a stretch without lying on it".
  On flat faces it cannot happen: two vertices of a straight edge
  within the tolerance of a plane make the edge lie on it.
* *`touches_in`* replaces `refuse_touches`: the same loop over both faces
  of a pair, recording for the touched face its operand, the face, the
  cluster, the touching mark (a vertex or a station of an edge) and the
  side (+1 when the touching boundary lies outside the touched body).
  Faces with an outline keep today's filter, the point in the closed
  face; flat faces skip it (it scans every edge of the face), because
  the split places the point exactly. `plane_meeting` calls it for every
  flagged pair, and the guard that skipped flat pairs is deleted: one
  path for flat and curved faces, the fix of the flat fault.
  `stretches` calls it instead of `refuse_touches`. `plane_meeting`'s
  filter of an interval of one cluster stays a filter: two faces'
  intervals meeting at one cluster are also an ordinary contact of two
  edges.

**The cut-pair filter.** On a curved pair, a touch at a vertex or a
station leaves two cuts at one cluster with opposite `up` under the
half-open rule that sees the zero mark, the cuts into and out of the
surface. `refuse_touches` used to stop such pairs before `assign` saw
them; without it, `assign` would refuse those cuts as "too shallow an
angle" or "nearly touch" (found by review). `cuts` therefore cancels,
per cluster, a pair of cuts with opposite `up`. A loop visits a cluster
once (`settle` refuses two vertices of one body in a cluster, and
`incidences` a point on two edges of one body), so a cluster carries
one cut, a crossing, or one cancelling pair; a debug assertion checks
it. A zero run longer than one cluster puts its two cuts in different
clusters, which the filter leaves alone. The cancelled pair adds nothing
to parity or to any span (`spans_on` already drops clusters whose cuts
net to zero). Where the face itself crosses the surface through that
cluster (a reflex corner, a ray through a cone's apex), removing the
pair also removes `assign`'s checks for it, so the near-apex goldens are
re-run and the fuzz outcomes diffed by class at that commit. For two
planes the filter changes nothing: `face_line` merges the two rules'
intervals, and an interval of one point is dropped.

**The Apart band.** `meet` can call two surfaces apart while a mark of
one face's boundary lies within the tolerance of the other: it snaps a
nearly parallel pair to parallel when the tilt times the reach is at
most ε/16; and two planes are apart when g's boundary lies strictly on
one side of f's plane, even if a vertex of f lies within the tolerance
of g's plane. With the vertex strictly inside g the second cannot
happen, since g's boundary then reaches both sides of any plane through
the vertex; it can where the vertex lies within the tolerance over the
sine of the planes' angle of g's boundary. A pair decided apart with a
zero flag set would still have to look for touches. How often it occurs
is measured before such a path is built; if the fuzz never reaches it,
it becomes an `InvalidResult` check instead.

**Settling a touched point.** After the imprints and before the split,
the touched points are sorted by (operand, face, cluster) and
deduplicated; one point touched from both sides is `InvalidResult`.
Then, for each touched cluster c in face g:

1. *A station of an edge e of T, where a segment of g carrying the
   meeting of a face of e ends.* The face crosses g's surface along a
   curve through c, tangent to e there: a crossing tangency, refused as
   `Unsupported` whatever else ends at c ("booleans where an edge of one
   body touches a curve where the faces meet (a crossing tangent to the
   edge; design note 0012 §18)"), since the planar split's exact tie of
   directions need not see it.
2. *A vertex v of T, where a segment of g carrying the meeting of a face
   at v ends.* An ordinary crossing through a vertex: a reflex corner of
   a face crosses the surface while all its edges stay on one side (an
   L-prism's inner corner on a slab tilted to (1, 1, 1)/√3, which main
   returns correctly in all four operations, measured). The point is
   dropped; the meeting is split at v in g, v is a corner of every face
   at it, and `check_vertices` decides. `leaves_to` cannot mark these
   crossed, because it reads edges, not face sectors.
3. *Any other point that ends a segment of g* lies on g's boundary or on
   a contact segment and is dropped: that segment's place, or
   `check_vertices`, decides it.
4. *Anything else* is an isolated point of g, passed to the split with
   its side.

**The split.** `planar::regions` takes the face's isolated clusters and
returns its regions and its *loose* items, each a bridge (by its
segment) or an isolated point, with the region that holds it:

* *Isolated nodes* are appended to the clusters, so the refinement and
  `consistent` test them against every chord: a point within its reach
  of a chord it is not on is the existing `Problem::Close`.
* *Bridges.* The walk records each half-edge's cycle. A chord is a bridge
  if and only if both its half-edges lie in one cycle: in a plane graph
  an edge has the same face on both sides exactly when its removal
  disconnects its component. This covers antennas, isolated segments and
  trees, and the connector between two loops that `prune` misses. A
  bridge with a side (an edge of the face's own boundary) is
  `Inconsistent`, "an edge of the face bounds no region". The bridges'
  half-edges are removed from the angle-sorted lists and the walk is
  run again; `prune` stays only as a debug assertion that no half-edge
  is followed by its twin. Components are joined over chords that are
  not bridges; cells, outsides and holes are computed as before.
* *Placing a bridge.* A half-edge of its first-walk cycle that is not a
  bridge lies, in the second walk, on the outer cycle of a cell or on
  one of its holes, and the bridge lies in that cell. A bridge whose
  first-walk cycle is made only of bridges (a tree) is placed by one of
  its clusters, as an isolated point is: in the innermost cell whose
  outer cycle contains it (by the bounding-volume tree and `contains`
  that place holes). A point in no cell inside the face, or in a hole of
  its cell, is off the face and dropped: a corner resting on the plane
  in the mouth of a through-hole touches nothing (found by review). A
  bridge off the face is `Inconsistent`, "an imprint edge lies outside
  the face". Containment is decided robustly, since `consistent` has
  put every such point farther than its reach from every chord.
* *Order.* Loose items are listed bridges first, by segment, then
  points, by cluster.
* *Bit-identity.* With no bridges and no isolated points the first walk
  is the only walk and the output is unchanged; a unit test asserts it.
  `merge_faces` passes no points, and every edge it passes has a side,
  so a bridge there is the existing `InvalidResult`.

**Places.** A bridge gives a place on its segment, plus one at each end
that is not a node of its region's loops (a free end). Its segment must
carry an edge of the other operand (slice 2 adds the touch source), or
it is `InvalidResult`: "a contact edge carries no edge of the other
body". The split used to drop a crossing curve that dangles (a numerical
miss) silently; this check now raises it, the fuzz is triaged before the
commit lands, and a legitimate class of such bridges, if one appears,
gets its own design. An isolated point gives a place at its cluster. A
face clear of the other body has no imprint, so a loose item in one is
`InvalidResult`. Each place records σ: for a point, from its run; for a
bridge, from the directions r₁ and r₂ into T's two faces across the
edge's middle, as the sign of n̂_g · rᵢ, the two agreeing or
`InvalidResult`.

**The rule** (`check_contacts`, before the pieces are consumed). With
`keep[i]` whether piece i is kept, each place collects, per operand,
whether a kept piece of it is at the place: its holders (the pieces
whose regions hold it), and the *incident* pieces of the other operand,
those whose loops have the place's segment or a corner at its cluster,
found in one pass over every piece's loops. The place pinches if and
only if both operands keep a piece there. The least pinching place is
reported, a segment before a point and then by index, with the texts of
`check_closed` and `check_vertices`: "the result would be non-manifold:
the bodies touch along an edge near {p}", "… at a point near {p}", at
the segment's middle or the point. Reports go through ordered maps, so
no result depends on hash order.

**Why the rule is sound and complete.** It fires only when a kept face
of one operand passes through x (the holder, x inside it) and a kept
face of the other reaches x: two sheets of different operands meet at x
without crossing, a pinch. When only one operand keeps material at x the
result is manifold there: a disc for the holder alone, or the operand's
own surface near its edge or vertex for the incident pieces alone. In
each of the cases (i) to (iv) the rule refuses exactly the operation the
table refuses. Every point of the selected surface lies inside a kept
piece away from contacts, on a segment of a kept loop (`check_closed`),
at a corner (`check_vertices`) or at a place (this rule): the
arrangement holds every imprinted point, loops holding segments and
corners and the loose items holding bridges, isolated points and free
ends. So completeness reduces to what reaches the split: `touches_in`
on every pair with a zero flag, lying edges imprinted (whole lying loops
along their own segments), the stretch refusal and the Apart band.
Surfaces that touch are refused at `relation` until slice 2. An
independent oracle in the tests (§18.5) guards this dependency.

**Classifying a holder.** A holder is classified at one test point
(`interior`, `Boundary::classify`), which decides by exact ray parity
and refuses within the tolerance of the other surface. The holder meets
T only along its loops and at its places, and those do not cross, so its
class is uniform and any point far enough from the places decides it.
Nothing kept the test point away from a contact: a convex piece takes a
central point, and a wedge's edge through the centroid of a plate's top
face puts the point on the contact. Places in the region become
*obstacles*: a straight bridge its segment, a curved bridge its chords
with their band, a point a point. `classify::interior_point_avoiding`
triangulates the region's polygons with the obstacles' points as free
vertices (sampled at most a 32nd of the piece's size apart; not
constraints, so the inside triangles are unchanged) and scores each
inside triangle's incircle, centre q and radius r, by min(r, d(q) −
band), d the exact distance to the obstacles. For large pieces
`probe_point` takes each obstacle as two extra sides, one each way, so
the winding number is unchanged and the clearance includes them. The
sliver refusal keeps the clearance from the loops alone, so a crowded
holder is never called a sliver; where the loops leave room but the
obstacles leave none beyond the tolerance, the holder is refused: "a
contact of the bodies lies within the tolerance of a narrow piece of a
face". At distance d in F from the contact, T's faces, crossing F's
tangent plane at an angle θ, are at least about d sin θ away, so the
test point needs d beyond about ε / sin θ: if `classify` refuses the
best point as near, the next candidates by score are tried, up to four,
and then the holder is refused: "the bodies touch at an angle too
shallow for the tolerance". A point too near can only cause a refusal,
never a wrong class.

**Cross-checks** (each `InvalidResult`), a second derivation of each
class that does not use the first:

* a place with no piece of the other operand at it;
* incident pieces whose classes differ or are on a coincident face, or
  a holder on a coincident face (a contact inside a coplanar overlap
  would need a vertex of T inside T's own face);
* an incident class other than `Outside` for σ = +1 or `Inside` for
  σ = −1;
* *the probe across*: q = x − 4ε σ n̂_g lies on the side of F away from
  T's touching boundary, at least 3ε from it on surfaces of radius over
  8ε, and classified against T must give the holder's class (skipped
  where q itself is within the tolerance of T);
* *the wedge*, for a bridge: with n₁ the outward normal of one touching
  face at the edge's middle and r₂ the direction into the other, the
  holder is `Outside` exactly when the edge is convex (n₁ · r₂ < 0),
  whichever side T's faces rise to; skipped where |n₁ · r₂| is below the
  angular tolerance.

**Origins** (ADR 0024, amending ADR 0020). Contacts inside a face leave
the holder's loops as they were, but a lying edge's segment carries the
meetings of T's faces with F, and an edge lying on an edge of the other
operand carries both operands' edges. `build` listed every source, so a
result with no trace would not equal the apart case (found by review):
W − P kept the knife edge with the meetings as merged origins, and A ∪ B
= B for a body A nested against an edge of B, A first, made A's edge the
primary origin of B's edge, the first operand's sources sorting first.
The rule is: **an origin names only operand entities that bound the
result entity.** In `build`, for an edge on segment s, let U be the
union of the sources of the result faces along s (merged faces and
coplanar partners included). `Edge(op, e)` is kept if and only if some
(op, f) in U has e in its loops, from a table of each operand's edges'
faces built once, and `Meet(fa, fb)` if and only if (First, fa) and
(Second, fb) are both in U; the least kept source is the primary, the
rest merged, and an edge left with none is `InvalidResult`. A vertex's
operand vertex is kept if and only if some result face at it comes from
a face of that operand holding the vertex, and a vertex keeping neither
is a meeting. Face origins are unchanged. Differences are confined to
results where an operand keeps no face at a shared entity, and every
golden that changes is reviewed at that commit. `origins_agree` checks
only that each listed origin lies on its entity, so the filter cannot
break it, and the naming crate does not yet read boolean origins, so no
stored reference changes.

**Splits rejoined.** A contact cluster on an operand's edge splits it,
and `merge_edges` joins the pieces again in a result with no trace. It
declines at vertices on a turned face's axis, so a box's corner on a
revolve cap's centre would leave an extra vertex: valid, but not equal
to the apart case. The remedy proposed is to allow the join when both
segments are pieces of one operand edge, rebuilding that edge's own
curve; a golden decides whether it is needed.

**Dropped contacts are not recorded.** A query for where two bodies touch
is a separate feature (ADR 0024). Design note 0011 §4.4's "isolated
tangent point returned as a vertex" is superseded: a contact inside a
face is never imprinted into a loop, following the precedent of not
subdividing polygons that only touch (Laidlaw, Trumbore & Hughes 1986).

**Refusals added by slice 1** (each followed by "near {p}", the contact's
own location: the least place or the cluster's point, never a face's
first corner, which is what `relation` reported):

| text | kind | where |
|---|---|---|
| an edge of one body runs within the tolerance of a face of the other over a stretch without lying on it | InvalidInput | `touches` |
| the bodies touch at an angle too shallow for the tolerance | InvalidInput | a holder's test point |
| a contact of the bodies lies within the tolerance of a narrow piece of a face | InvalidInput | a holder's test point |
| booleans where an edge of one body touches a curve where the faces meet (a crossing tangent to the edge; design note 0012 §18) | Unsupported | settling |
| boolean: a contact edge carries no edge of the other body; a contact with no piece of the other body at it; a contact's pieces disagree on their side; one point touched from both sides; an edge left with no origin | InvalidResult | the stage that finds it |

The non-manifold texts and "a point where the bodies meet lies within
the tolerance of an edge it is not on" are the existing ones.

**Cost.** `touches_in` runs only on pairs with a zero flag, over the
coedges the face's index already finds near the surface; flat faces skip
locating. The second walk runs only when a bridge exists, and isolated
nodes join the refinement's point tree. `check_contacts` is one pass
over the loops; obstacle-aware test points are computed for holders
only. The number of places is capped like the other work limits and
refused as `TooLarge` beyond it.

### 18.3 Slice 2: surfaces that touch (design)

Faces leave no zero marks where their interiors touch, so detection
comes from `meet`; places, the rule, obstacles and the origin filter are
slice 1's, with a neighbourhood view for contacts on face boundaries.

**Contact geometry.** `Meet::Touch` gains the contact: a point, a line, a
ray from an apex, or a circle, with the two sides (the case of §18.1)
and the relative normal curvature Δκ across it, over a band of
half-width h = √(2ε/Δκ) about which the surfaces stay within ε of each
other (about 32 µm for a rod of 5 mm on a plane at the default
tolerance, 45 µm for a ball of 10 mm, 14 µm for one of 1 mm; along a ray
h grows from zero at the apex). Each site of `meet` that returns a touch
fills in what it already computes, by the closed forms of Miller (1987)
and Miller & Goldman (1992), deciding on lengths and angles, never on
discriminants: a plane, cone or cylinder about a sphere's axis (a point
at the pole or a circle), two spheres (a point on the line of centres),
meridians meeting on the axis (the apex or pole), a plane along a
cylinder (a line), a plane through a cone's apex (a ray along a ruling,
or the apex alone), and parallel cylinders (a line). Two corrections
come with it. *Nearly coincident* surfaces, where the band is wider than
the faces (Δκ ≤ 2ε/L², L the faces' extent across the contact), are
`Same` when the largest gap over the faces is at most ε and otherwise
`InvalidInput`, "two faces of the bodies nearly coincide …", as for
planes: spheres of nearly equal radii, and parallel cylinders of equal
radius whose axes are between ε/16 and ε apart, return a touch today. A
*varying gap* that leaves the band inside the faces is `InvalidInput`,
"the faces are nearly tangent over a stretch".

**Localised to the faces.** After the hull and wrap exemptions, the
refusal of touching surfaces in `relation` becomes a localisation,
padded by ε + h so that a real contact is never dropped: a point
farther than that outside either face, or a line, ray or circle missing
either face's box padded the same, is dropped, and a pair with nothing
left is apart. This removes the refusals made on surfaces alone, where
the surfaces touch outside the faces: by the census, 168 of 740 surface
touches in the random curved pairs, 149 of 1,096 in the random
revolves, and one first cut of a cone.

**Stretches of contact.** `spans_on` cannot find stretches on a tangent
track: its direction rule needs the surfaces to cross, touch cuts net to
zero, and an open contact line gets no spans at all (found by review).
So the boundary stations on the contact become clusters (`add_crossings`
both ways); per face, a line, ray or circle's ends are the face's
boundary clusters on it, each interval between ends is decided by
locating its midpoint, a contact with no ends by the probes of
`spans_on`'s branch for that case, and the two faces' intervals are
intersected (`both`); a point item is located in both faces. Each
stretch is added with a new touch source, ordered after the meetings,
merging with an operand edge or a crossing track along it (the
operand's edge stays primary); a point item becomes a raw touch point
before clustering, and its cluster joins both faces' isolated points.

**Inside both faces** (a rod on a plate, a ball on a plate, two balls, a
ball in a bore) the stretch is a bridge, or the point an isolated node,
in both faces; both pieces are holders and the slice-1 rule applies
unchanged, a bridge's source being an edge of the other operand or the
touch. A closed contact (a ball in a bore along a circle) separates the
faces, so `check_closed` decides it and `merge_faces` rejoins the
pieces where there is no trace. A sphere's pole and a cone's apex are
vertices of their faces: there the point is a vertex contact, and
`smooth_at` no longer exempts it, since the point now comes from `meet`.
*Inside one face, on the other's boundary*: a rod lying on a plate along
its top edge is tangent to the top face along the edge while the
plate's side face crosses the rod there. In the rod's face the contact
is a bridge with sources the plate's edge, the side's meeting and the
top's touch, and the rule refuses the union although `check_closed`
sees only two uses (found by review).

**On both boundaries.** Where the contact is a segment of both faces'
graphs, a view of the rays in the plane normal to the contact decides
it at its middle and a quarter from each end: the radial ordering of
faces about an edge (Hoffmann, Hopcroft & Karasick 1989), with
second-order terms breaking the ties first order leaves. A face with
the point inside gives two opposite rays, a face with it on its edge
one ray into the face; each ray carries its operand, its material side
and its normal curvature in the section (from a new
`Quadric::normal_curvature`: 0 for a plane, 1/R for a sphere, (w·c)²/r
for a cylinder, (w·c)² cos α / ρ for a cone, c the circumferential
direction). Rays sorted by angle tie by construction for the two faces
of a touch, and curvature orders them; other rays within an angle φ with
φρ ≤ ε (ρ the reach) are refused as meeting too shallowly, and
curvatures equal within the band as nearly coincident. Each operand's
rays alone alternate in and out of it (else `InvalidResult`), each
sector of the merged order takes the operation's membership, and k, the
number of rays with different membership on their two sides, gives no
surface (0), a manifold (2) or a pinch (4 or more). For k = 2, rays at
angle π are a G1 edge, accepted; at angle 0 a *cusp*, refused with
`InvalidInput` ("the result would have an edge where its faces meet at
zero angle (a cusp)"; a cusp keeps a knife thinner than the tolerance
over about √(2rε), which §6 refuses); otherwise a crease. If both of an
operand's faces at an edge are tangent to the other surface the view
does not apply and the contact is `Unsupported`. For the obround, box ∪
cylinder of equal height with the box's sides tangent: at the box's
edge the rays are the box's side S along +x, its end E along −y, and the
cylinder's C₊ along +x curving down and C₋ along −x, sorted C₊, S, C₋, E.
The union has S and C₋ at π, a G1 edge (the end face crosses the
cylinder along the same ruling, so the edge exists already); box −
cylinder has S and C₊ at 0, a cusp, refused; cylinder − box and the
intersection have creases. The G1 edge's origin is the box's edge, with
the touch of the side merged; the filter drops the end face's meeting,
which bounds nothing there.

**Near a boundary.** A stretch inside a face within ε + h of the face's
boundary is allowed only along an operand edge or ending across a
boundary edge transversally (a node `check_vertices` decides); otherwise
the band reaches an edge the contact does not run along, refused with
"a point where the bodies meet lies within the tolerance of an edge it
is not on". Every touch segment, separating ones included, is an
obstacle with band 2h, so the gap at a test point is at least 4ε (for a
ray, h at the stretch's far end).

**Still refused after slice 2:** crossing tangencies (a figure-eight, an
apex on a surface it crosses, a touch point on a crossing curve); nearly
coincident surfaces and gaps that vary across the band; cusps; the
boundary-band rule; a face whose loop would touch itself at a tangency
(a pocket's rim tangent to a face's edge from inside), `Unsupported`
until the walk, tessellation and the Euler check are tested on it;
curved edges tangent to curved edges in coplanar faces; and quartic
tangencies. The width of the band in which a crossing turns into a
refusal and then a touch is measured on curved goldens moved in steps
of ε/4 and recorded as built; a reference across that flip may fail, as
design note 0013 §7 accepts.

### 18.4 Slice 3 (outline)

1. **Tangencies refused today as quartics.** A separating plane or the
   gap is tested first, then the closed forms of Miller (1987) feed
   slice 2's touch items: skew rods with axes r₁ + r₂ apart touch at a
   point, with a band of √(2ε(r₁ + r₂)) / |sin θ| along it; a ball
   against a rod's side at d = R + r or d = r − R touches at a point (d
   = R − r is a figure-eight, refused); a ball against a cone off its
   axis; a cylinder along a cone's ruling outside it gives a ray (inside
   it a residual conic crosses the ruling, refused); cones sharing an
   apex give a ray, with the nearly coincident test where their angles
   are close.
2. **Crossing tangencies**, traced through the singular point with its
   branch directions (Bajaj, Hoffmann, Lynch & Hopcroft 1988) as a node
   of four branches in each face's arrangement, ordered by curvature:
   a station on a crossing curve (settling, step 1), a vertex where a
   curve leaves tangent to one of its edges (left today to the split's
   ordering and refusals), figure-eights.
3. **Nearly coincident internal contacts:** the threshold between `Same`
   and a refusal, by its own ADR.
4. **Curved edges tangent to curved edges in coplanar faces:**
   `add_edge_contacts` handles straight edges only; closed forms for arc
   and arc, and arc and line, would give the arrangement the node.
5. **A self-contact check in validation:** faces within ε that share no
   entity, found through a tree over faces and computed exactly for
   analytic faces, with its own time budget.
6. **Cusps**, if they are ever accepted: zero-angle corners in
   tessellation and in the split's ordering of first chords.

### 18.5 Tests

Goldens go in a new `tests/touching.rs` (`tests/contacts.rs` holds the
sketch curves' contacts). Each asserts the outcome of all four
operations in both orders of the operands, the exact text and location
of each refusal, and the counts, volume and area of each result, with
results that leave no trace compared with the apart case (§18.1).

**Flat operands.** P = `cuboid((−2, −2, −1), (2, 2, 0))`, (6, 12, 8),
volume 16. W, the wedge: the triangle (0.3, 0), (0.8, 0.5), (−0.2, 0.5)
in the frame {origin (0, 0.5, 0), x = X, y = Z, z = −Y}, extruded 1:
(5, 9, 6), volume 0.25, its edge x = 0.3, z = 0, |y| ≤ 0.5. C, the cube
on its corner: side 0.5, corner at (0.3, 0.2, 0), axes (√(2/3), 0, 1/√3),
(−1/√6, 1/√2, 1/√3), (−1/√6, −1/√2, 1/√3), all three edges rising:
(6, 12, 8), volume 0.125. The goldens:

1. W on P (case i): P ∪ W refused, "touch along an edge near
   (300.000000, 0.000000, 0.000000) mm"; P − W = P and W − P = W; P ∩ W
   empty.
2. C on P (case i): P ∪ C refused, "touch at a point near (300.000000,
   200.000000, 0.000000) mm"; P − C = P, C − P = C; P ∩ C empty.
3. Goldens 1 and 2 on P with a round hole of radius 0.2 through
   (−1.5, −1.5), whose top face is not flat: the answer no longer
   depends on an unrelated hole.
4. W with its edge through the centre of P's top face: P − W = P, which
   needs the test point that avoids the contact.
5. W inverted inside P (case ii), the triangle (0.3, 0), (0.8, −0.5),
   (−0.2, −0.5): P ∪ W = P, P ∩ W = W, P − W refused at (300, 0, 0) mm,
   W − P empty.
6. A square pyramid nested at its apex (case ii): base 1 × 1 at z = −1,
   apex (0.3, 0.2, 0), in `cuboid((−2, −2, −2), (2, 2, 0))`: union the
   box, intersection the pyramid, box − pyramid refused at the apex,
   pyramid − box empty.
7. A groove's root (case iii): A = `cuboid((−0.5, −0.5, −1), (0.5, 0.5,
   0))`, B = `cuboid((−2, −2, −2), (2, 2, 1))` minus the prism of the
   triangle (−0.5, 1), (0.5, 1), (0, 0) in xz along y: A ∪ B = B,
   A ∩ B = A, A − B empty, B − A refused along the root.
8. A notch onto the face (case iv): N = `cuboid((−1, −1, −0.5), (1, 1,
   0.5))` minus the prism of (−0.25, −0.5), (0.25, −0.5), (0, 0) in xz
   along y: P ∩ N refused along the root; P ∪ N (11, 24, 16), volume 18,
   with no edge on the root's line; P − N volume 14.25 and N − P (6, 12,
   8), volume 2, as main gives.
9. The slit: H = P minus `cuboid((−1.2, −1.2, −2), (−0.8, −0.8, 1))`,
   (10, 24, 16), volume 15.84; K the knife, the triangle (−1, 0),
   (−0.8, 0.3), (−1.2, 0.3) in {origin (−2.5, 0, 0), x = Y, y = Z,
   z = X} extruded 1.5, lying from beyond P's edge to over the hole:
   H − K = H with (10, 24, 16) (main gives (10, 27, 18)), H ∪ K refused
   at (−1600, −1000, 0) mm, K − H = K, H ∩ K empty.
10. Contacts on P's boundary: W's edge running out past P's edge (a
    bridge attached at one end); C's corner at the middle of P's top
    edge and on P's corner; W's edge on P's top edge, the difference
    keeping that edge's one origin.
11. `bodies_touching_along_an_edge_or_at_a_point` gains counts and
    origins for both differences.
12. Nested, edge on edge: a wedge inside the box of golden 6, its edge
    along the box's top front edge, A ∪ B = B with A first, B's edge
    keeping only B's edge as origin.
13. Review's counterexamples: CE1, a block with a through V-groove whose
    root (x = 0.7) lies inside a plate's top face, A − B refused at the
    root as on main; CE2, a wedge hanging from the top of P's box, its
    end flush in the side, B − A refused; CE3, a knife edge from a
    post's re-entrant vertical edge lying on the plate, the union
    refused at the knife; CE3′, the same knife between two posts, P − A
    two clean holes, equal to P minus the posts with the knife raised.
14. The reflex corner: an L-prism's inner corner at the centre of a slab
    face with normal (1, 1, 1)/√3, all four operations valid with main's
    counts and volumes.
15. Edges crossing: W's edge across an inverted wedge's edge at one
    point; the union refused at the point, the differences the operands.
16. A crossing and a touch in one face: a post through P carrying a beam
    and, under the beam, W resting on P; each operation compared with
    the beam and wedge raised.
17. Many contacts in one face: five wedges and a 4 × 4 grid of corners.
18. Goldens 1, 2 and 9 at millimetre and kilometre scale.

**Curved goldens** (slice 1's curved pairs; rod R of radius 0.5 along x
over [−1, 1], axis at y = 0.3, z = 0.5): a box turned 45° about y
resting by its lowest edge across R's top ruling (a station touch);
C's corner on a ball off its poles and seams (the cut-pair filter must
keep `assign` from refusing); a box's edge tangent to the ball; W with
its edge along R's top ruling; a coin tilted 30° resting by its rim on
P; a cylinder cut at an angle (its rim a conic NURBS, ADR 0023) resting
on P; a box's edge in a rod's cap plane tangent to its rim from outside;
a box's corner on a turned cap's centre (the `merge_edges` question); a
crossing tangency (a plane tangent to a rod's rim while crossing its
side in an ellipse tangent to the rim), `Unsupported` in all four
operations; the knife ring (a plate's hole rim resting on a ball); the
turned knife rim on P; and a coin whose rim seam lies 5 µm from its
lowest point, a stretch within the tolerance, refused, with the seam
1 mm away accepted. Tests that state today's refusals are reviewed and
given their new verdicts.

**Unit tests.** The split: an antenna, an isolated segment, a tree, a
connector between the outer loop and a hole, points inside, in a hole
and outside, a point within reach of a chord (`Close`), a bridge with a
side (`Inconsistent`), and bit-identity with no bridges. The cut-pair
filter: a cancelling pair, a crossing, a run over two clusters. The
table of §18.1: all sixteen cases and operations, from the holder's and
incident pieces' classes to the keep flags and the verdict.

**Fuzz.** A contact family in the flat fuzz: wedges, corners, pyramids,
V-notches, grooves, reflex corners, knives and posts on boxes' faces,
edges and vertices from either side, the generator knowing each
contact's σ, τ and convexity and moving it by g from −2ε to 2ε in steps
of ε/4: apart for g > ε, the table for |g| ≤ ε, and for g < −ε the
volume identities or a named sliver or closeness refusal, outcome
classes changing monotonically in g. The existing fuzzes are
re-baselined (silent pinches now refuse as `InvalidInput`; any
`Unsupported` in the flat fuzz stays a bug), and at every commit the
outcomes are diffed by class, operation by operation, against the
parent commit, with every change explained in that commit's message.

**Oracles** (through the public interface only). *Silent pinch*: for
every accepted result in the fuzz and goldens, every vertex against
every face that does not hold it, and every edge against every face
that does not hold it beyond a quarter of the edge's length from the
vertices they share, at most ε/2 apart fails; exact for planar faces
from their meshes (their polygons) and for curved faces from the
surface the test built, through the result's face origins. *Spurious
refusal*: every refused pinch, moved 4ε apart along the contact's
normal toward the side that separates, must succeed with the shells the
table predicts.

**Mutations**, each failing a named test, recorded as built:
`check_contacts` returning `Ok`; bridges pruned again; the flat guard
restored; `touches_in` finding nothing; obstacles off; the origin
filter off; the side flipped; the settling rule's station step off; the
cut-pair filter off; the stretch refusal off; the hole test in placing
points off; whole lying loops taking the seam path again.

**Time** (`timing.rs`, release): a comb of 512 teeth lying on P, once as
bridges ending inside the top face and once as chords across it, and
256 corners on a 16 × 16 grid, each timed for P − comb and the refused
P ∪ comb, within twice the time of the same comb sunk 0.1 mm into P (the
crossing case with as many features), with a scaling check from 128
teeth to 512; the existing budgets must not move. Numbers are taken
with callgrind first and recorded in docs/PERFORMANCE.md.

### 18.6 Order of work

Each step is a commit that passes the four checks on its own:

1. This section and ADR 0024, with pointers from §3, design note 0011
   §4.4 and ADR 0020.
2. The origin filter in `build`, with goldens 10 to 12 and every change
   of origin in the existing tests reviewed.
3. Bridges: the split's loose items and second walk, places,
   `check_contacts`, the cross-checks, obstacles and the avoiding test
   point, the triage of bridges carrying only meetings; goldens 1, 4, 5,
   7 to 9 and 13 to 17. This fixes the edge half of the flat fault and
   the slit.
4. Points on planes: the zero flags, runs' sides and clusters, the
   stretch refusal, `touches_in` without the flat guard, settling and
   isolated points; goldens 2, 3 and 6 and the flat contact family; the
   re-baselined rates.
5. The cut-pair filter, with its debug assertion and the near-apex
   goldens reviewed.
6. Curved pairs: `refuse_touches` deleted, `stretches` calling
   `touches_in`, whole lying loops along their own segments, the curved
   goldens, the `Kernel::boolean` documentation (an edge or a vertex
   touching a face is decided, with no trace or refused as non-manifold;
   surfaces touching and crossing tangencies are still `Unsupported`),
   and the curved fuzz counts.
7. The Apart band, measured first.
8. The sweeps, time guards and mutations; the as-built section; STATUS,
   ROADMAP (slice 1 ticked only when built and tested) and
   docs/PERFORMANCE.md.

Slice 2 follows the same way: the touch geometry at every site with the
two corrections and localisation (a localised touch still refused, but
at the contact); the touch relation, contact stretches, touch sources
and points, touch places and obstacles (goldens: rod on plate, ball on
plate, two balls, ball in a bore, rod along the plate's edge, a plane
along a cone's ruling, a cone's apex on a plane); the ray view, the
boundary-band rule, G1 edges accepted and cusps and self-touching loops
refused (the obround with its origins in all four operations, a slot
tangent to a wall); and the curved sweeps, the measured band, the oracle
on touch results and the time guards.

**Literature cited in this section.**

* A. A. G. Requicha, "Representations for rigid solids: theory, methods,
  and systems", *ACM Computing Surveys* 12(4):437–464, 1980.
* R. B. Tilove & A. A. G. Requicha, "Closure of Boolean operations on
  geometric entities", *Computer-Aided Design* 12(5):219–220, 1980.
* M. Mäntylä, "Boolean operations of 2-manifolds through vertex
  neighborhood classification", *ACM TOG* 5(1):1–29, 1986.
* D. H. Laidlaw, W. B. Trumbore & J. F. Hughes, "Constructive solid
  geometry for polyhedral objects", *SIGGRAPH* 1986, 20(4):161–170.
* C. M. Hoffmann, J. E. Hopcroft & M. S. Karasick, "Robust set
  operations on polyhedral solids", *IEEE CG&A* 9(6):50–59, 1989.
* J. R. Miller, "Geometric approaches to nonplanar quadric surface
  intersection curves", *ACM TOG* 6(4):274–307, 1987.
* J. R. Miller & R. N. Goldman, "Using tangent balls to find plane
  sections of natural quadrics", *IEEE CG&A* 12(2):68–82, 1992.
* N. M. Patrikalakis & T. Maekawa, *Shape Interrogation for Computer
  Aided Design and Manufacturing*, Springer, 2002, §6.4.1.
* C. L. Bajaj, C. M. Hoffmann, R. E. Lynch & J. E. Hopcroft, "Tracing
  surface intersections", *CAGD* 5(4):285–307, 1988.

## 19. As built: slice 1 on flat faces (2026-10-16)

Steps 2 to 4 of §18.6 are built: the origin filter, bridges and places,
and points on planes. Steps 5 to 7 are not: the cut-pair filter, curved
pairs (a contact where a face is curved, an edge or a vertex touching a
curved face or a curved edge resting on a plane, is still refused as
`Unsupported`, "a tangent contact", as on main) and the Apart band. Of
step 8, the contact fuzz and its sweeps are built and this section
records them; the contact time guards (ADR 0024: a comb of 512 teeth, a
grid of 256 corners) are not, and the mutation checks were run on the
backstop's rules only (below). There are 31 numbered goldens, golden 11
in `tests/booleans.rs` and the rest in `tests/touching.rs`: §18.5's
eighteen flat ones and thirteen found while building (19 to 31), by
review and by the contact fuzz. Every fix below has a test that fails
without it, but for one backstop that no test reaches (*Bridges with no
edge*, (a)).

**Where the build differs from §18.2.**

* *Bridges off the face are dropped*, not refused as "an imprint edge
  lies outside the face": the imprint of a coincident face takes the
  pieces of its edges within the tolerance of this face, and one can
  reach just outside it (a copy of a body moved by just over the
  tolerance, found by the fuzz of motions). Points off the face were
  already dropped.
* *"One point touched from both sides"* is checked only for the isolated
  points that settling keeps. A point that ends a segment of the touched
  face is dropped first (rules 2 and 3), and at a vertex where a face
  crosses through a reflex corner the runs on either side of the
  crossing report both sides, correctly.
* *Holders* take up to four test points (`HOLDER_TRIES`) from
  `classify::interior_points_avoiding`, large pieces as well as small:
  the triangulation of the region with the obstacles' points as free
  vertices, each inside triangle's incircle scored against the
  obstacles. `probe_point` with each obstacle as two extra sides is not
  built.
* *No separate limit on places.* A place is a segment or a cluster, each
  once, so the existing limits on segments and points bound them; the
  cap refused as `TooLarge` that §18.2's "Cost" names is not built.
* *Operand faces whose loop passes a vertex twice* (a face of an earlier
  union whose boundary touches itself at a vertex, a valid face) are
  classified: their triangulation takes the vertex's two visits as one
  point, where main refused "two vertices at one point of the face's
  plane". Two vertices at one point are still refused.
* *Bridges that carry no edge* of the other body, which §18.2 makes a
  kernel error, are decided where the fuzz found them (below), and a
  bridge whose faces leave the face to both sides, a kernel error in
  `wedge`, likewise. Contacts whose pieces disagree on their side are
  refused as input where the tolerance puts the place on two faces of
  the holder's body away from their edge.
* *A backstop on every result* (`check_apart`, below), which §18.2 does
  not have.
* *Golden 6* nests a cube by its corner (`corner_cube`) in place of
  §18.5's square pyramid by its apex; the outcomes asked for are the
  same.

Refusals added beyond §18.2's table, each followed by "near {p}":

| text | kind | where |
|---|---|---|
| an edge of one body lies on a face of the other and leaves it at a point the tolerance cannot place | InvalidInput | `imprint::plane_meeting`, `edgeless_contact`, `assemble::check_contacts` |
| a corner of one body lies within the tolerance of two faces of the other away from their edge | InvalidInput | `check_contacts` |
| the bodies leave a sliver thinner than the tolerance (the existing text) | InvalidInput | `edgeless_contact`, `merged_chord` |
| an edge of one body runs within the tolerance of a face of the other over a stretch without lying on it (§18.2's text) | InvalidInput | `astride` |
| the result would be non-manifold: the bodies touch along an edge / at a point (the existing texts) | InvalidInput | `check_apart` |

**Faults on main found by the new tests.** The contact fuzz at a
hundred times its default run (6,000 cases, 264,000 operations), with
the final oracle, on main (b032fc6): 24,571 operations fail. On this
slice none does.

1. *Silent pinches* (the fault §18 was written for): 20,457 results, in
   3,693 of the 6,000 cases, that validate and touch themselves (the
   oracle), 15,813 of them with B within half a tolerance of the
   contact and 3,549 with it exact. The oracle reports a point within
   half the tolerance of a face; nearer the edge of the band main
   returns more that it does not see: on golden 22's blade nested 0.75
   tolerances under P's top face, P − blade validates with the blade's
   edge 0.75 tolerances under that face, and goldens 20 and 21 come
   back as one solid whose two parts are joined by a sliver 0.75
   tolerances thick.
2. *Origins that are not along their entity*: 2,148 results naming an
   operand entity that does not bound the result entity (2,014 "edge …
   is not along … edge …", 134 "edge … leaves the plane of face …").
   Two causes: the sources a cluster on a lying edge adds (the origin
   filter of §18.2 fixes those), and joined edges. `merge_edges` joins
   two segments at a point where only they meet, and the joined segment
   kept the sources of both pieces, so an edge four tolerances long that
   ended at a corner, joined to the long meeting beyond it, named the
   meeting as a piece of that edge. The joined segment now keeps only
   the sources that run along all of it (`join_sources`): those of both
   pieces, and those of one piece that the other piece's far end and
   middle are within the tolerance of (a straight edge's line, or the
   surfaces of a curved edge's faces or of both faces where the bodies
   meet). Where no kept source bounds the faces along the segment, the
   point stays a vertex. Golden 19, `a_wedge_leaving_an_edge_near_a_corner`.
3. *Shell volumes far from the origin*: 1,391 kernel errors "a void …
   is inside no solid" and 535 "outer shell is inside-out", 1,644 of
   them with B pushed four tolerances into A or out of it.
   `assemble::shell_volume` summed the signed volumes of the faces'
   cones from the world origin. Its faces, projected onto their planes,
   fail to close by up to the tolerance times their perimeter, and the
   sum was off by that times the distance from the origin, so for a
   shell a few tolerances across, a metre or less from the origin, the
   sign came out wrong: a sliver was taken for a void, or a shell for
   one turned inside out. A fan over the places, tried first, gave the
   pillow slivers between a face and its moved copy (two faces, four
   edges) no volume at all. It now takes each face's plane relative to a
   corner of the shell: v = Σ (o_f − c) · n_f A_f / 3, exact for a
   closed flat shell and free of the origin's distance. The validator's
   `shell_volumes`, from a mesh, lost the sign to rounding there (a
   triple product of points a metre away) and takes a local origin
   likewise.
   Tests: `corners_dipping_four_tolerances_through_a_face` (contacts 1,
   2, 15 and 42), which fails on main, and
   `slivers_between_a_face_and_its_moved_copy` (motions 836, 4539, 5883
   and 6751 of the fuzz of motions run a hundred times longer), which
   main passes and the fan failed.
4. 40 kernel errors "the kept pieces do not close along an edge", in 9
   cases, at steps −4, −1, 0.75 and 1. None remains on this slice; they
   were not traced one by one.

The fuzz of planes and cylinders at twenty times its run found one more
(`corners_of_boxes_resting_on_flat_faces_of_prisms`, in `tests/curved.rs`):
a box's corner resting on the flat side of a D-shaped prism from inside
it (pair 3493), where main returned the prism minus the box touching
itself at the corner. The curved fuzzes do not run the oracle (below),
so main passed it; the difference is now refused as touching at a point.

**Loud failures the contact fuzz found in the slice**, at a hundred
times its run, each class fixed and kept as tests that fail without the
fix:

* *Bridges with no edge* (185 operations in 33 cases, at 50 of their
  steps, "a contact edge carries no edge of the other body"). Three kinds. (a) An edge E of
  the other body lying on the face F within the tolerance and leaving F
  across an edge K of F at a slant, its faces at a shallow angle to F:
  the curves where F meets E's two faces merge into one segment, which
  ends where K crosses their planes, a few tolerances along K from where
  E's own pieces end, so the bridge is E's contact moved by the
  tolerance. `lying_pieces` recognises it (the bridge carries the curves
  of two faces of the other body that share a straight edge lying on F,
  with one end of the bridge on that edge) and places the contact on
  E's own pieces from that end; where an end of E lies on the side of F
  opposite to the side E's faces leave it to, the faces cross F a few
  tolerances from E, and the boolean is refused as an edge leaving a
  face where the tolerance cannot place the point (a backstop: the check
  of the meeting's ends below refuses the cases found so far before a
  bridge forms, and no test reaches this one). On main the nested
  and crossing forms come back touching themselves (golden 22). (b) The
  curves where F meets two faces of the other body with no common edge
  on F, merged over a stretch: a sliver of F, refused. (c) A corner of
  the other body within the tolerance of F next to F's boundary, where
  one of the curves at it merged into the corner and the bridge dangles
  from it: the same sliver refusal (`lone_corner`). All of a face's
  edgeless bridges are judged before any of its contacts is placed
  (`edgeless_bridges`), so the bridge that shows the cause decides (a
  tree of bridges at a corner, contact 161); a bridge none of these
  explains stays a kernel error, raised after the face's other contacts.
  `contacts_whose_bridges_carry_no_edge`, golden 22.
* *Meetings along a lying run that end off it* (57 operations in 24
  cases, at 31 of their steps, "an edge with no origin"). Where one face's interval of the
  planes' common line is a run of its edges lying on the other's plane,
  the stretch is that run (§18.2), so an end of it that the other face's
  boundary gives must be a point of the run. A run within the tolerance
  of a plane can be a few tolerances from the planes' exact common line
  where they meet at a slant, and the other face's boundary crosses that
  line, not the run: the meeting then ran beside the run, a copy of the
  edge on the other body's faces that no face of the first bounds.
  `plane_meeting` now refuses such an end on the run's face's side as an
  edge leaving a face where the tolerance cannot place the point; off
  that side the copy lies off the face and is dropped. On main (where
  the copy was kept) goldens 20 and 21 come back as one solid joined by
  a sliver 0.75 tolerances thick, which validates. Review found three
  more forms of the same check: the side measured against the line from
  a run's first vertex to its last, where the run bends within the
  tolerance (golden 28, now against the run's edge nearest the end); an
  end of a meeting that goes on past the run's reflex corner into the
  face (golden 29); and both faces' intervals being runs, the faces then
  within the tolerance of each other over a band, which classifying
  their pieces refuses (golden 30).
* *Lying edges whose faces leave the face to both sides* (24 operations
  in 7 cases): a box's corner within the tolerance of an edge of the
  other body that crosses one of the box's faces at a slant a few
  tolerances from the corner, the piece of that edge between them
  imprinted as lying on the face. `astride` refuses it as an edge running
  within the tolerance of a face over a stretch (golden 23); where one
  order of the operands reached such an edge before a face found within
  the tolerance of another over a band, the error now waits until every
  face is split, and the band's refusal comes first in all four
  operations (contact 2080). A third form the triage found, a prism's
  edge parallel to an edge of the face and about the tolerance from it,
  is refused by the run check above (golden 24).
  `contacts_whose_lying_edges_leave_the_face_to_both_sides`.
* *A contact's pieces disagreeing on their side* (4 operations, contact
  3003): an edge lying on a face within the tolerance that leaves the
  body through a neighbour meeting the face at an acute angle before it
  reaches their edge, where it pierces the neighbour within the tolerance
  of the face but more than the tolerance from the edge. The tolerance
  cannot say which face the contact is on; refused as an edge leaving a
  face where the tolerance cannot place the point, and a corner touching
  two such faces as "a corner of one body lies within the tolerance of
  two faces of the other away from their edge" (golden 25). On main the
  same input's difference validates and touches itself along the edge.
  Review found the same disagreement where a blade's edge leaves a face
  three tolerances from its corner and the contact, placed to the point
  of the edge nearest the corner, took in a piece outside the face; the
  contact now stops where the edge leaves the face (golden 26).
* *Sections thinner than the tolerance*: a blade crossing a face near
  its corner whose two faces' curves on that face end within the
  tolerance of each other, so one segment carries both (found by
  review); refused as a sliver (`merged_chord`, golden 27), where main
  answered with the face's piece dropped and an edge named after a face
  that does not bound it.

**The backstop on every result** (`check_apart`). The rule after
selection decides the contacts the imprint records; this check decides,
from the result alone, those that decisions at the edge of the tolerance
let through. It tests the result's *exact* faces, those on planes whose
loops are all straight (faces with a curved side wait for the curved
slice): no end of a segment is on such a face without a corner there,
and no middle of a straight segment on such a face the segment does not
bound, else the result would touch itself and is refused as
non-manifold there. A point is on a face when it is within the tolerance
of it: within the tolerance of the face's plane and, in the plane,
inside the face's loops or beside them with the distances off the plane
and in it together within the tolerance.

Faces with a corner at an end of a segment are let off for its middle
where they are a corner the segment leaves (an edge leaving a corner of
a face at an angle θ stays within the tolerance of the face for ε / sin θ
from it):

* (a) a segment shorter than four tolerances with a face at both ends:
  only the face's inside counts, where the face holds the segment lying
  on it;
* (b) a segment shorter than four tolerances with a face at one end:
  let off;
* (c) a longer segment with a face at one end whose other end is
  farther than the tolerance from the face's plane: let off. A segment
  with an end on a face's corner and the other within the tolerance of
  its plane lies along that plane, and its middle counts.

Each was found by disabling the rules over longer runs of the contact
fuzz, and disabled one at a time, each refuses the cases that pin it
(`sharp_corners_the_check_of_places_lets_off`): (a) a piece of an edge
shorter than four tolerances between two corners of a face, its middle
outside the face but within the tolerance of its boundary (contacts 481
at 0.5 and 2896 at 0.75); (b) one shorter than four tolerances leaving a
face's corner (84 at 1, and golden 26); (c) a longer one leaving a
corner whose far end is more than the tolerance off the face's plane
(22 at −1, 281 at 1, 1335 and 5358 at −4, crossing the
face's plane at 14° and 3.6°, and 2714 at −1). The last is at the edge
of the band: an edge a tolerance under the face, whose far end the
imprint calls off it (1.0000000002 tolerances) and whose other end it
places on the face's corner; it may be answered or refused, and is
answered. A fourth rule, letting off a face at either end of a segment
no longer than two tolerances, was dropped: with it removed, no test
fails and the contact fuzz at a hundred times, the flat fuzz, the
motions and the curved fuzzes at twenty times refuse exactly the
operations they refused with it.

The check was first written letting off a face at either end of any
segment shorter than four tolerances; the fuzz run a hundred times
longer found three results where a piece of an edge 2.98 to 3.60
tolerances long, from a corner within the tolerance of the edge to
where the edge crosses a face at that corner, lies inside that face
(contacts 86, 2504 and 5245, `edges_lying_on_a_face_for_a_few_tolerances`;
golden 31 is contact 86 by hand). Two narrower rules were measured on
that run before (a): letting off only faces within the tolerance of the
middle refused 107 more operations, 66 of them a box's corner pushed
four tolerances through a face, whose triangle there has sides a few
tolerances long and sharp corners; testing a face at both ends for its
boundary as well as its inside refused three more, segments of 2.0 to
2.3 tolerances lying in a face's plane just outside it, 0.94 to 0.95
tolerances from its boundary. The test of a point's two distances
apart, first written, took points up to √2 tolerances from a face beside
its loops for on it, and refused 51 operations of the hundredfold run
that the test of the two together answers, each for a point between the
tolerance and √2 tolerances from a face
(`corners_just_beyond_the_tolerance_of_a_face`: contacts 3827 at −4 and
369 at −1, corners 1.08 and 1.31 tolerances from a face). Disabled, the
test of the two together fails that test. The check first ran only on
results made wholly of exact faces, so the same contact escaped where
any face of the result was curved (golden 31 with a round hole through
the cube's sides).

Cost: each place (a cluster, or a straight segment's middle) is tested
once. Each point is looked up in a tree of the exact faces' boxes, and
each face it may be on in that face's outline, built when first needed:
its sides in the plane, a tree of them, and a tree of its loops' boxes,
with a loop of more than 32 sides given a tree of its own, so a face
with many loops is not walked edge by edge for every point near it
(found by review: a plate with 1,600 flush pins took a quadratic time,
each point's parity test crossing every hole in its row). The faces
with a corner at each cluster and the faces along each segment are
grouped by a counting sort. A union of a box with a corner of a plate of 400 flush
pins is guarded (docs/PERFORMANCE.md).

**The oracle as built** (`tests/self_contact`). It finds a result's
edge within half the tolerance of a face that does not bound it, other
than the faces at the edge's ends near those ends. Flat faces bounded
by straight edges (exact faces) are measured exactly: the mesh of each
flat face is first moved onto its plane, since
its corners are the result's vertices, which may be up to the tolerance
off the plane, and a mesh leaning towards an edge a tolerance from the
face reported contacts that were not there (four reports of the contact
fuzz, the edge exactly a tolerance under the face). A straight edge is
measured whole against each exact face with no corner at its ends
(segment to triangle distances), so a crossing between sample points is
found (`the_oracle_measures_straight_edges_whole`, two diamond ridges
three tenths of a tolerance apart, made one body by a kernel a hundred
times finer); every other point of an edge (its mesh polyline's points
and the middles of their pieces) is tested against the faces near it,
other than a face at an end of the edge within the tolerance of the
point (every point that near a vertex is that near each face at it), and
other than a flat face that a straight edge leaves at a corner, its
other end more than the tolerance from the face's plane: as for
`check_apart`'s rule (c), such an edge stays within the tolerance of the
face for ε / sin θ from the corner, which is not a contact. Without the
vertex rule, fourteen of the twenty-one reports at the hundredfold run
were points within the tolerance of a vertex of the face; without the
corner rule, it reported contacts 1335, 2714 and 5358 of the pins in
`sharp_corners_the_check_of_places_lets_off`, and motion 3527 of the
fuzz of motions run a hundred times longer. Its own tests:
`the_oracle_finds_a_wall_thinner_than_half_the_tolerance` (a C whose
slot ends three tenths of the tolerance from its side, reported, and
one three fifths of the tolerance from it, passed) and the ridges
above. The half tolerance leaves contacts between half the tolerance
and the tolerance to `check_apart` on exact faces (main's golden 22
above). Curved faces, and flat faces with a curved side, are checked
through their mesh and confirmed on the surface; run on the curved fuzzes, the oracle reports 965 results,
the same on main, and its curved check counts a point of the surface
just outside a face, within the chord of its boundary, as on it, so
those reports are not yet sorted into contacts and false reports.
Sorting them, with an exact test of a point in a curved face, is part
of the curved slice; until then the curved fuzzes do not run the
oracle.

**The contact fuzz as built** (`contacts_swept_across_the_band`). A
random polyhedron A (a random prism of the flat fuzz) and a body B
resting on a random point of one of A's faces (inside it, on its edge or
at its corner): a box on a corner, its three edges rising more than
8.6° (the sine of each edge's rise above 0.15), or a triangular prism
on an edge, along an edge of the face or across
it, its faces rising at 9° or more, from outside A or inside it, moved
along the face's normal by −4, −1 to 1 by quarters, and 4 tolerances
(§18.5 planned −2 to 2). Every operation must pass the identities, the
origins check and the oracle, or be refused; where B rests from outside
inside A's face, A ∪ B must be refused as non-manifold at step 0, and
from inside, A − B; moved four tolerances clear, that operation must
not be refused as non-manifold. Sixty cases at the default run (2,640
operations), six thousand at the hundredfold run (264,000); CI runs it
at ten times, with the flat fuzz and the motions.

*What it does not assert.* §18.1's "no trace" (a result equal, entity
by entity, to the boolean with B moved clear) is asserted by the
goldens, not by the fuzz: B rests on a random point, often near an edge
of A, and leaves A's face across that edge, a second contact that moving
B clear also removes, so the counts with B moved clear differ for
reasons of their own. Within the band the counts do change: of the
83,065 operations answered both at step 0 and at a step within three
quarters of a tolerance of it, at the hundredfold run, 704 differ from
step 0's counts, 619 of them by one edge and one vertex (an operand's
edge left split near where B's edge meets it); 673 of the 704 differ
the same way on main, 29 are operations main failed, and 2 differ on
main by other counts. The results are valid and pass every check; they
are not yet canonical under motions within the tolerance.

**Known fault: a fin along a lying run between reflex corners.** A
prism's face F tilted 10° to a box's top face G, crossing it along a
line, with an edge of F lying on G's plane 5 tolerances beside that
line, between two reflex corners of F past which F's boundary turns
away and the line runs on inside F (the profile (−0.5, −0.5), (0, −5ε),
(1, −5ε), (1.5, −0.5), (2, −0.5), (2, 0.5), (−0.5, 0.5), of `tilted` in
`tests/touching.rs`, against `cuboid((−0.25, −1, −1), (1.25, 1, 0))`):
the strip of F between the line and the edge lies within the tolerance
of G, and A ∩ B is returned as the two wedges beyond the corners joined
by a fin 5 tolerances wide and under a tolerance thick along the edge.
It validates and passes the oracle (no edge of the fin comes within the
tolerance of a face it does not bound). Main returns the same
intersection, and in the union and the box minus the prism names
origins that do not bound their entities; this slice answers those two
correctly. Where the box's edge crosses the line beside the edge (the
box from x = 0.2 to 0.8, or from −0.25 to 0.5) every operation is
refused, as in golden 29. The end check of `plane_meeting` looks at the
meeting's ends beside a run; the meeting here passes the run with both
ends beyond it. Detecting it means deciding, for a meeting that passes
beside a run on its face's side, between this fin and golden 28's
reflex bend, whose triangle of F between the line and the run is
absorbed into the exact contact's outcomes; §20 does so and refuses
the fin.

**Rates.** Refusals, at a hundred times the flat fuzzes' default runs
and twenty times the curved ones', on main and on this slice:

| fuzz | operations | main: refused | main: failed | slice: refused | slice: failed |
|---|---|---|---|---|---|
| random polyhedra | 120,000 | 47 | 0 | 47 | 0 |
| motions by the tolerance | 28,000 | 10,004 | 0 | 10,004 | 0 |
| contacts swept across the band | 264,000 | 53,168 | 24,571 | 80,627 (305.4 per thousand) | 0 |
| random solids of planes and cylinders | 16,000 | 334 | 0 | 329 | 0 |

The contact fuzz refuses 5,573 of its 24,000 operations at step 0, the
exact contact (1,991 on main, which fails 3,582 more there), and 4,213
and 2,841 at four tolerances either way (−4 and 4). The other curved
families (cylinders at an angle, balls and cones about one axis, planes
across cones, second cuts, chains of cuts) refuse the same operations as
on main.

**Outcome changes in the other fuzzes**, against main: in the flat fuzz,
none; in the motions, the same 10,004 operations are refused, 70 of
them (in 19 motions) now as an edge leaving a face where the tolerance
cannot place the point, where main gave points closer than the
tolerance (34), a sliver (32) or a face within the tolerance of the
other's surface (4); with the final oracle, main's motions fail none.
In the planes and cylinders at twenty times:
pairs 1072 (a box's corner on the bottom of a rounded-square prism) and
1733 (a box's corner on the bottom of a disc with a square hole), from
inside, refused by main in all four operations as tangent contacts, now
give three results and refuse the prism minus the box as touching at a
point; pair 3493 is the fault on main above.
`corners_of_boxes_resting_on_flat_faces_of_prisms` asserts all three.

**Cost.** Instructions (callgrind, one run each) against main: a union
of two 256-sided star prisms, 295.5 million to 311.9 million (+5.6%:
`check_apart` 11.3 million, `build` 5.3 million against 3.8 million);
64 pins cut one at a time, 2,168 million to 2,344 million (+8.1%:
`check_apart` 95.5 million, `build` 68.5 million against 49.6 million).
The guard and scaling check on the plate of flush pins are in
docs/PERFORMANCE.md. Found while measuring, on main too and not changed
here: `shape::shape` groups faces by surface with a linear search over
the groups found so far, so its time grows with the faces times the
distinct surfaces (14.5 times from 100 to 900 pins, against 9 times
the pins; 82.8 million of the union's 1,739 million instructions at
900).

**What remains of slice 1**: the cut-pair filter (step 5), curved pairs
(step 6), the Apart band (step 7), the contact time guards and the
mutation checks beyond `check_apart` (step 8), and the oracle's exact
test of a point in a curved face. The fin above is refused since §20.

## 20. A meeting beside a lying run (2026-10-16)

§19's known fault, taken up before the rest of slice 1.

**The fault.** Two faces on planes meet along their planes' common line
L at a slight angle θ. An edge E of face F lies on the other face's
plane P (its vertices within the tolerance ε of P) a distance d from L
in F's plane, between two reflex corners of F, past which F's boundary
turns away from P and L runs on inside F. F's interval on L is then
longer than E, so the meeting's stretch, the segment between the places
of its two end clusters, runs along L (§18.2 makes the stretch the run
only where the run is the whole interval) and passes beside E. The
strip of F between E and the stretch, d wide, lies within the tolerance
of P: its edges do (E within ε, the stretch on P) and F is flat. Where
P's face G holds the strip, F and G run within the tolerance of each
other over it, and the result holds a fin under a tolerance thick: A ∩
B with the box of §19 is one solid, two wedges beyond E's corners
joined along E. Probed with the profile of §19 (θ = 10°, the box from
x = −0.25 to 1.25), all six operations:

| E's offset d | outcome |
|---|---|
| 0.5 ε | answered as the contact: A ∪ B refused as touching along an edge, A ∩ B two solids |
| 2 ε, 5 ε | answered with the fin: A ∩ B one solid of 0.000918 m³ |
| 8 ε (E 1.39 ε off P: not lying) | refused as a sliver thinner than the tolerance |
| 5 ε, the box ending beside E or at a corner | refused where the stretch's end lies beside E (the end check, golden 29) |
| 5 ε, F's boundary starting at E's first vertex (one reflex corner) | answered with a fin from E's start, joined to the one wedge |

So the outcome jumps from the contact (d within ε, where the clusters
make E and L one place) to the fin and back to a refusal at 8 ε, and an
end of the box beside E already refuses.

**Why nothing sees it.** `plane_meeting`'s check (§19, *Meetings along
a lying run that end off it*) looks only at the stretch's ends, and here
both lie beyond E. Classification tests one point well inside each piece
of a face (§9), and the strip joins the two wedges into one piece of F
and one of G, whose test points fall in a wedge, far from the other
body. The result's faces meet only along their shared edges (E bounds
both the strip of F and that of G once the clusters merge them), so
`check_apart`, which tests places, and the oracle, which tests edges
against faces, find nothing: the fin is two faces within the tolerance
of each other away from any edge, the contact of slice 2 (§18.3).

**Options.**

1. *Imprint E as the meeting where it lies beside the stretch*, with
   short edges from L to E's corners, so that the outcomes are d =
   0.5 ε's (the contact). Rejected: those edges, d long, are on no curve
   where the surfaces meet; the strip is within the tolerance of P on
   both sides of them, so which side of the band the meeting runs on is
   what the tolerance cannot decide; and §18.2's precedent (the run as
   the stretch) holds only where the run is the whole interval, its
   ends being the interval's own.
2. *Test each piece of a face for a neck narrower than the tolerance.*
   The strip is d wide in F's plane; it is thin only in space, between
   F and G, so no test in a face's plane sees it.
3. *Test the result for faces within the tolerance of each other away
   from their shared edges.* The general test slice 2 needs (§18.3), on
   exact distances between faces; larger, and it would refuse only
   after the whole result is built.
4. *Test, in `plane_meeting`, the strip of F between E and L for the
   other face holding it.* Chosen: it decides the fault where it arises,
   from the faces' planes and boundaries, and refuses every operation
   before any face is split. (It began as a test of the stretch's ends
   and E's vertices beside the stretch; review showed that points
   cannot decide it, below.)

This keeps §6: a face thinner than the tolerance left by near-coincident
geometry is refused, never dropped and never returned (Segal 1990 makes
the same choice for polyhedra: what the tolerance merges is merged
consistently, and what it cannot is reported).

**Review of the first rule.** The rule as first written measured each
vertex of E inside the stretch's span and refused where the stretch
passed more than ε inside F beside E, except where the other face's
interval was itself a run. Three adversarial reviews (over-refusal,
under-detection, the plan) and probes of their cases found:

- *The exception let the fin through.* With the box's top edge on L
  (the box to y = 0) or 4 ε beside it (to y = 4 ε), that edge lies on
  F's plane and is G's whole interval, so the stretch is G's run, and
  G's material lies on E's side of it. Main and the first rule answer
  both with the fin (A ∩ B one solid, 0.000918 m³), and the one-corner
  profile against golden 29's box ending at y = 0 likewise. Narrowing
  the exception to both runs being whole, as the review suggested,
  would refuse the mirror image, the stretch passing 5 ε inside F beside
  E's vertices: the box from y = 0 (or 4 ε, or −0.5 ε) to 1, its edge on
  L with its material away from E, where F only rests on that edge, and
  main gives the contact's outcomes (unions refused as touching along an
  edge, intersections empty, differences the operands).
- *Nothing tested that G holds the strip.* Where E lies over a hole of
  G whose back edge is on L, or over a notch of G whose floor is on L
  or 0.5 ε beside it, or over the slot of a U, G's interval on L still
  reaches the stretch (the hole's edge lies on F's plane, so the face's
  line passes over the hole), the strip lies over no face of A, and
  main answers B as resting on that edge (the unions refused as
  touching along an edge; the hole: intersections empty, differences
  the operands; the notch and the U: intersections the two wedges in
  the box's prongs); the first rule would have refused them.
- *Stretch ends from F's own boundary.* The argument that vertices
  suffice assumed that the end check decides every end of the stretch
  beside E, but it checks only ends the other face gives. A profile
  with notches from F's top down to 2 ε beside L, 3 ε above E's
  corners, makes a stretch between F's own notch floors. Probed with
  six boxes, main refuses each (by classification, by the end check
  where the box ends beside E, or as a bridge leaving a lying edge
  where the box ends in a notch); the second rule measured both ends
  anyway.
- *The plan.* The new check must follow both faces' end checks, so that
  golden 29 keeps its refusal; the stretch is the segment between its
  clusters' places, which runs at a slant where an end is a vertex off
  L; the distance test is implied by the side test (the distance to an
  edge's line is at most that to the edge); and finding E's edge nearest
  each vertex by a scan is quadratic in the run's length.

Probing the reworked rule found more configurations that main answers
with B's material under G along the strip: a hole whose front edge is
on L, G's material on E's side (A ∩ B one solid, the wedges joined
along E); a hole from x = 0.2 to 0.8 on either side of L, G holding the
strip beside E's corners; and a notch of G whose floor lies 2 ε beside
L, between L and E (A ∩ B one solid: B's material between L and the
floor, 2 ε wide, joins the wedges, where with the floor on L it is two);
and the box's top ending on L, the box from x = −0.25 to 0.5, so that
the stretch ends beside E (A ∩ B one solid: a wedge and the strip to
x = 0.5).

**The second rule and its review.** The rule first built measured
points: each vertex of R within the stretch's span, at the stretch's
point s nearest it, and each end of the stretch within R's span; at
each it refused where s lay more than ε inside F across R's edge
nearest s and the point q of that edge nearest s lay in G or within ε
of G's boundary. It passed goldens 32 and 33 and the fuzzes, but a
second review (a code review, two fact-checks, and two searches for
fins let through and contacts refused, each finding re-checked by
another agent) showed that measured points cannot decide a region:

- *The side across the nearest edge.* How far the stretch passes inside
  F is linear along it across one edge's line, not across whichever
  edge of R is nearest. With short steps at E's corners nearly square
  to L (θ = 3°, E 18 ε off L, the steps from 1.5 ε to 18 ε), a step is
  the edge nearest each measured point and its line passes within ε of
  the point: A ∩ B is the fin, as on main.
- *The hold at one point.* Whether G holds the strip was tested only at
  q. Triangular nicks in G under E's corners (0.4 wide and 0.1 deep,
  their apexes 2 ε beside L toward E, and others), with G under the
  strip between them, put each q in a nick: every operation returns,
  and A ∩ B is one solid joined along E, as on main.
- *A hold within ε of G's boundary.* Counting q within ε of G's boundary
  as held refuses B resting over a notch or a hole of G whose walls are
  flush with E's corners and whose back edge is on L, or over a pocket
  whose far lip is within ε of E: the strip lies over the gap, and main
  answers these as the contact (unions refused as touching along an
  edge, A ∩ B the two wedges, the differences valid, the oracle clean).
- *G's boundary under the strip away from the stretch.* An edge of G
  crossing L lies within ε of F's plane wherever it is within ε / sin θ
  of L, so it can pass under the strip far from where it meets L. A
  single edge of G crossing L at a slant at x = 1.2, beyond E, with G
  on E's side of it, makes a stretch that does not reach E, and A ∩ B
  is the fin, as on main. The argument that G's boundary passes under
  the strip only at measured points was wrong.

Where the boolean is refused at more than one place, by this rule or
another check, the place reported can depend on the operands' order,
as it already did on main.

**The third rule and its review.** The third rule took each run R of F
lying on P that is not one of F's whole intervals on L, swept it to L
(each point to its foot there), left out the points within ε of L, cut
the rest to F, and refused where G held a point of that strip more than
ε inside it, found by the search below. Two reviewers checked the rule
as designed, one the search's completeness and one its outcome on every
configuration found so far, and the contact fuzz ran on the rule as
built. They found:

- *A meeting between the run and L.* Where an edge of G lies on F's
  plane between R and L and is G's whole interval, it is the meeting
  the imprint makes (§18.2): F is split along it, not along L, and the
  part of the sweep beyond it is a part of F's piece across the
  meeting, which classification judges with that piece. Swept on to L,
  the strip took in G's band between that edge and L, and the rule
  refused contacts that main answers: golden 32's box mirrored, from
  3 ε or 2 ε short of L (its edge 0.52 ε or 0.35 ε from F's plane) to
  the far side; a box from 20 ε short against the prism at 1°, E 40 ε
  off; and golden 33's notched block, its floor 5 ε off L, against
  golden 30's wedge, its edge 3 ε or 2 ε off (there the wedge's edge is
  the meeting and the floor is R). Main answers each as the contact:
  the unions refused as touching along the meeting, the intersections
  empty, the differences the operands. The contact fuzz, run at a
  hundred times its runs on the rule as first built, changed the
  outcome of five contacts from main's: contact 2100 at −0.75 ε is
  this case (refused in every operation where main answers it), and
  contacts 1456 at −0.75 ε and 4416 at −1 ε are answered as on main
  with either this fix or the restriction of each segment's sweep to the
  side where F lies toward L (found in building the fixes, below), and
  refused without both.
- *A hole lying whole.* A hole of F whose loop lies within ε of P
  throughout is not a run (a run ends where its loop leaves the plane),
  so the rule did not look at it. A plate at 1° with a hole from 50 ε
  to 20 ε off L (each vertex within 0.87 ε of P) gives with golden 32's
  box an A ∩ B of one solid, as on main, holding a flap under ε thick
  between the hole and L (a cut across the flap is refused as a sliver);
  with the hole reaching 200 ε off L, its near edge is a run and the
  rule refused. The outcome turned on a dimension that has nothing to
  do with the flap.
- *A lone reflex corner.* A vertex of F within ε of P whose two edges
  leave the tolerance on its side, at a reflex corner of F, is not a
  run: F's boundary touches P there without crossing it, and F lies on
  both sides of the corner's way to L. With the corner 5 ε off L at 10°
  (the apex of a notch in F) and the box's top under that way, A ∩ B is
  one solid on main, the wedges on either side joined by a neck under ε
  thick; cut flat by a lying edge 12 ε, 20 ε or 100 ε long, the same
  apex is a run and was refused. The rule had let the lone corner
  through because its results validate and the oracle finds nothing,
  which holds for golden 32's fin as well.
- *A curved edge at a run's end.* The filter for curved edges near the
  strip, as designed (boxes, and the distance from a circle's centre),
  took every arc leaving a run's end for one entering the strip. It was
  built on the convex hulls of the Bézier control points instead (below).
- The search held, with its text made exact: the curves tested are cut
  to the strip, not to the sweep; the edges of G near the strip are
  found by the strip's box grown by ρ; F's edges include those of a loop
  lying wholly inside the sweep, and a part is split wherever F's
  boundary meets it; L is placed from a point of P near the run, never
  from the planes' constants, which far from the origin would lose more
  than the rounding margin; and two circles at most 2ρ apart (with a
  margin) have their crossings tested themselves.

Building the fixes found two more. Golden 35's pocket was refused at its
lip: the pocket's side wall gives the box's top an edge square to L
lying on P, whose sweep is a sliver as wide as rounding, which the
search took for a region G holds. And the sweep of a segment on the side
away from F (F lying beyond it from L, as past the far edge of a hole of
F, or past the floor of a notch in F) lies outside F but for the segment
itself, which counts as within ε of F's boundary, and the rule refused
where G held the segment. An audit of every configuration probed so
far, run on the rule without the fix for this, found seven that main
answers and the rule refused in all six operations, each against a box
from x = 0.2 to 0.8: in golden 32's prism at 10°, a hole from 5 ε below
L to 0.3 above it and one from 0.3 below to 5 ε above, a U-notch whose
floor lies 5 ε below L, a tab whose top lies there, and the face beyond
an edge 5 ε above L, crossing L only in a tab past the box; and at 3°,
holes from 15 ε and from 5 ε below L to 0.3 above. The contact fuzz's
contact 3448 at 4 ε was refused in all four of its operations, where
main refuses only B − A, as non-manifold. A segment is now swept only
on the side where F lies toward L; the seven give main's outcomes in
every operation (golden 39 asserts the first hole, the U-notch and the
tab), and the contact main's refusal alone.

**The side.** Probing the holes and corners on the other side of L, and
with B turned over (its body below its face, so that F's outward normal
points away from the box), showed that the side of L decides what a
strip's layer is. A point of F's plane on one side of L lies off P on
one side of it, so the thin layer between F and P over a strip lies in
F's body on one side of L and outside it on the other; and it lies in
G's body on one side of L and outside it on the other. With golden 32's
prism above its face, the layer below L is in both bodies (A ∩ B holds
it, golden 37) and the layer above L in neither; with the prism below
its face, the layer above L is in B alone (B − A holds it) and the layer
below L in A alone. Main answers every case with the layer in F's body
with a fin, a flap or a neck in the result that holds the layer. Where
the layer lies outside F's body, main answers a hole or a corner
correctly: the result that would hold a crack (the unions, the layer in
neither) or a sliver (A − B, the layer in A alone) is refused as
touching along the strip's edge, and the others validate with the
oracle clean (golden 38). The rule as first extended tested every
chain on both sides and refused those. For a run on that side, main
fails every operation ("the faces at a lying edge leave the face it
lies on to both sides", an invalid result): there the refusal replaces
the failure, and runs are still tested on both sides.

**The rule (as built).** For each pair of faces on planes that cross,
after the end checks of every stretch of the pair, and for each face F
of the pair against the other, G, on the plane P:

- *The chains* are the parts of F's boundary lying on P that can bound a
  region of F between them and L: each run R of F's straight edges lying
  on P that is not one of F's whole intervals on L (§18.2 takes such a
  run itself as the meeting, so that no line is imprinted beside it, and
  classification judges what lies between: golden 28's reflex bend,
  golden 30's two runs); each loop of F lying on P whole, closed (a
  hole: an outer loop lies there whole only where F does, which the
  faces' relation takes for coplanar or refuses first); and each vertex
  V of F within ε of P whose two edges are straight and run away from P
  (the other end of each farther from P than V, on V's side), at a
  reflex corner of F: F's boundary touches P there without crossing it,
  and F lies on both sides of V's way to L.
- *The strip* of a chain is the region of F's plane swept by each of its
  segments to L, each point to its foot there, on the side of the
  segment where F lies toward L (on the other side the sweep leaves F at
  once; a segment whose ends are no farther apart along L than the
  rounding margin sweeps nothing), less its points within ε of L, and
  cut by the pair's stretches: where a stretch lies between the chain
  and L over a part of the chain's span along L, the strip there ends
  ε short of the stretch's line, the meeting the imprint makes, and
  where the chain lies on that line it has none. A vertex V sweeps the
  segment from V to its foot. Every point of a strip lies within ε of
  P: the chain does, L and each stretch are on P, and the distance from
  P is linear along each segment. Under the tolerance a strip is a
  region of F lying on P.
- *The sides tested.* A run's strip is tested on both sides of L; a
  hole's or a vertex's only on the side where the layer between F and P
  lies in F's body (above): with `across` = n_F × the direction of L in
  F's plane, a point t across L lies t (across · n_G) off P along n_G, so
  the layer runs from it along −sign(t (across · n_G)) n_G, into F's
  body where sign(t (across · n_G)) (n_F · n_G) is positive.
- The boolean is refused, with the text classification gives a piece
  within the tolerance of the other body ("a face of one body runs
  within the tolerance of the other's surface without lying on one of
  its faces"), near the point found, if G holds a point of a strip on a
  side tested more than ε inside it: the point's foot on P lies in G
  farther than ε from G's boundary. The faces would then lie flush over
  a region of the strip that neither the meeting's imprint nor
  classification separates: the fin of a run (golden 32), the flap
  between a hole and L, or the neck from a reflex corner to L (goldens
  34, 37 and 38).

Where G holds no point of the strip that deeply, the strip lies over a
gap of G or within ε of its boundary (the flush walls above), and the
contact is decided as before. R need not be beside a stretch: with
G's slanted edge above, the stretch is away from E and G still holds
the strip. Points within ε of L or of a stretch are left out because F
and G meet there, and points within ε of G's boundary because that
boundary is one place with the strip's edge or with R. Points near R
are not left out: a run 1.1 ε to 2 ε from L has a strip all within ε
of R, and main returns the fin there. The regions of F's plane within
ε of P beside L that are no strip are not refused by this rule. Each
kind below is answered as on main (goldens 39 and 40), and a notch tip
across L whose two edges cross L fails or is refused as on main, by the
operands' order (*A fault on main*, below): past a segment where F lies
away from
L, the region between the segment and L is outside F; at a convex
corner of F pointing toward L, its edges running away from P on its
side, F lies beyond the corner from L, and the wedge between F and P
ends at the corner in a tip under the tolerance thick, as every
crossing at a slight angle ends along L; at a convex corner across L
whose edges cross L, F's part across L tapers to the corner, a tip of
the same kind; and at a reflex corner one of whose edges crosses L, the
region between the corner and L is the end of F's part across L,
bounded by that edge, which joins nothing.

*Finding a held point.* The search runs in P's plane, each point of
F's plane taken to its foot on P: a linear map, so F's straight edges
stay straight, and G's boundary keeps its distances. Let D be the points
of G farther than ρ from its boundary, ρ = ε (1 + 2⁻¹⁰), just over ε so
that every point found is more than ε inside past rounding. If D meets
the strip's interior, take a component K of the common part. Its
boundary lies on the strip's boundary (the chain's segments, the band's
edge ε from L, the ends of the sweep, the lines ε from the meetings, F's
boundary inside the sweep) or on D's boundary, which is made of
segments at ρ from G's edges (their offsets) and arcs of radius ρ about
G's vertices. So K's boundary has a piece on one of those segments, or K
is bounded by arcs alone and has a corner where two of the circles cross
(its boundary turns the wrong way along each arc, so it must turn at
corners). The rule tests:

1. Each side of each convex piece of the strip, each edge of F that can
   enter the strip, and each offset, on both sides, of each edge of G
   near the strip: its parts in the strip that lie farther than ρ from
   every edge of G near the strip (each edge's capsule of radius ρ),
   split where an edge of F crosses or touches it, each at its middle.
   No edge of G comes within ρ of a part, so G holds all of it or none
   of it, and F likewise: if K's boundary has a piece on the segment,
   the part holding that piece is held throughout.
2. For each two vertices of G near the strip at most 2ρ apart, the two
   crossings of their circles of radius ρ: if K is bounded by arcs alone,
   one of its corners is such a crossing, which lies in K's closure at
   least ρ from every edge of G.

A tested point is a find if it lies in the strip (to a margin of 2⁻¹⁰ ρ
in F's plane, for rounding), in F or within ε of its boundary, and in G
farther than ε from G's boundary. So every strip G holds deeper than ρ
is found, and the boolean is refused only where G holds the strip
deeper than ε (with F's boundary taken to the tolerance): the difference
is rounding. Each piece of a strip is the part of one segment's sweep
beyond the band, a trapezoid or a triangle whose sides are the band's
edge, two lines square to L and the segment's own line, so no side
takes its direction from a short edge (a vertex's is a rectangle about
its way to L, twice the rounding margin wide), cut by the meetings over
its span: split square to L at each meeting's ends and where the chain
crosses the meeting's line, and each part cut by the line ε from the
meeting on the chain's side, so every piece stays convex. The edges of
G near the strip are those within the strip's greatest distance from P
plus 2ρ of F's plane whose boxes meet the strip's box grown by 2ρ (ρ
and a margin of ρ for rounding), found by one walk of G's tree over its
edges' boxes that keeps to both; the edges of F that can enter it,
those within that distance of P whose boxes meet the strip's box grown
by ε, by F's tree likewise (a face of few edges has no tree, and each
of its edges is tested by its box). Where a tested point lies against
F and against G is decided, for a face whose edges are all straight, by
its loops in its plane with trees over their sides and their loops'
boxes (the outline `check_apart` uses, §19), built when a strip first
needs that face and kept for the rest of the boolean: a point meets the
sides within ε of it, the loops whose boxes hold it, and of those loops
all the sides of one of at most 32 and, of a longer one, those its ray
to +x may cross, not every side of every loop. The
first version walked every edge of both faces for each point, and an
edge of 1,024 segments lying beside the line over golden 35's notch
took more steps than the limit below and was refused as too large. For
a face on a plane with arcs among its edges, it is decided by the face's
outline in its chart (§11), whose bands about the chords hold the arcs,
as classification decides it. The outline of `check_apart` takes an arc
by its chord, and the first build used it for every face: a round bar's
top held no point of the strip, so the bar's A ∩ B with golden 32's
prism was the two wedges joined by the fin, and a round hole under the
strip held every point the slab's top held, so that contact was refused
(found by review; golden 41).

*Not decided.*

- A curved edge of G that may come within ρ of the strip, or a curved
  edge of F that may enter it deeper than the rounding margin, judged by
  the convex hulls of the control points of its rational Bézier pieces,
  which hold it (each against each piece of the strip, by separating
  axes): the offsets of curves are not segments, so the search above
  does not cover them. Such a boolean is refused as not supported
  (golden 41).
- A boundary of F within ε of P that is no chain: a curved edge lying
  within ε of P, or a corner one of whose edges is curved. This rule
  does not see it.
- A run on the side of L where the layer lies outside F's body is
  refused where main failed every operation; answering it, as a hole or
  a corner there is answered, would need §19's check of a lying edge to
  accept the faces leaving it to both sides there.
- The refusal does not depend on the operation, as classification's
  refusal of a piece within the tolerance does not: golden 32's A ∪ B,
  which holds the layer inside both bodies, is refused with the rest.
- Faces within 10⁻⁶ of square to each other (the cosine of the angle
  between them): a run within ε of P lies at most ε / sin θ from L, so
  their strips are under 10⁻¹² ε wide, and they are not tested.

*A fault on main found here, not fixed.* A notch of F from its far side
whose tip lies 3 ε or 4 ε across L from F's part beyond L, within ε of
P, its two edges crossing L (golden 32's prism at 10°, and at 1° with
the tip ten times farther off, against golden 32's box): with the box
first, every operation fails ("a contact's pieces disagree on their
side", an invalid result); with the prism first, every one is refused
as edges passing within the tolerance of each other. Main does the
same. The tip is no chain (its edges cross L), so the strips are not
involved; the failure is in deciding the contact's side (§19), and is
left for the next step.

*Cost and limit.* A face whose boundary has no mark within ε of the
other's plane (the faces' relation knows it) has no chain and is not
searched. A chain within ε of L in F's plane at every vertex, as every
run of an exact contact is, has no strip and costs a pass over its
vertices. Holes and corners are found among F's coedges within ε of P
(from the face's tree), each loop tested whole at most once. Otherwise
a chain's strip is at most one convex piece per segment and side tested
(a vertex's one), cut by the meetings over its span: the pair's
stretches are sorted once along L, with the greatest end so far, when a
chain first needs them, and each chain takes those over its span by two
searches, so a meeting beside no chain costs nothing beyond the sort.
Each meeting splits a part at most three times (at its ends and where
the chain crosses its line), so k meetings over a piece leave at most
3k + 1 parts. Each line is clipped against the pieces whose boxes it
meets, cut by the capsules and split by the edges of F whose boxes it
meets, each found by a tree over their boxes; each part tested costs a
lookup in each face's outline (a step for each side and loop it meets,
or, in a chart outline, one and one for each chord it looks at). The
searches of one boolean count their steps together (a stretch, an edge
or a cut taken, a box tested, a clip, or a side, a loop or a chord met
is one; building a face's outline is one for each side) and refuse the
boolean as too large past 2²² of them (`Limits::strip_steps`).

Profiled with callgrind on the zigzag of the time guards (an edge of n
segments lying beside L over golden 35's notch, zigzagging between 2 ε
and 4 ε off L: a run, and a reflex corner at each of its n / 2 vertices
2 ε off, n / 2 + 1 chains), the strip checks took 25.6 million
instructions at n = 256 and 224 million at 1,024, 8.8 times for four
times the segments: each chain's search walked every edge of both faces
near P, a walk that grew 15 times, with the square of n. With each walk
kept to the strip's box (by the face's tree), they take 16.5 and 67.4
million, 4.1 times, and the intersection (`Kernel::boolean`) 92 and 399
million. The time guard (`tests/timing.rs`) holds the intersection at
256 segments to 120 ms (it takes 11.6 to 13.1 ms), and its growth from
128 to 1,024 segments to 24 times (it grows 7.5 to 8.1 times); with the
walk undone it grew 9.7 times (measured in review) and would pass, so a
unit test of the boolean holds its steps to 300 a segment at 256 and
1,024 segments, where the walk undone takes about n² / 2 steps and is
refused at 1,024. The end checks of §19 found the faces' lying runs again for
each stretch of the pair and tested every pair of the faces' intervals
on L, as on main; they now find the runs once per pair, sorted along L
with their ends, and meet each interval of one face only with those of
the other within ε of it. On a face crossing the other's top 1,024
times (the comb of the time guards), the meeting of each pair, the
faces' intervals on L included, took 216 million instructions, 206
million of them in finding the runs, and takes 12.8 million, 2.5
million of them in finding the runs. A unit test of the boolean counts
the coedges looked at in finding them, and the sides looked at in
placing points in regions' cycles (`planar::regions`, which on the comb
walked whole cycles running along all the crossings; docs/PERFORMANCE.md,
eighth pass): each grows at most six times from 256 to 1,024 crossings
(with the runs found for each stretch, 15 times; with long cycles walked
whole, 24 times).

*Tests.* Goldens 32 to 41 (`tests/touching.rs`) assert all six
operations on: the fin of a run and its cases from the reviews (32, 33,
34); the gaps beside E: a notch and a hole flush with E's corners, a
hole flush at one end, and a pocket whose far lip lies within ε of E
(35); a meeting between the run and L (36); a hole and a reflex corner
on the side where the layer lies in F's body (37); runs, holes and
corners on either side of L, with the prism above and below its face
(38); five configurations with F beyond its edge from L, three of them
among the seven the audit found refused (39); corners beside L that are
not apexes, convex across L or pointing at it and reflex with an edge
crossing L (40); and faces with arcs (41): a round bar's top holding
the strip (refused) and a round hole under it (answered as the
contact), F turning away from P past the run's end along a quarter
circle that leaves the strip (refused as golden 32), a round hole of
G whose rim passes 0.9 ε beyond the strip (refused as not supported,
or with the prism first as a tangent contact, found first) and 1.5 ε
beyond it (refused as golden 32), and a round hole of F whose
rim dips 2 ε across L into the strip (refused as not supported) or
half an ε across, staying in the band (answered). Unit tests check the
search against a grid (above), and the strip itself: its band, its
sides, its cut by a meeting (ε short of the meeting's line, split where
the chain crosses that line, nothing where the chain lies on it), and
that cutting each piece by only the meetings over its span gives the
same pieces as cutting it by all of them, on 400 random chains and
meetings, at the cost the sort promises; a unit test checks the walk of
a face's tree kept to a box against filtering the slab's edges by their
boxes, and one the edges of F the check takes as entering the strip:
with a hole of F in the strip over G's region, every line tested but
those edges has its middle outside G or in the hole, and the check finds
the region's corners only through the hole's edges. Tests of the
boolean check that its strip checks draw on one budget (the fewest
steps within which k strips are answered are at least k times those of
one) and count the work of the three costs profiling found (above).
Of 27 changes to the rule and its searches (each test of a chain, a
side, the band, a cut, a window or a search left out or turned round,
alone or two together), the first run's tests caught all but three:
the skip of a face with no mark within ε, which only skips work; a
chain on a meeting's line on the other side of L, now caught by the
unit test of the cut; and the walk kept to a box, now caught by its
unit test. The placing of points in long cycles was caught then only by
the time guards' comb, whose boolean then fails, and a face crossing a
block's top 64 times along one line (`tests/booleans.rs`) now catches
it too. Run again on the final code with nine more changes (the search
the stretches over a chain and the runs beside an end now share, left
empty; and one for each fix of the last review and of the profiling:
arcs placed by their chords, the test of a curved edge near a strip
forced either way, F's edges entering the strip left out, a budget for
each strip, the runs found again for each stretch, long cycles walked
whole, and the walk's box taking every edge), the tests catch 35 of the
36, all but the skip of a face with no mark.
The flat fuzz, the motions and the contact fuzz at a hundred times their
runs, and the curved fuzzes at twenty and at fifty, refuse exactly what
main refuses (90,686, 8,290 and 20,522 refusals, the same operations
with the same texts), with no failures. Of the five contacts the third rule had changed in
the contact fuzz, four are refused only where main refuses them (their
other results pass the fuzz's checks), and contact 2080 at
4 ε, which main refuses as a face running within the tolerance (found
by classification), is refused with main's text at main's place (the
third rule found it at another).

References: M. Segal, "Using tolerances to guarantee valid polyhedral
modeling results", *Computer Graphics* (Proc. SIGGRAPH) 24(4), 1990;
Jackson 1995 (§2).
