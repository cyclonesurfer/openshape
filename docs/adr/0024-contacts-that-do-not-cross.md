# ADR 0024 — Contacts that do not cross: decided in the face's split, pinches refused, origins named by the faces that bound them

Status: accepted, 2026-10-16 (design note 0012 §18; amends ADR 0020 §2,
and settles for touching bodies the "later decision" on non-manifold
results in design note 0012 §3).

## Context

Two bodies can touch without crossing: an edge or a vertex of one
resting on a face of the other (a cube on its corner, a wedge on its
edge, a rim on a plate, a groove's root on a face), or faces touching
faces (a rod on a plate, a ball in a bore, an obround's tangent sides).
The booleans had no rule for any of these. Curved faces refused such
contacts as unsupported "tangent contacts". Two flat faces skipped that
refusal and the planar split dropped an edge lying inside a face as an
antenna, so a plate's union with a cube standing on its corner, or with
a wedge on its edge, returned two solids that touch, and validated; a
round hole elsewhere in the plate made its face not flat and the same
input refused. A lying edge from a face's outer loop to a hole's loop
was kept as a slit: a holed plate minus a knife returned (10, 27, 18)
instead of (10, 24, 16), and validated. Both are faults on main, found
by review and confirmed on 2026-10-16, against `Kernel::boolean`'s
documentation, which says non-manifold results are refused.

Four questions had no answer in earlier records: where contacts are
detected; what a result that would be non-manifold becomes; which
origins a result's entities carry where a contact leaves no trace; and
how the smaller cases are classed (crossing tangencies, cusps, G1 edges,
dropped contacts). The exact result is not in question: booleans compute
regularised sets (Requicha 1980), closed under the regularised
operations (Tilove & Requicha 1980), and at a contact that does not
cross each operation gives no trace, nothing, or a pinch, each of the
four ways the bodies can meet pinching in exactly one operation (design
note 0012 §18.1).

## Options

### 1. Where contacts are detected

Three designs were drawn up and reviewed against the code.

* **A. At the imprint, from runs along loops**: read each loop's marks
  within the tolerance of the other surface, confirm a contact when
  every face at the touching entity shows a run, and filter a confirmed
  lying edge out of the imprint. It sees contacts where they are first
  known, but rests on the runs being found and confirmed, and review
  found three ways it fails silently: filtering the lying edge left a
  groove's root to topology, which returned two prisms touching along
  the root where main refuses; runs through a vertex where the body
  crosses the surface are skipped, so their lying edge reaches the split
  as an antenna and is dropped; and a knife between two such vertices
  leaves a slit.
* **B. In the face's planar split**: lying edges stay imprinted and
  touched points become nodes with no segments; the split reports
  bridges (chords with the same region on both sides, an exact fact of
  plane graphs) and isolated nodes, each held by the piece whose region
  contains it; one rule after selection decides them all, a pinch if
  and only if kept pieces of both operands are at the contact. Contacts
  on a face's boundary are shared nodes or segments, which the existing
  counts of uses and fans decide. Detection does not depend on which runs
  are confirmed, and with no contact the split is bit-identical to
  before; but completeness rests on point touches being recorded on
  every pair that can have them.
* **C. A local model of the neighbourhood, before the split**: cells
  about each contact from first-order rays, the outcome decided from
  them, a premise it cannot model refused. It gives two derivations of
  each class, and for surfaces touching on face boundaries its
  second-order view decides G1 edge, crease, cusp and pinch in one count,
  as a radial ordering of faces about an edge (Hoffmann, Hopcroft &
  Karasick 1989; Mäntylä 1986 for vertex neighbourhoods). But its
  premise refused a reflex corner crossing a plane, which main handles
  correctly; its absolute angle threshold added refusals; it decided a
  second time what the split already decides; and it kept A's run-based
  miss.

### 2. What a non-manifold result becomes

* **Refuse it** with `InvalidInput` at the contact, as ADR 0020 already
  says for bodies sharing only an edge or a vertex.
* **Represent it**: a radial-edge structure (Weiler 1986) or selective
  Nef complexes, which treat touching and coincident parts by one rule
  (Hachenberger, Kettner & Mehlhorn 2007). Every stage, the validator,
  tessellation and naming would need non-manifold entities.
* **Perturb the operands into general position** (expand for a union,
  shrink for an intersection; Ouchi & Keyser 2004). The result is
  manifold but has invented features: a bridge about 2√(2ε/Δκ) wide
  between tangent rods, or a slot through a wall, where the contact in a
  model is usually intended (a part resting on another).

### 3. Origins where a contact leaves no trace

* **A. Literal**: every source a segment collects is listed. A lying
  edge's segment carries the meetings of the touching faces with the
  touched face, and an edge on an edge of the other operand carries both
  operands' edges, so a result with no trace differs from the apart
  case; in a nested contact the other operand's edge can even become the
  primary origin, the first operand's sources sorting first.
* **B. Filter by the faces that bound the entity**: an origin names an
  operand entity only where it bounds a result face along the entity.
* **C. Filter only next to contacts**: needs each source to record
  whether a contact made it, which is fragile and adds state to every
  segment.

## Decision

**1. B, with grafts.** Contacts are detected in the face's split and
decided by one rule after selection. From A: gating by the zero flags
`relation` already computes, so bodies with no contact pay nothing;
cancelling, per cluster, a pair of opposite cuts, so touches on curved
pairs do not reach the track assignment's refusals; telling a station
from a vertex where a meeting of its own face reaches a touched point (a
station there is a crossing tangency, `Unsupported`; a vertex there is
an ordinary crossing through a reflex corner, left to the vertex fan
test, which keeps main's correct results); a guard for pairs called
apart while a mark lies within the tolerance, built only if measurement
finds it reached; and "no trace" defined as equality with the same
boolean with the touching body moved clear, since the stages normalise
some operands even when the bodies are apart. From C: independent
cross-checks (a probe across the touched face, and wedge convexity for a
lying edge), raised as `InvalidResult`; for slice 2, the ray view with
curvature breaking first-order ties; and the rule that a contact stretch
near a face's boundary runs along an operand edge or crosses the
boundary transversally. Neither A's lying filter nor C's premise is
taken. Added to B: an isolated point in a hole of its cell is off the
face; localisation of surface contacts is padded by the tolerance plus
the contact's band; and every contact segment keeps test points away,
including segments that separate.

**2. Refuse.** A result that would be non-manifold fails with
`InvalidInput`, "the result would be non-manifold: the bodies touch
along an edge near {p}" or "… at a point near {p}", at the contact. A
union of touching bodies is refused, not returned as a body of several
solids that touch: such a body touches itself, and later booleans assume
operands do not. Disjoint operands still unite into several solids.
Representing non-manifold bodies stays a later decision with its own
ADR; perturbation is not a default, and could at most be a named,
explicit operation ("fuse touching bodies").

**3. B, amending ADR 0020 §2:** an origin names only operand entities
that bound the result entity. For an edge, with U the union of the
sources of the result faces along it (merged faces and coplanar
partners included), `EdgeOrigin::Operand(op, e)` is kept if and only if
some (op, f) in U has e in its loops, and `Meet(fa, fb)` if and only if
(First, fa) and (Second, fb) are both in U. For a vertex, `Operand(op,
v)` is kept if and only if some result face at it comes from a face of
op holding v; otherwise the vertex is `VertexOrigin::Meet`. An edge left
with no source is `InvalidResult`, never silent.

**4. The smaller cases.** Crossing tangencies (a curve where faces meet
tangent to an edge at the contact, figure-eights, an apex on a surface
it crosses) stay `Unsupported` until slice 3 traces them through their
singular points (Bajaj, Hoffmann, Lynch & Hopcroft 1988); no flat input
reaches this refusal, so the flat fuzz keeps counting `Unsupported` as a
bug. G1 edges are accepted; cusps are refused with `InvalidInput`, since
a cusp keeps a knife thinner than the tolerance over about √(2rε).
Near-tangencies are `InvalidInput`: an edge within the tolerance of a
face over a stretch, gaps that vary across the band, nearly coincident
surfaces, contacts too shallow to classify. A face whose loop would
touch itself at a tangency is `Unsupported` until the walk, tessellation
and the Euler check are tested on it. Dropped contacts are not recorded;
a query for where two bodies touch is a separate feature. Design note
0011 §4.4's "an isolated tangent point is returned as a vertex" is
superseded: a contact inside a face never enters a loop, following the
precedent of not subdividing polygons that only touch (Laidlaw,
Trumbore & Hughes 1986).

## Consequences

* Unions that returned touching solids on flat faces are refused, the
  slit is gone, and an unrelated hole no longer changes the answer.
  Contacts of edges and vertices with curved faces, refused as tangent
  contacts, give results with no trace or a refusal as non-manifold. In
  slice 2, surfaces that touch only outside their faces are no longer
  refused (by the census, 168 of 740 surface touches in the random
  curved pairs and 149 of 1,096 in the random revolves). Bodies stored
  from before the fix that touch themselves fail to regenerate with the
  pinch message at the contact, which is correct.
* Origins change only in results where an operand keeps no face at a
  shared entity; every golden that changes is reviewed in the commit
  that changes it. The naming crate does not yet read boolean origins,
  so no stored reference changes.
* A curve where faces meet that dangles inside a face, dropped silently
  by the split before, is now an `InvalidResult`; the fuzz is triaged
  before that lands, and a legitimate class, if one appears, gets its own
  design.
* Records to change: the `Kernel::boolean` documentation; design note
  0012 §3 and §6; design note 0011 §4.4; a pointer from ADR 0020;
  STATUS.md and ROADMAP.md, with slice 1 ticked only when built and
  tested.
* Tests: goldens in `tests/touching.rs` with exact counts, origins and
  texts; an oracle that finds faces of a result within half the
  tolerance of each other away from shared entities, and one that moves
  refused pinches apart; contact fuzz families swept across the band in
  steps of ε/4; mutation checks; refusal rates re-baselined, with a
  diff of outcomes by class at each commit.
* Cost: contact work is gated by flags `relation` already computes, the
  split walks a second time only where a bridge exists, the rule is one
  pass over the kept loops, and release guards hold a comb of 512 teeth
  and a grid of 256 corners to twice the time of the same bodies
  crossing, with a scaling check and a limit on the number of contacts.

## As built: slice 1 on flat faces (2026-10-16)

Design note 0012 §19 records the build, with its departures from
§18.2. For edges and corners resting on flat faces, decisions 2 and 3
hold as written and decision 1 as built there; the cut-pair filter,
curved pairs and the Apart band (§18.6 steps 5 to 7) are not built, so
a contact where a face is curved is still refused as a tangent contact.
Three additions, each from a fault the tests found:

* **A backstop on every result** (`check_apart`): a place of the result
  within the tolerance of an exact face (planar, its loops straight) it
  is not on is refused as non-manifold there, so a contact that a
  decision at the edge of the tolerance let through is refused, never
  returned. Sharp corners a segment leaves are let off by three rules
  that §19 lists, each pinned by a test that fails without it.
* **Refusals for contacts the tolerance cannot place**: an edge lying on
  a face and leaving it where the tolerance cannot place the point, and
  a corner within the tolerance of two faces away from their edge, both
  `InvalidInput` with the place.
* **Origins of joined edges** keep only the sources that run along all
  of the joined edge (`join_sources`), which decision 3 implies but did
  not say.

Of the consequences, the time guards of contacts (the comb and the
grid) and the mutation checks beyond `check_apart` are not yet built;
the rates are re-baselined (§19: at a hundred times its run the contact
fuzz refuses 80,627 of 264,000 operations, where main refuses 53,168
and fails 24,571). One fault is known and not fixed: a fin under the
tolerance thick along an edge that lies beside a crossing between two
reflex corners, returned as on main (§19).

## References

* A. A. G. Requicha, "Representations for rigid solids: theory, methods,
  and systems", *ACM Computing Surveys* 12(4), 1980.
* R. B. Tilove & A. A. G. Requicha, "Closure of Boolean operations on
  geometric entities", *Computer-Aided Design* 12(5), 1980.
* M. Mäntylä, "Boolean operations of 2-manifolds through vertex
  neighborhood classification", *ACM TOG* 5(1), 1986.
* D. H. Laidlaw, W. B. Trumbore & J. F. Hughes, "Constructive solid
  geometry for polyhedral objects", *SIGGRAPH* 1986.
* K. Weiler, *Topological Structures for Geometric Modeling*, PhD
  thesis, Rensselaer Polytechnic Institute, 1986.
* C. M. Hoffmann, J. E. Hopcroft & M. S. Karasick, "Robust set
  operations on polyhedral solids", *IEEE CG&A* 9(6), 1989.
* C. L. Bajaj, C. M. Hoffmann, R. E. Lynch & J. E. Hopcroft, "Tracing
  surface intersections", *CAGD* 5(4), 1988.
* K. Ouchi & J. Keyser, "Handling degeneracies in exact boundary
  evaluation", *Proc. ACM Symposium on Solid Modeling and Applications*,
  2004.
* P. Hachenberger, L. Kettner & K. Mehlhorn, "Boolean operations on 3D
  selective Nef complexes: data structure, algorithms, optimized
  implementation and experiments", *Computational Geometry* 38(1–2),
  2007.
