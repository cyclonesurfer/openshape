# Design note 0013 — Persistent naming of faces, edges and vertices (Phase 3)

Status: design, before implementation. Details ADR 0003 (history-based
generative names with a geometric fallback). Booleans (design note 0012)
and every feature that refers to earlier geometry stand on this.

## 1. Problem

Features refer to faces, edges and vertices of the bodies earlier
features made: sketch on this face, extrude up to that face, fillet these
edges, mate on this face (Phase 5). A reference is made once, by a click,
and must keep pointing at "the same" entity when the model regenerates
after an upstream edit: a dimension change, a sketch edit that adds or
removes curves, a reordered, suppressed or deleted feature, a pattern
whose count changes. Bodies are rebuilt from scratch on every
regeneration, so entity indices mean nothing across regenerations.

"The same" is not always defined. A face can be split in two by a cut
that now crosses it, merged with a neighbour that became coplanar, or
vanish. The design's job is to follow references wherever the meaning is
clear, and to fail loudly, naming the reference and the feature, wherever
it is not: a reference that silently jumps to a different face is the
worst failure a parametric modeller can have, because the model still
regenerates and looks plausible.

## 2. Literature

* V. Capoyleas, X. Chen & C. M. Hoffmann, "Generic naming in generative,
  constraint-based design", *Computer-Aided Design* 28(1), 1996: names
  built from the generating operation and the names of the entities that
  produced an entity; disambiguation of split entities by the names of
  their neighbours.
* X. Chen & C. M. Hoffmann, "On editability of feature-based design",
  *Computer-Aided Design* 27(12), 1995: the editing operations under which
  references must survive, and why purely geometric references fail.
* J. Kripac, "A mechanism for persistently naming topological entities in
  history-based parametric solid models", *Computer-Aided Design* 29(2),
  1997: a topological ID system recording, per operation, which entities
  were created, split, merged and deleted, and matching after
  regeneration through that history.
* D. Marcheix & G. Pierra, "A survey of the persistent naming problem",
  *Proc. 7th ACM Symposium on Solid Modeling and Applications*, 2002: the
  taxonomy used below (naming, matching, and the split/merge/delete
  cases) and the observation that no scheme resolves every edit, so a
  scheme must say which it refuses.
* L. Agbodan, D. Marcheix, G. Pierra & S. Thomas, "A topological entity
  matching technique for geometric parametric models", *Proc. Shape
  Modeling International*, 2003: matching by structural similarity of
  names when exact matching fails.
* R. Bidarra, P. J. Nyirenda & W. F. Bronsvoort, "A feature-based
  solution to the persistent naming problem", *Computer-Aided Design and
  Applications* 2(1–4), 2005: naming through feature faces and their
  semantics rather than through raw topology.
* Our ADR 0003 (direction), ADR 0008 (periodic faces split at π), ADR
  0011 and design note 0006 §5 (region keys: the sketch-level instance of
  this scheme, already in use).

## 3. Principles

1. **Names are generative.** An entity's name says how it was made: the
   feature, the role in that feature's operation, and the names of the
   entities it was made from. Nothing in a name is a body-local index or
   a floating-point value, so names are identical on every replica and
   every regeneration (op-log replay must be deterministic, ADR 0004).
2. **The kernel reports origins, not names.** The kernel does not know
   features or sketch ids. Every operation returns, beside the body, an
   *origin table*: for every face, edge and vertex of the result, how it
   arose in terms of the operation's inputs (profile curve 3, cap at the
   end, face 12 of input A, the intersection of face 4 of A with face 9
   of B). The `naming` crate composes these with the feature's own maps
   (profile curve 3 is sketch edge key K; face 12 of A has name N) into
   names. Bodies stay opaque (CLAUDE.md), and the kernel is testable
   without the document model.
3. **Logical entities, not kernel pieces.** Faces split at π (ADR 0008)
   and edges split with them are a kernel representation choice, not
   something the user sees or selects. Names identify logical entities:
   the maximal set of pieces that share a name apart from the π-split
   qualifier. Selection, highlighting and references all work on logical
   entities, so changing a revolve from 170° to 190° (one piece to two)
   changes nothing a reference can see.
4. **References are queries with a hint.** A feature stores the name of
   the entity it refers to *as it was when selected*, plus a geometric
   hint (the picked point and the surface normal or curve tangent there,
   in model coordinates). Resolution matches names first; the hint only
   breaks ties among candidates the names cannot separate, and only when
   it lies on exactly one of them (§6).
5. **Refuse rather than guess.** Every case §6 cannot resolve fails the
   referring feature with a message that names the reference ("the face
   'Extrude 1, side of line 3' no longer exists"), and the features after
   it regenerate as far as they can, as for any failed feature.

## 4. Names

```text
Name        := Origin [Qualifier]
Origin      := Feature(fid) Role (Source*)
Role        := Box(±x|±y|±z)
             | Side | Cap(start|end) | Lateral        (extrude)
             | Side | Cap(start|end) | Axis           (revolve)
             | Kept(A|B)                               (boolean)
             | Instance(i)                             (patterns, mirror)
             | ...                                     (each later feature)
Source      := SketchEdge(EdgeKey) | SketchVertex(VertexKey)
             | Name                                    (an input entity)
Qualifier   := Split(Distinguisher) | Pi(k)            (§5)
```

Faces carry an origin from the kernel. **Edges and vertices** are named,
following Capoyleas et al. and Kripac, by the faces they bound: an edge
by the unordered pair of its adjacent faces' names, a vertex by the set
of names of the faces around it. Where this is ambiguous (two edges
between the same two faces, as on a full revolve), the kernel's edge
origin (piece boundary j, lateral edge at profile vertex V) is added; if
that is also ambiguous, an ordinal along a canonical direction is added
and marked *weak*, which §6 treats like a geometric tie-break.

Examples, for an extrude of sketch region R with distance d:

* the face swept by the sketch line with edge key K:
  `Feature(f) Side (SketchEdge K)`;
* the top cap: `Feature(f) Cap(end)`; caps of a multi-region extrude add
  the region key as source;
* the edge between that side face and the top cap:
  `{ Feature(f) Side (SketchEdge K), Feature(f) Cap(end) }`.

After a boolean, a kept face of A is `Feature(g) Kept(A) (N)` where N is
its name in A, so a reference to "the top cap of Extrude 1" still finds
it after a later cut, and keeps finding it when the cut is suppressed or
deleted (the name then is plain N again; §6 matches a stored name through
the `Kept` wrapper, which carries no information of its own). Intersection
edges are named by their two faces, as every edge is.

Sketch keys are reused unchanged: edge and vertex keys (design note 0006
§5) already name sketch geometry by sketch item ids that derive from op
ids, so they are stable across replicas.

## 5. Splits and merges

When an operation splits one input face into several result faces (a cut
crosses it), each piece's name gets a `Split` qualifier distinguishing it
from its siblings by the names of the faces that bound it and do not
bound the other pieces (Capoyleas et al.'s neighbourhood rule): the piece
of the top cap on the far side of a slot is the piece bounded by the
slot's far wall. If two pieces have the same distinguishing neighbours
(a face cut into strips by a pattern of identical slots, before pattern
instances are named), the qualifier is an ordinal and weak.

When an operation merges faces (a boolean's simplification, design note
0012 §4.5, or a fillet consuming a face), the result face carries the
names of all its constituents; a reference to any of them resolves to it.

π-splits (ADR 0008) get a `Pi(k)` qualifier that resolution ignores: all
pieces with the same name apart from `Pi` are one logical entity (§3).

## 6. Resolution

A reference is resolved against the bodies as they are just before the
referring feature, in this order:

1. **Exact.** The logical entities whose name equals the stored name. One
   match: done.
2. **Split.** No exact match, but entities whose names equal the stored
   name plus a `Split` qualifier exist: the referenced entity has been
   split since the reference was made. What happens depends on what the
   referring feature needs, which each feature declares:
   * a *carrier* reference (sketch plane on a face, extrude up to a face,
     a mate on a face, a direction along an edge) uses only the surface or
     curve the entity lies on, which all pieces share: resolved, and
     reported as nothing at all;
   * a *set* reference (edges to fillet, faces to delete or shell) takes
     all pieces, which is what the user selected: resolved;
   * a reference that must be one bounded entity takes the piece the hint
     lies on, if exactly one; otherwise it fails as ambiguous.
3. **Changed sources.** No exact or split match, but entities whose
   origin has the same feature and role and whose sources match the
   stored sources *by sketch item* (same entity id and direction, with
   different end-vertex keys because a sketch edit split or extended the
   curve): the sketch-level analogue of step 2, treated the same way
   (several candidates are pieces of one former entity).
4. **Merged.** An entity carrying the stored name among its merged
   names: resolved to it.
5. **Otherwise** the reference fails: "no longer exists". A candidate
   that merely shares the feature and role (another side face of the same
   extrude) is never substituted, however close the hint lies.

Steps 2 and 3 only ever choose among entities descended from the
referenced one, so the hint never moves a reference to an unrelated
entity. The hint is checked even when names resolve uniquely, only to
report: a resolved entity that the stored hint no longer lies on gets an
informational note in the feature's status ("Sketch 3's plane moved with
its face"), never a failure, since moving with the model is the point.

## 7. What changes are expected to keep references

| edit | effect on names | result |
|---|---|---|
| dimension change, same topology | none | exact |
| dimension change, a cut now crosses a face | `Split` added | carrier/set: resolved; single: hint or fail |
| sketch curve split or extended | edge keys change | step 3 |
| sketch curve deleted | its faces vanish | fail |
| later feature suppressed or deleted | `Kept` wrappers vanish | exact through §4's rule |
| referenced feature suppressed or deleted | names vanish | fail |
| reorder so the referenced entity is made after the referring feature | names not yet present | fail (the reorder itself is allowed; the failure names the cause) |
| pattern count increased | new `Instance(i)` names | old references unchanged |
| pattern count reduced below a referenced instance | name vanishes | fail |
| revolve angle crosses π | `Pi(k)` pieces change | exact (logical entities) |
| faces become coplanar and merge | merged names | step 4 |
| an edge vanishes (its faces became tangent or coplanar) | name vanishes | fail |

## 8. Data and APIs

* Kernel: every body-producing operation returns its origin table with
  the body; the table is indexed by the body's face, edge and vertex
  indices, which become part of the `Kernel` trait as opaque, per-body
  handles (`FaceIx`, `EdgeIx`, `VertexIx`) so the UI can pick and the
  naming layer can describe what was picked. Origins of the existing
  operations follow the construction order design note 0001 already
  fixes (one side face per profile edge per ≤π piece, caps last).
* `naming` crate (new): `Name`, origin composition, logical-entity
  grouping, `Reference { name, hint, use: Carrier | Set | Single }`,
  `resolve`, and the matching rules of §6. No kernel dependency beyond
  the trait's types.
* `features`: per feature, the maps from profile curves to sketch keys
  (already present as `SketchRegions::owners` and region keys) and the
  name table of every body it produced, kept in `StudioResult` so the UI
  can turn a pick into a stored name.
* Storage: references are part of feature definitions in the op log and
  the native format, versioned with it. Names are trees and grow with the
  number of features that touched an entity; they are serialised with
  structural sharing and their size is measured on the test models.

## 9. Tests

* The ADR 0003 bracket: a base extrude from an L-shaped sketch, a slot
  cut, a boss sketched on the top face, a hole through the boss, an
  extrude up to a face; references of each kind (carrier, set, single).
  For every edit in §7's table, and for random dimension sweeps over
  valid ranges, every reference must resolve to the entity an independent
  geometric oracle identifies (written per test: "the face with normal +z
  and the largest z"), or fail with the expected message. A reference
  resolved to the wrong entity is a test failure that blocks the commit.
* Property tests: random feature histories with random edits (dimension
  changes, suppress, reorder, delete); names are identical across two
  regenerations and across op-log replicas that merged the same ops in
  different orders; resolution never panics and never returns an entity
  that does not descend from the referenced one.
* Kernel: origin tables are complete (every entity has one), and
  consistent with topology (an edge's origin agrees with its two faces).

## 10. Order of work

1. Origin tables and entity handles for `make_box`, `extrude` and
   `revolve`, with completeness tests.
2. The `naming` crate: names, logical entities, resolution; unit tests on
   hand-built tables.
3. The first referencing features: sketch on a face and extrude up to a
   face (carrier references), with the picking UI.
4. Origins through booleans (design note 0012) and the split/merge
   qualifiers, then the bracket suite and the property tests.

## 11. As built (steps 1 to 3, 2026-10-12 and 13)

Steps 1 and 2 are done: ADR 0018 settles the origin tables and handles,
and the `naming` crate names box, extrude and revolve bodies and
resolves references against them. Where building it refined this note:

* **Where the types live.** `Name`, `Reference` and their parts are data
  stored in feature definitions, so they live in `core-model`
  (`opendatum_core_model::naming`) beside the sketch keys they are built
  from; the `naming` crate holds the algorithms and depends on the
  kernel only for the origin and handle types.
* **Caps** always carry the selected regions they cover, by the outer
  loop of the region key the feature stores (design note 0006 §5), not
  only in multi-region sweeps: a cap then keeps its name while the stored
  key keeps finding its region, whatever else the sketch edit changed,
  and adding a region to the feature leaves the other caps' names alone.
  Regions that share an edge become one cap carrying both; resolution
  treats caps (and side faces) that share a source as descendants, which
  is §5's merge rule and its inverse. A region selected anew whose
  boundary has other sketch curves has no recorded link to the old one,
  so references to its caps fail rather than guess.
* **Logical entities** are what the user sees: kernel faces joined across
  seams (the kernel's own judgement of where pieces of one surface meet
  with no crease, which the display already uses), so the `Pi(k)`
  qualifier is never stored and a sketch line split by a vertex where
  nothing else meets it still gives one face, named by both sketch edges
  (§5's merged name). Edges are the kernel edges between the same two
  logical faces joined end to end where nothing else meets them;
  vertices are where logical edges end. A body whose table would give
  one name to two entities is refused.
* **Tie-breaks** (§4) are the kernel origin in sketch terms: the sketch
  edge and end of a profile edge, and for a profile vertex the lesser of
  the two piece ends meeting there by sketch key, so they do not depend
  on the order of the profile's curves. Only where those also collide is
  the weak ordinal used.
* `Split` qualifiers and `Kept` wrappers come with booleans (step 4).

Step 3 (ADR 0019), model side, 2026-10-13: every part carries its names;
a sketch can sit on a planar face and an extrude can end on a planar face
parallel to its sketch, both carrier references resolved against the
parts made before the feature. The region profile reports the sketch
edge behind each curve, so side and cap names come from the same keys
the regions are named by; inline profiles get made-up keys from the
feature and the curve index. Hints are checked against a fine mesh
within 0.1 mm until the kernel can classify points against trimmed
faces (ADR 0019 §6).

Step 3, application side, 2026-10-13: the viewport picks the faces,
edges and vertices of the parts under the pointer (a vertex within 7
points first, then an edge within 5, then the face hit by the pointer's
ray; a candidate hidden behind a nearer face is skipped, as are seams
and anything the part's names do not cover), highlights the whole
logical entity while hovered and when selected, starts a sketch on a
selected planar face, and picks the face an extrude ends at. Step 3 is
complete.
