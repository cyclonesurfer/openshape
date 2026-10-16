# ADR 0020 — Booleans in the kernel interface, their origins and their first scope

Status: accepted, 2026-10-13 (Phase 3; design note 0012 §8 step 1,
design note 0013 §10 step 4).

## Context

Design note 0012 plans union, intersection and difference of two solids
in stages, the first being the whole pipeline on solids bounded by planes
(polyhedra: extrudes of polygons, boxes). Before code, three questions
need answers that the note leaves open:

1. **The interface**: the operation's signature, what an empty result
   is, and how the kernel says "not yet" for inputs a stage does not
   handle.
2. **Origins** (ADR 0018): every body carries a table saying how each of
   its faces, edges and vertices arose from the operation's inputs. A
   boolean's inputs are bodies, so its origins refer to their entities;
   the vocabulary must cover pieces (a face cut in two), merges (two
   coplanar faces joined, design note 0012 §4.5) and new entities (where
   the two bodies meet).
3. **What step 1 accepts**, and what it refuses.

## Options and decisions

### 1. Interface

* **A. One call per pair**, `Kernel::boolean(a, b, op)` with
  `op ∈ {Union, Intersect, Subtract}` (`a − b`), as design note 0012 §3
  writes it. Features that combine a new body with several parts call it
  once per part; the result's origins refer to exactly two operands.
* **B. N-ary** (one body against a list of tools). Saves repeated work
  on large patterns, but every origin then names an operand index and the
  first implementation gains nothing from it.

**Decision: A.** N-ary can come later as a wrapper or a separate call
when patterns need it, without changing this one.

An **empty result** (the intersection of disjoint bodies, or a difference
that removes everything) is not a body: the call fails with
`KernelError::InvalidInput` saying the result would be empty, which a
feature reports as its failure. Inputs a stage does not handle yet fail
with a new `KernelError::Unsupported` naming what is missing, so callers
and tests can tell "not yet" from "invalid" and from a kernel bug
(`InvalidResult`). Results that would be non-manifold (bodies touching
only along an edge or at a point) fail with `InvalidInput` saying so
(design note 0012 §3; for every contact that does not cross, design note
0012 §18 and ADR 0024).

### 2. Origins

The operands are called **first** and **second** (`Operand::First`,
`Operand::Second`): `a` and `b` of the call.

* **A. Per entity, one origin naming the operand entity it is (a piece
  of), plus a side table for entities made of several.** Faces:
  `FaceOrigin::Operand { operand, face }`. Edges:
  `EdgeOrigin::Operand { operand, edge }` for an edge that is (a piece
  of) an operand's edge, `EdgeOrigin::Meet { first, second }` for an
  edge where face `first` of the first operand meets face `second` of
  the second. Vertices: `VertexOrigin::Operand { operand, vertex }`, or
  `VertexOrigin::Meet` for a new vertex where the bodies meet (vertices
  are named by their faces, design note 0013 §4, so the kernel need not
  say more). `Origins::merged` lists the further origins of entities
  made of several: faces joined by the simplification of §4.5 (and a
  kept face that coincides with a face of the other operand), edges
  joined end to end, an edge lying on edges of both operands, a vertex
  of both operands.
* **B. Every entity carries a list of origins.** Uniform, but changes
  `Origins::faces`, `edges` and `vertices` from `Vec<Origin>` to
  `Vec<Vec<Origin>>` for every operation, sweeps included, to serve the
  minority of merged entities.
* **C. A history graph** (operation → entity → entities). What naming
  ultimately composes, but that composition is the naming layer's job
  (design note 0013 §4: `Kept(A) (N)`), and the kernel would duplicate
  it.

**Decision: A** (amended by ADR 0024: an origin names only operand
entities that bound the result entity). The existing tables and their
consumers are unchanged;
a boolean's origins say exactly which operand entities each result
entity comes from, and the side table carries the merges §5 of design
note 0013 needs. The naming layer reports boolean origins as unsupported
until naming through booleans is built (design note 0013 §10 step 4).
Origins are recomputed at every regeneration and are not stored, so the
vocabulary can still grow without a format change.

### 3. Scope of step 1

* Both operands must be bounded by planar faces (`Surface::plane`) with
  straight edges; otherwise `Unsupported`. That covers boxes and extrudes
  of polygons with holes; revolves and arcs come with step 2 (planes
  with cylinders, cones and spheres in closed form).
* Coincident faces (on one plane within the model tolerance, either
  orientation), shared edges and vertices, and edges lying in the other
  body's faces are part of step 1: they are the common cases in CAD
  (a boss sketched on a face, a cut flush with a side), not exotica.
* Refused with a message and a location, never guessed: geometry of
  the two bodies closer than the tolerance without coinciding (two
  faces within the tolerance over part of their extent but not on one
  plane; points that would be distinct vertices closer than the
  tolerance), and any internal inconsistency the stages detect.

## Consequences

* `Kernel` gains `boolean(&self, a, b, op: BooleanOp) -> Result<Body,
  KernelError>`; `KernelError` gains `Unsupported`; `Origins` gains the
  `Operand` and `Meet` variants and `merged`.
* The result passes the validator like every kernel result, and its
  shells and solids are rebuilt from scratch: a difference can leave
  several solids, or a solid with a void.
* Tests: the volume identities of design note 0012 §7 on random
  polyhedral pairs and rigid motions (including motions by fractions of
  the tolerance), golden cases with exact volume, area and topology
  counts, origin tables that cover every entity and agree with the
  operands' geometry, and a printed refusal rate.
* Features do not use booleans yet; that is design note 0012 §8 step 4,
  after the curved cases.
