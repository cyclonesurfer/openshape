# ADR 0018 — Origin tables and entity handles in the kernel interface

Status: accepted, 2026-10-12 (Phase 3; design note 0013 §8 and §10 step 1).

## Context

Persistent naming (ADR 0003, design note 0013) composes names from
*origins*: for every face, edge and vertex of a body, how the operation
that made it produced it, in terms of the operation's inputs. The kernel
reports origins and knows nothing of features or sketch keys; the
`naming` crate turns origins into names. Design note 0013 leaves three
questions to this ADR:

1. **Where origins travel**: beside the body or inside it.
2. **How entities are addressed** across the `Kernel` trait, which keeps
   bodies opaque.
3. **The origin vocabulary of the sweeps** (extrude, revolve, and the box,
   which is an extrude of a rectangle).

## Options

### 1. Where origins travel

* **A. Returned beside the body**: `extrude` returns `(Body, Origins)`.
  Explicit, but every body-producing operation and every caller changes
  signature, and a caller that drops the table cannot name the body's
  entities later.
* **B. Carried by the body**, read with `Kernel::origins(&body)`. The table
  is immutable data about the body's own entities, made with it and never
  out of step with it; operations keep their signatures. A body that
  enters a boolean brings its table, which the boolean's origins refer to.
* **C. Callbacks** from the kernel into a naming sink during the
  operation. Couples the kernel to the naming layer's types and order of
  work; rejected.

### 2. Addressing entities

* **A. Indices as handles**: `FaceIx`, `EdgeIx`, `VertexIx` (newtypes of
  `u32`), the positions of the body's faces, edges and vertices. They are
  already the order of `BodyMesh::faces` and `BodyMesh::edges` that picking
  uses, mean nothing outside their body, and reveal no kernel internals.
* **B. Opaque tokens** (hashes, generational keys). Nothing to gain while
  bodies are immutable values: an index into an immutable body cannot
  dangle.

### 3. Sweep vocabulary

A sweep moves a profile through *positions*: an extrude has two (0 at the
start, 1 at the end), a revolve one more than its pieces (ADR 0008 splits
it so that no piece turns more than π; a full revolve's last position is
its first). Profiles are curves in the caller's order; the kernel splits
some (arcs into pieces of at most π, closed splines in two). So:

* a **profile piece** is `(curve, piece)`: the curve's index in
  `Profile::curves` and the piece's index along the curve's own direction;
* a **profile vertex** is the pair of piece ends `(curve, piece, end)`
  that meet there, sorted (so it does not depend on loop direction or on
  which curve the walk started from; the naming layer orders the pair by
  sketch keys, so names do not depend on the order of the curves);
* faces: `Side { piece, span }` (swept by a piece between positions
  `span` and `span + 1`) and `Cap { at, region }` (the profile region
  `region` at the first or last position);
* edges: `Profile { piece, at }` (a copy of a piece at a position; for a
  partial revolve, a piece lying on the axis is one edge at position 0)
  and `Swept { vertex, span }` (the line or arc a vertex sweeps between
  two positions);
* vertices: `Profile { vertex, at }` (a vertex on the revolve axis is one
  vertex at position 0);
* regions are numbered in the kernel's order and described by the pieces
  of their loops, outer loop first, so the naming layer can relate a cap
  to the sketch regions it covers.

The alternative, feature-shaped roles (`Box(+x)`, `Cap(end)`, `Lateral`,
`Axis`), puts the naming layer's vocabulary into the kernel. The sweep
vocabulary is exact and small; design note 0013's roles are a direct
translation of it (a box side is the side of rectangle edge `i`; a
revolve's interior positions are its π-split `Pi(k)` qualifiers).

## Decision

1B, 2A and 3 as described. `Kernel::origins(&body)` returns the body's
`Origins` (faces, edges and vertices in handle order, the number of
positions, and the regions' pieces); `Kernel::adjacency(&body)` returns
each edge's faces, each vertex's edges and which edges are seams (where
pieces of one surface meet with no crease, as the display already
judges them), from which the naming layer joins kernel faces into the
faces the user sees and names edges by their faces and vertices by
theirs. `BodyMesh` gains each vertex's position, for picking vertices. Every operation fails
rather than return a body whose table is incomplete. Boolean, pattern and
later operations extend the origin enums with their own variants (kept
from input A or B, intersection of two input faces); the enums are
`#[non_exhaustive]` so that adding variants is not a breaking change for
the naming layer's catch-all handling.

## Consequences

* `Body` holds its origin table beside its B-rep; operations that copy or
  transform bodies carry or remap it.
* Tests (design note 0013 §9): the table is complete and unique (one
  origin per entity, no two entities alike), consistent with topology (an
  edge's origin agrees with its two faces, a vertex's with its edges),
  and right by an independent geometric check (a side face contains its
  piece swept to mid-span; a cap lies in the plane of its position).
