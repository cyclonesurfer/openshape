# ADR 0003 — Persistent topological naming (direction; full design before Phase 2/3)

Status: accepted as direction; detailed design note due before Phase 3 work starts

## Context

Features reference faces, edges and vertices of earlier results (fillet these
edges, sketch on this face). Those references must survive upstream edits:
dimension changes, reorder, suppress, pattern-count changes. This is the
classic persistent naming problem (Kripac 1997; Capoyleas, Chen & Hoffmann
1996; Marcheix & Pierra 2002; OpenCascade's TNaming; FreeCAD's toponaming
work by Realthunder).

## Options

1. **Index-based** references (face #7): trivially broken by any change. Rejected.
2. **Geometric** references (closest face to a point/normal): survive some
   edits, ambiguous and fragile under large changes.
3. **History-based generative names**: each kernel operation tags every
   output entity with how it was generated (feature id + generator role +
   source entity names, e.g. `sweep(f42).side(profile-edge e3, piece 0)`,
   `sweep(f42).cap(end)`); splits/merges add disambiguating qualifiers;
   references are stored as queries over those names with a geometric
   tiebreak. Used (in different forms) by most commercial systems.

## Decision

Option 3 with geometric disambiguation as a fallback, implemented in a
`naming` crate. The kernel's sweep operations already build topology in a
fixed, documented order from profile edges (one side face per profile edge
per ≤π piece, caps last), which is the hook for generator tags. Feature IDs
derive from op IDs, so names are globally unique and stable across replicas.

## Consequences

* The kernel must return per-entity generation records from every operation
  (Phase 3, together with booleans and the `regen` crate).
* Tests: reorder, suppress, pattern-count change and dimension sweeps on a
  bracket model, asserting no broken references (Phase 3 milestone).
