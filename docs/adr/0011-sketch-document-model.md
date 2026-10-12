# ADR 0011 — Sketches in the document model

Status: accepted (Phase 2)

## Context

Sketches are the first feature whose content users edit interactively and
incrementally: draw a line, add a constraint, change a dimension, drag a
point. Each of these must be an undoable step, merge sensibly with a
collaborator's concurrent edits, and replay deterministically. Solving is
iterative, and its result depends on the starting configuration (a dimension
change has several solutions; users expect the one closest to what they
see).

## Options

1. **Whole-sketch replacement**: every edit is `EditFeature` with the full new
   sketch. Simple, but two collaborators editing different parts of a sketch
   always conflict, and logs grow with sketch size per edit.
2. **Fine-grained sketch ops** with stable item IDs, applied atomically in
   batches: add/delete entities and constraints, change values, update
   positions.
3. **Store only constraints, never positions**, and recompute geometry from
   scratch on replay. Breaks the "closest solution" behaviour, and any solver
   change would silently move every user's geometry.

## Decision

Option 2.

* `Feature::Sketch(SketchDef)` holds the plane, entities and constraints. The
  entity, constraint and geometry types are defined in
  `opendatum-sketch-solver` (pure data, serde) and instantiated with
  `SketchId` keys; `core-model` depends on that crate for the types only.
* `SketchId { op: OpId, index: u32 }`: the op that created the item and its
  index within that op's batch. In an edit, `op: OpId::NULL` means "the op
  being applied", so one batch can create a line and its end points and
  constrain them together. The reducer substitutes the real op id.
* `ModelOp::EditSketch { element, feature, edits: Vec<SketchEdit> }` applies
  its edits atomically (all or nothing) and is one undo step. Edits: add
  entity, add constraint, set a dimension value, delete items (cascading to
  dependent curves and constraints), set construction flags, and set point
  positions and radii.
* **Positions are stored.** They are the solver's starting point (the
  "seed"). A UI edit commits the change and the solved positions in one
  batch, so the stored geometry already satisfies the constraints.
  Regeneration solves again from the stored seed (a no-op for UI edits; it
  matters when upstream inputs change), so the solution is always the one
  nearest the geometry the user last saw.
* **Profiles from sketches.** Extrude and revolve take a `ProfileSource`:
  either inline curves (as today) or `Sketch { sketch: FeatureId, regions:
  Vec<RegionKey> }`. A `RegionKey` names a region of the sketch's planar
  arrangement by the sketch items on its boundary (design note 0006), with
  an interior point as a geometric tie-break, which is the sketch-level
  instance of ADR 0003's generative naming.

## Consequences

* Concurrent edits to different items of one sketch merge; edits to the same
  item follow the op log's last-writer-wins rule, and deleting an item
  another replica referenced makes that later constraint fail visibly
  (reported, never silently dropped).
* The serialized form of sketch types is part of the native format
  (docs/NATIVE_FORMAT.md) and versioned with it.
