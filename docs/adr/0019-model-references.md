# ADR 0019 — The first references to model geometry

Status: accepted, 2026-10-12 (Phase 3; design note 0013 §10 step 3).

## Context

Persistent names exist for box, extrude and revolve bodies (design note
0013 §11, ADR 0018). Step 3 of the note makes features use them: a
sketch placed on a planar face and an extrude that ends at a face, both
*carrier* references (they need only the plane the face lies on), with
picking of faces, edges and vertices in the viewport. Five questions
need answers before code:

1. how the references are stored in feature definitions;
2. which coordinate frame a sketch on a face gets;
3. what "up to a face" means before the kernel can intersect surfaces;
4. how a feature learns the sketch edge behind each profile curve;
5. what a body whose names cannot be made does to its feature;
6. how a stored hint is checked against the model.

## Options and decisions

### 1. Storage

* **A. New variants and optional fields.** `PlaneRef::Face`, holding a
  `Reference`, for a sketch's plane; `Feature::Extrude` gains `up_to`, an
  optional `Reference` (serde default `None`); both are boxed, as names
  are large. Documents saved before this change read as they were: a
  missing field is `None`, and old documents have no `Face` plane.
* **B. A new `Extent` enum replacing `distance` and `symmetric`.** Cleaner
  for the end conditions still to come (up to vertex, through all, up to
  next), but every stored extrude changes shape and needs a migration.

**Decision: A.** When `up_to` is set the stored distance is not used (the
face gives it) and `symmetric` must be false; the extrude dialog keeps the
last typed distance so switching back restores it. The end-condition enum
of B can replace both fields with a format migration when the third end
condition arrives.

### 2. The frame of a sketch on a face

The face must be planar; a sketch on a curved face fails with a message
naming it. The frame's normal is the face's outward normal, so the sketch
is seen from outside the part. For its origin and x axis:

* **A. The kernel surface's own frame.** It depends on how the operation
  built the plane (a profile's start point, a line's start), which the
  user never sees, and it moves sideways when those change.
* **B. From the plane alone**: the origin is the world origin projected
  onto the plane; the x axis is world X projected for a face that faces
  mostly up or down (|n·Z| at least |n·X| and |n·Y|), and Z × n for one
  that faces mostly sideways, which keeps the sketch upright (its y axis
  points up). On faces parallel to the standard planes this gives exactly
  the Top, Front and Right frames (and their opposites), so a sketch on a
  box's front face has the coordinates of a sketch on the Front plane,
  offset. A face that moves along its normal moves the sketch with it
  without shifting it sideways.

**Decision: B.** Like any rule that picks a direction on every plane, it
must change abruptly somewhere: here where a face tilts through 45°
between facing up and facing sideways. The sketch then turns with the
rule; the solver keeps its constraints, so only unconstrained geometry
moves. This is recorded as a known limit.

### 3. Up to a face before surface intersection

Ending an extrude exactly on an arbitrary face needs the side faces cut
by that face's surface: curve/surface intersection (design note 0011)
and, in general, a boolean. Until those exist:

* **A. Planes parallel to the sketch only**: the extrude distance is the
  signed distance from the sketch plane to the face's plane, measured on
  the sketch normal, and the end cap lies on the face's plane exactly. A
  face that is not a plane parallel to the sketch (within the angular
  tolerance) fails with a message saying so; so does a face at zero
  distance.
* **B. Approximate other faces** by the nearest parallel plane through the
  face's farthest point. Silently wrong geometry; rejected.

**Decision: A**, documented as a limit, the ROADMAP item unticked until
other faces work.

### 4. Profile curves and sketch edges

`Kernel::region_profile_edges` (beside `region_profile`, which is
unchanged) returns, with each curve, the arrangement edge it was made
from and the direction the profile runs along it (a closed edge cut in
two gives two consecutive curves of one edge). The features crate turns
those into the `EdgeKey`s the naming layer needs, with the same function
that names regions (design note 0006 §5), so a cap's region keys and its
sides' edge keys agree by construction. Inline profiles (curves given
directly: the application's cylinder, plate and ring templates, and
tests) get made-up keys: curve `i` is an edge of its own entity, numbered
`i` in the feature's op, and their regions are the profile regions the
kernel found, keyed by those edges. Their names last as long as the curve
order does, which only the template dialogs set.

### 5. A body that cannot be named

The geometry is valid; only references to it are affected. The feature
**succeeds with a warning** that its faces cannot be referenced, giving
the naming error, and any reference to the body then fails with that
message. Failing the feature would hide valid geometry; silently skipping
names would break references without a cause.

### 6. Checking a hint

Resolution asks whether the stored hint (the picked point) lies on a
candidate entity: to break the tie between entities only an ordinal tells
apart, to choose the piece a single reference needs, and to report that
an entity moved (design note 0013 §6). The kernel has no exact
point-to-face distance yet.

* **A. Distance to the entity's mesh** (chord tolerance 20 µm, as fine as
  the display mesh picks are made on), within 0.1 mm: both meshes'
  deviation with room to spare.
* **B. Project onto the exact surface and test the trimming loops.** Exact,
  but a new kernel query of its own (point classification against
  trimmed NURBS faces), for a tie-break.

**Decision: A**, as `HINT_TOLERANCE`. A hint is only a tie-break: one within
the tolerance of two candidates is ambiguous, an error, never a guess. The
known limit is models with features finer than 0.1 mm, where ties between
pieces closer than that fail as ambiguous. B replaces it when the kernel
gains point classification (needed for booleans anyway, design note 0012).

## Consequences

* `StudioResult` keeps every body's names; references resolve against the
  bodies made before the referring feature (each feature still makes its
  own part until booleans), first by exact name, then by the descent rules
  of design note 0013 §6. A reference to a part that is missing,
  suppressed, failed or later in the list, to an entity that no longer
  exists, to a face that is not planar, or (up to) to a face not parallel
  to the sketch or in its plane fails the feature with a message saying
  which; a face split since it was picked is used when all its pieces lie
  in one plane. A symmetric extrude cannot end at a face.
* The kernel gains `Kernel::face_plane` (the plane and outward normal of
  a face whose surface is a plane: a plane, the extrusion of a line, or
  the revolution of a line perpendicular to the axis) and
  `Kernel::region_profile_edges`.
* The viewport picks faces (ray against the mesh), edges (screen distance
  to the polylines) and vertices (screen distance to the points), the
  nearest kind within a few pixels winning in that order: vertex, edge,
  face. A pick becomes a stored `Reference` with its name and the hit
  point and normal or tangent as the hint.
* Tests: references to box, sketch-extrude and inline-profile faces
  (caps, sides, a face in two pieces), through dimension changes, a
  region gaining a hole, a sketch edge redrawn, suppression, reordering
  and deletion, checked against the parts' positions; the frame rule on
  200 tilted planes; old documents regenerate unchanged (an op-log
  fixture saved by the code before this ADR).
