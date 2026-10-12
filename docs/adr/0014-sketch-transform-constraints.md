# ADR 0014 — Constraints under sketch transforms

Status: accepted (Phase 2)

## Context

The sketcher's Move, Rotate and Scale tools (design note 0007 §13) apply a
similarity T(p) = c + s·R(θ)(p − c) + v to the selection. The sketch is
constraint-driven, so the tools must say what happens to every constraint
the selection takes part in. A similarity keeps every relation between the
selected entities (coincidence, tangency, parallelism, equal lengths…) and
every angle, and multiplies every length by s. Three kinds of constraint
are not kept by construction:

1. constraints to the **sketch frame**: horizontal, vertical, horizontal and
   vertical distances, fixes (a rotation breaks the first three; any motion
   moves a fixed point);
2. **boundary relations** between selected and unselected geometry;
3. **boundary dimensions** between selected and unselected geometry.

The choices are ours (the clean-room rule: behaviour from public
documentation and our own reasoning). Publicly documented sketchers offer a
move/rotate/scale command with an option to keep or break relations to
unselected geometry; the details below are decided here.

## Options

Boundary relations:

* **A. Keep them**, solving like a drag: the rest of the sketch moves to
  follow (a line attached to a moved rectangle stretches). The selection
  may then be unable to reach the requested place.
* **B. Remove them** (detach): always reaches the place, but silently opens
  closed profiles if it is the only behaviour.
* **C. Refuse** any transform that breaks one. Safe but rarely useful.

Boundary dimensions:

* **A. Measure them anew** after the transform (the moved hole's distance to
  the edge becomes the new distance).
* **B. Keep driving** the unselected geometry (the edge moves with the
  hole). Surprising: dimensions would move things the user did not select.

Frame constraints under a rotation that is not a whole number of quarter
turns:

* **A. Remove them all**: a rotated rectangle drawn with horizontal and
  vertical lines loses its right angles.
* **B. Convert them**: keep what the rotation keeps — lines that were
  horizontal or vertical stay parallel or perpendicular to one another; one
  of them loses its constraint, which frees exactly the rotation.
* **C. Refuse** such rotations.

Fixed points in the selection:

* **A. The fix moves with its point** (its position is transformed, like a
  dimension measured anew).
* **B. Refuse** to move a fixed point. **C. Remove** the fix.

## Decision

* Boundary relations: **A by default, B as an option** ("Keep constraints"
  in the sketch dialog). With A, a transform the rest of the sketch cannot
  follow (a point attached to a fixed one) is refused with the reason,
  never approximated: the selection lands exactly where asked or nothing
  happens.
* Boundary dimensions: **A** in both modes; a dimension that would measure
  nothing (≤ 1e-7 m) is removed and reported.
* Frame constraints: **B**. Quarter turns swap horizontal and vertical
  exactly; other angles convert lines as above. Horizontal and vertical
  point pairs, and horizontal or vertical distances between points that
  were not level, have no equivalent in the rotated frame: they are removed
  and reported (principle 4 of design note 0007: no silent loss of intent).
* Fixes: **A**. The transform is an explicit instruction about where the
  point goes; the fix keeps the selection as constrained as it was.
* Copies (the Copy option) get the constraints among the copied entities,
  adjusted the same way, and nothing tying them to the rest of the sketch.

## Consequences

* Every outcome is exact or refused: the moved geometry is pinned where T
  puts it while the rest of the sketch solves (each moved point lands at
  T(p) to rounding, each moved circle has radius s·r), and the tools commit
  through the checked edit (no new redundancy, the sketch still solves).
  A transform that would collapse a curve joined to the selection is
  refused too.
* Degrees of freedom change only by what the status message reports
  (removed constraints, the freed rotation).
* Revisit if users want boundary dimensions to drive (option B) for some
  workflows; it could become a third mode without changing the others.
