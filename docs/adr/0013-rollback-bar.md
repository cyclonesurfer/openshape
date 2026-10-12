# ADR 0013 — The rollback bar in the op log

Status: accepted (Phase 2)

## Context

History-based modelers let the user roll the feature list back to an
earlier point: the features after a marker are kept but not regenerated, so
the model can be inspected as it was there and new features can be inserted
there. ADR 0012 puts that marker, the rollback bar, in the feature list. The
document is derived from the op log (ADR 0004), so the bar's position, and
how it reacts to every other edit, must be defined by the reducer:
deterministic on replay, convergent when collaborators edit concurrently,
and undoable.

## Options

1. **View state only**: the bar lives in the application, outside the log.
   Not undoable, lost on save, and a collaborator's insert cannot be placed
   relative to it.
2. **An index into the feature list**. Any concurrent insert or delete above
   it silently changes which features are active.
3. **An anchor feature in the log**: the bar is "after feature F", "at the
   start" or "at the end", set by its own op, and every op that changes the
   list adjusts it by fixed rules.

## Decision

Option 3. `PartStudio::rollback` is an `InsertAt`; `ModelOp::SetRollback`
moves it. Features `features[..rollback_index()]` are regenerated; the rest
report `FeatureStatus::RolledBack`.

Rules (in `Document::apply`):

* **Normal form.** After every op the bar is `End` when it is after the last
  feature (also when the list is empty), `Start` only when features follow
  it, otherwise `After(f)` for an existing `f`. So "at the end" stays at the
  end when features are appended, including a collaborator's.
* **Insert.** A feature inserted exactly at a bar that is not at the end
  goes above it: the bar advances past it. The application inserts new
  features at the bar, so the new feature is active and the bar follows it.
  An insert above the bar is active; an insert below it is rolled back.
* **Delete.** Deleting the anchor moves the bar to the feature before it (or
  the start), so the same features stay active.
* **Move.** As a delete followed by an insert: moving the anchor first moves
  the bar to the feature before it; a feature moved exactly to the bar goes
  above it.
* **SetRollback** to a feature that no longer exists is a recorded conflict;
  the bar does not move.

`crates/core-model/tests/rollback.rs` checks these rules against an
independent index-arithmetic model on random edit sequences, together with
normal form and replay equality.

Element deletion is in the same change: `ModelOp::DeleteElement` removes a
tab; deleting a document's last element is refused as a conflict, so
concurrent deletes of the last two tabs keep one.

## Consequences

* Moving the bar is an undoable step and is saved with the document.
* Two collaborators share one bar per Part Studio. Per-user bars would need
  per-replica view state in the log; not needed yet.
* Editing a feature below the bar requires moving the bar first; the
  application says so instead of regenerating hidden features.
* While a feature above the bar is being edited (its dialog or sketch mode
  is open), the application shows the history up to and including that
  feature and greys the features after it, as CAD users expect when they
  edit in place. This is view state (`Model::set_edit_point`), not an op:
  the document's bar does not move, nothing is recorded, and closing the
  editor shows the whole active history again. A toggle to preview the
  final result while editing is not implemented yet.
