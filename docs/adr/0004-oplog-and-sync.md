# ADR 0004 — The operation log is the document; undo is an operation

Status: accepted (Phase 0, implemented in `crates/oplog`)

## Context

Undo, full history, versions, branches, merge, real-time sync and offline
editing must all work on one structure, designed for sync from day one.

## Decision

* A document is an append-only set of `Op<P>` entries. `P` is the model
  payload (`ModelOp` in core-model).
* **IDs:** `OpId = (lamport: u64, replica: ReplicaId)`. A replica is one
  editing session on one device (random 64-bit id). IDs are never reused.
* **Causality:** each op lists the heads it was created against (`deps`).
  Ops are integrated only when their deps are present; early arrivals are
  buffered. Lamport time of an op exceeds its deps' (checked).
* **Replay order:** ascending `OpId`. Because Lamport order extends the
  causal order, this is a deterministic total order every replica agrees on.
* **Undo/redo:** `OpKind::Undo { target }` cancels its target; redo is an
  undo of that undo. An op is active iff no active undo targets it,
  computed in one descending pass (an undo always sorts after its target).
  Each replica keeps local undo/redo stacks of its own ops; they are
  reconstructible from the log (`rebuild_undo_stacks`).
* **History/versions/branches:** `effective_at(heads)` replays the causal
  closure of any set of heads. A version is a named, immutable head set; a
  branch is a mutable head set (Phase 5).
* **Sync:** `missing_for(remote_heads)` returns what a peer lacks;
  `integrate` accepts ops in any order. Convergence is property-tested.
* **Merge semantics** live in the reducer: concurrent ops are applied in
  `OpId` order; ops whose targets no longer exist become recorded conflicts
  (feature-level merge UI in Phase 10).

## Alternatives considered

* Snapshot + command stack (classic CAD undo): no sync story. Rejected.
* Full sequence CRDT for the feature list (RGA/Fugue): ordering of feature
  inserts uses "insert after feature X" anchored on stable feature IDs with
  OpId tiebreaks instead, which is sufficient for a tree edited at feature
  granularity. Revisit if fine-grained concurrent reordering demands it.
* Hash-chained (content-addressed) op IDs à la Automerge: better integrity,
  larger IDs. Can be layered later (sign/hash op batches) without changing
  the replay model.

## Consequences

Undo of an op that later ops depend on (e.g. undoing a feature another
replica edited) yields a well-defined state; the regeneration engine reports
the dependent feature as failed rather than corrupting anything.
