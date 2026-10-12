# OpenDatum native document format (draft, container version 1)

Status: draft. The container is implemented (`OpLog::encode`/`decode`) and
used by the desktop app's Save/Open; payload schemas will be frozen and
versioned when Phase 5 (documents and history) completes.

## Principles

* The file **is the operation log**. No derived geometry is stored; opening
  a document replays the log through the reducer and regenerates.
* Open and documented: CBOR (RFC 8949) with the structures below, so any
  language can read and write documents.
* Identical bytes are used for sync between replicas.

## Container

File extension `.odatum` (documents saved as `.osdoc` before the rename of
ADR 0017 open too: opening does not look at the extension). One CBOR map:

| key | type | meaning |
|---|---|---|
| `magic` | text | `"OpenDatum-oplog"`; readers also accept `"OpenShape-oplog"`, written before the rename (ADR 0017) |
| `version` | uint | container version, currently `1` |
| `ops` | array of `Op` | every operation, in replay order (ascending `OpId`) |

`Op` = map `{ id: OpId, deps: [OpId], kind: OpKind }`

* `OpId` = map `{ lamport: uint, replica: uint }` (replica is a 64-bit id;
  `0` is reserved for genesis ops such as a new document's first Part Studio).
* `deps` lists the heads of the log when the op was created. Every dep must
  appear in the file; `lamport` must exceed every dep's `lamport`.
* `OpKind` is the serde externally-tagged enum:
  `{ "Apply": ModelOp }` or `{ "Undo": { "target": OpId } }`.

An op is **active** unless an active `Undo` targets it. The document is
obtained by applying the payloads of active `Apply` ops in ascending `OpId`
order.

## ModelOp (payload, schema 1)

Externally tagged serde enums; lengths in metres, angles in radians.

* `CreateElement { name, kind: "PartStudio" }` — element id = this op's id.
* `RenameElement { element, name }`
* `InsertFeature { element, at, name, feature }` — feature id = this op's id;
  `at` is `"End"`, `"Start"` or `{ "After": FeatureId }`.
* `EditFeature { element, feature, definition }`
* `RenameFeature { element, feature, name }`
* `SetSuppressed { element, feature, suppressed }`
* `MoveFeature { element, feature, to }`
* `DeleteFeature { element, feature }`

`ElementId` and `FeatureId` are newtypes over `OpId`.

`Feature` is one of:

* `Box { corner: [x,y,z], size: [x,y,z] }`
* `Extrude { profile: ProfileDef, distance, symmetric }`
* `Revolve { profile: ProfileDef, axis_origin: [x,y], axis_dir: [x,y], angle }`

`ProfileDef = { plane: PlaneRef, curves: [CurveDef] }` where `PlaneRef` is
`"Top" | "Front" | "Right" | { Custom: { origin, normal, x_dir } }` and
`CurveDef` is `Line { start, end }`, `Arc { center, radius, start_angle,
sweep_angle }`, `Circle { center, radius }` or `Spline { through: [[x,y]…] }`.

Ops that no longer apply on replay (for example an edit of a feature deleted
concurrently) are skipped and reported as conflicts; they never make a file
unreadable.

## Planned (Phase 5)

* Zip container holding `log.cbor`, content-addressed attachments
  (imported files, images), thumbnails, and a manifest with schema versions.
* Named versions and branches as head sets stored in the log.
