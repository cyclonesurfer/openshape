# ADR 0005 — Native document format: versioned CBOR op log

Status: accepted (Phase 0 draft; file I/O and attachments in Phase 5)

## Decision

The native format **is the op log**: a CBOR (RFC 8949) map
`{ magic: "OpenDatum-oplog", version: 1, ops: [Op, …] }` with ops in replay
order, written by `OpLog::encode`. File extension `.odatum`. (ADR 0017
renamed them from `"OpenShape-oplog"` and `.osdoc`; old documents still
open.) Payload schemas are the serde representations of
`core-model::ModelOp`, documented in `docs/NATIVE_FORMAT.md` and versioned
with the container. Large attachments (imported files, images,
thumbnails) will live alongside the log in a zip container in Phase 5,
referenced by content hash.

## Why

* No derived state is persisted, so files cannot disagree with their history.
* CBOR is compact, self-describing, has a public standard and many
  implementations, keeping the format open.
* Same bytes are used for sync, so there is one serializer to test and fuzz.

## Rejected

JSON (too large for long histories, but available for debugging via serde);
a custom binary format (needless spec burden); SQLite (C dependency).
