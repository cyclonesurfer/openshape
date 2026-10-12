# OpenDatum architecture

OpenDatum is a history-based parametric CAD system. Its single source of
truth is a document's **operation log**; everything the user sees is derived
from it by deterministic replay and regeneration.

```
 user action ──► ModelOp ──► OpLog (append-only, causal, replicated)
                                 │  effective ops in OpId order
                                 ▼
                         core-model reducer ──► Document (elements, feature trees)
                                                     │  per Part Studio
                                                     ▼
                             features (regen) ──► Kernel trait ──► Bodies (opaque B-rep)
                                                     │
                                       tessellate ◄──┘
                                           │
                       render (wgpu) ◄─────┘──► io (STL, …)
                           │
                         ui (egui shell)
```

## Layers and crates

| Layer | Crate | Status |
|---|---|---|
| Op log, sync primitives | `oplog` | Phase 0 ✔ (undo/redo, history, CBOR container, convergence tests) |
| Document model | `core-model` | Phase 0 ✔ (part studios, primitive/extrude/revolve features, units) |
| Geometry kernel | `kernel` (+ private `nurbs`) | Phase 1 ✔ for B-rep, NURBS, tessellation, extrude, revolve |
| Persistent naming | `naming` | Phase 3: names of box, extrude and revolve entities from the kernel's origin tables, and reference resolution ✔ (design note 0013, ADR 0018); references from features and booleans' splits and merges next |
| Regeneration | `features` today; `regen` in Phase 3 | Phase 0 sequential evaluator ✔ |
| Sketcher | `sketch-solver` | Phase 2 |
| Assemblies | `mate-solver` | Phase 6 |
| Custom features | `script` | Phase 9 |
| Viewport | `render` | Phase 0 ✔ (shaded + edges, orbit camera, offscreen) |
| Application shell | `ui` (binary `opendatum`) | Phase 0 ✔ (egui, feature tree, undo/redo) |
| Formats | `io` | STL export ✔; rest per ROADMAP |
| Drawings | `drawing` | Phase 7 |
| Simulation | `sim` | Phase 11 |
| CAM | `cam` | Phase 12 |
| Path tracer | `render-studio` | Phase 13 |
| Collaboration | `server`, `client-sync` | Phase 10 |
| Automation API | `api` | Phase 9 |
| Plugins | `plugins` (wasmtime) | Phase 9 |

## Key decisions

* **One data structure for history.** Undo, history restore, versions,
  branches, merge and sync all operate on the op log (ADR 0004). Undo is an
  appended operation, never a deletion, so it composes with collaboration.
* **Deterministic replay.** `OpId = (lamport, replica)` is a total order that
  extends causality. The reducer in `core-model` is a pure function, so every
  replica that holds the same op set computes the same document.
* **Stable IDs from op IDs.** Elements and features are identified by the
  `OpId` that created them, so IDs are unique across replicas without
  coordination.
* **Kernel isolation.** The rest of the system sees only the `Kernel` trait,
  neutral inputs (`Profile`, `Frame`, specs) and neutral outputs (`BodyMesh`,
  `MassProperties`, `TopologyCounts`). Bodies are opaque (ADR 0001).
* **Validated geometry.** Every kernel operation runs the B-rep validator
  (closed loops, manifold edge use, edges on surfaces, Euler–Poincaré,
  outward orientation) before returning; failures become feature errors.
* **Renderer independence.** `render` consumes neutral `MeshData` and owns
  its wgpu pipelines. The UI renders the viewport into a texture shown by
  egui, and tests render the same scene offscreen (ADR 0007).

## Deployment modes

1. **Solo desktop** (`opendatum` binary): local, offline; documents are op-log
   files on disk (ADR 0005, `docs/NATIVE_FORMAT.md`). Never links the server.
2. **Client/server** (Phase 10): `opendatum-server` hosts the same op logs;
   clients exchange ops using `OpLog::missing_for` / `integrate`, which already
   handle out-of-order delivery and concurrent edits.

## Kernel internals (private)

* `math` — frames, tolerances, rotations.
* `geom` — `Curve3` (line, circle, NURBS) and `Surface` (plane, extrusion,
  revolution, NURBS) with evaluation, derivatives, projection, classification.
* `topo` — arena B-rep: vertex, edge, coedge, loop, face, shell, solid.
* `profile` — profile curves → oriented loops → regions (outer + holes),
  with robust-predicate intersection checks.
* `sweep` — extrude and revolve building B-reps directly.
* `tess` — watertight, tolerance-driven tessellation (design note 0002).
* `check` — validation, Euler–Poincaré counts, mass properties.
