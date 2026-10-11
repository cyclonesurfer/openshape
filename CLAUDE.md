# CLAUDE.md — working on OpenShape

OpenShape is a native Linux, pure-Rust, history-based parametric 3D CAD system
(clean-room; Onshape is a functional reference only — never copy its code,
assets, icons, branding or formats, and never use its name in the product).
License: GPL-3.0-or-later.

## Start of every session (and after any context compaction)

1. Read this file, `STATUS.md` and `docs/ROADMAP.md`.
2. Continue from STATUS.md's "Next steps".
3. Keep these files current as you work; they are the project's memory.

## Build, test, lint

```sh
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo deny check            # licenses, advisories, bans, sources
cargo run -p openshape      # desktop app (needs a GPU or Mesa)
cargo run -p openshape -- --headless-demo out/   # offscreen proof, writes PNGs
cargo run --release -p openshape -- --bench out/ 600   # offscreen frame times
scripts/gpu-check.sh        # build, tests, demo, bench on a real GPU (docs/RUNNING.md)
```

All four checks must pass before every commit. `main` is never broken.

The render tests and the headless demo need a Vulkan or GL driver. In a
container without a GPU, install Mesa's software Vulkan (lavapipe):
`apt-get install mesa-vulkan-drivers` — wgpu then reports `llvmpipe`.
`OPENSHAPE_WGPU_BACKEND=vulkan|gl|any` forces a backend.

## Crate map

| crate | package | role |
|---|---|---|
| `crates/oplog` | `openshape-oplog` | append-only op log: IDs, causal order, replay, undo/redo, sync primitives, CBOR container |
| `crates/core-model` | `openshape-core-model` | documents, elements (tabs), feature tree as data, units, `ModelOp` reducer |
| `crates/nurbs` | `openshape-nurbs` | NURBS curves/surfaces (Piegl & Tiller); **kernel-private** — only `kernel` may depend on it |
| `crates/kernel` | `openshape-kernel` | B-rep kernel behind the `Kernel` trait |
| `crates/linalg` | `openshape-linalg` | deterministic dense linear algebra (pivoted QR, COD, SVD) for the solvers (ADR 0010) |
| `crates/sketch-solver` | `openshape-sketch-solver` | 2D sketch constraint solver with DOF/redundancy/conflict diagnostics (design note 0005) |
| `crates/features` | `openshape-features` | evaluates a part studio's feature list through `Kernel` |
| `crates/render` | `openshape-render` | wgpu viewport renderer + offscreen renderer; no kernel/egui deps |
| `crates/io` | `openshape-io` | import/export formats (STL today) |
| `crates/ui-kit` | `openshape-ui-kit` | the look: theme tokens, original icon set, view cube, offscreen egui rendering (ADR 0012) |
| `crates/ui` | `openshape` (bin) | egui/eframe application shell |

Planned crates (not yet created; see docs/ARCHITECTURE.md): `naming`, `regen`,
`mate-solver`, `script`, `drawing`, `sim`, `cam`,
`render-studio`, `server`, `client-sync`, `api`, `plugins`.

## Rules (from the project brief — non-negotiable)

- Rust stable, edition 2021, f64 throughout. No C/C++ geometry dependencies;
  kernel, solvers and file formats are pure Rust. FFI only for GPU drivers
  (wgpu), windowing, and `spacenavd`.
- Every dependency must be GPL-3.0-compatible; `deny.toml` enforces it.
- No kernel internals outside the `Kernel` trait. Bodies are opaque.
- Correctness over breadth: every kernel operation validates its output or
  fails with a clear error. Never silently corrupt a model.
- Never stub a feature and report it done. ROADMAP items are ticked only when
  implemented **and tested**. If something can't be done (e.g. a proprietary
  format without a public spec) say so in STATUS.md / docs/FORMATS.md and
  implement the best legitimate alternative.
- Write a design note in `docs/design/` (citing literature) before hard
  algorithms: persistent naming, NURBS surface intersection, booleans,
  fillets, offsets, constraint solving, merge.
- Ambiguity → write an ADR (`docs/adr/NNNN-title.md`) with options and a
  recommendation, then proceed with the recommendation. Stop and ask only for
  decisions that are irreversible and evenly balanced.
- Document every public API item (`missing_docs` is on).
- `unsafe` is forbidden workspace-wide.
- Solo mode must never link the server crate (it will live behind a feature).
- Clean room (docs/CLEANROOM.md): build only from public standards,
  literature, public documentation and our own reasoning. Never consult or
  use another product's code, binaries, leaked material, proprietary
  formats, icons, colours, screenshots or name; every icon is drawn in code
  by us. The product never names the reference product
  (`crates/ui/tests/cleanroom.rs`). The project record may become legal
  evidence: describe your own decisions and sources.
- Golden-model tests assert volume, area, topology counts and validity.
  Fuzz solvers, booleans, parsers and op-log merge.
- Commit at every green milestone with a clear message.

## Conventions

- Internal units: metres and radians. UI displays millimetres by default.
- Model tolerance: `Tolerance::default()` = 1e-7 m linear, 1e-9 rad angular.
- Periodic faces never span more than π (ADR 0008).
- Feature IDs and element IDs derive from the `OpId` that created them.
- Tests: unit tests next to code; integration/golden tests in `crates/*/tests/`.
- Use subagents for independent crates with a precise interface contract and
  acceptance tests; review their output before merging.
