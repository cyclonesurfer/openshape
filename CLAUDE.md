# CLAUDE.md — working on OpenDatum

OpenDatum (named OpenShape until 2026-10-12; ADR 0017) is a native Linux,
pure-Rust, history-based parametric 3D CAD system (clean-room; Onshape is a
functional reference only — never copy its code, assets, icons, branding or
formats, and never use its name in the product). License: GPL-3.0-or-later.

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
cargo run -p opendatum      # desktop app (needs a GPU or Mesa)
cargo run -p opendatum -- --headless-demo out/   # offscreen proof, writes PNGs
cargo run -p opendatum -- --screenshots out/     # screenshots of the whole app
cargo run --release -p opendatum -- --bench out/ 600   # offscreen frame times
scripts/gpu-check.sh        # build, tests, demo, bench on a real GPU (docs/RUNNING.md)
cargo test --release -p opendatum-kernel --test timing -- --nocapture   # time guards
```

All four checks must pass before every commit. `main` is never broken.

The render tests and the headless demo need a Vulkan or GL driver. In a
container without a GPU, install Mesa's software Vulkan (lavapipe):
`apt-get install mesa-vulkan-drivers` — wgpu then reports `llvmpipe`.
`OPENDATUM_WGPU_BACKEND=vulkan|gl|any` forces a backend.

## Crate map

| crate | package | role |
|---|---|---|
| `crates/oplog` | `opendatum-oplog` | append-only op log: IDs, causal order, replay, undo/redo, sync primitives, CBOR container |
| `crates/core-model` | `opendatum-core-model` | documents, elements (tabs), feature tree as data, units, `ModelOp` reducer |
| `crates/nurbs` | `opendatum-nurbs` | NURBS curves/surfaces (Piegl & Tiller); **kernel-private** — only `kernel` may depend on it |
| `crates/kernel` | `opendatum-kernel` | B-rep kernel behind the `Kernel` trait |
| `crates/linalg` | `opendatum-linalg` | deterministic dense linear algebra (pivoted QR, COD, SVD) for the solvers (ADR 0010) |
| `crates/sketch-solver` | `opendatum-sketch-solver` | 2D sketch constraint solver with DOF/redundancy/conflict diagnostics (design note 0005) |
| `crates/naming` | `opendatum-naming` | persistent names of faces, edges and vertices from the kernel's origin tables, and resolution of stored references (design note 0013, ADR 0018) |
| `crates/features` | `opendatum-features` | evaluates a part studio's feature list through `Kernel` |
| `crates/render` | `opendatum-render` | wgpu viewport renderer + offscreen renderer; no kernel/egui deps |
| `crates/text` | `opendatum-text` | font layout for sketch text: bundled Liberation fonts (OFL), shaping, exact outlines (ADR 0015) |
| `crates/io` | `opendatum-io` | import/export formats (STL today) |
| `crates/ui-kit` | `opendatum-ui-kit` | the look: theme tokens, original icon set, view cube, offscreen egui rendering (ADR 0012) |
| `crates/ui` | `opendatum` (bin) | egui/eframe application shell (ADR 0012), offscreen screenshots and end-to-end UI tests |

Planned crates (not yet created; see docs/ARCHITECTURE.md): `regen`,
`mate-solver`, `script`, `drawing`, `sim`, `cam`,
`render-studio`, `server`, `client-sync`, `api`, `plugins`.

## Rules (from the project brief — non-negotiable)

- Linux only, by design (project decision, 2026-10-13): no Windows or
  macOS ports, CI or platform work. Others may port it under the GPL, so
  the cross-platform crates (wgpu, winit) stay, but no effort goes there.
- Rust stable, edition 2021, f64 throughout. No C/C++ geometry dependencies;
  kernel, solvers and file formats are pure Rust. FFI only for GPU drivers
  (wgpu), windowing, and `spacenavd`.
- Every dependency must be GPL-3.0-compatible; `deny.toml` enforces it.
- No kernel internals outside the `Kernel` trait. Bodies are opaque.
- Correctness over breadth: every kernel operation validates its output or
  fails with a clear error. Never silently corrupt a model.
- Efficiency is part of the bar (docs/PERFORMANCE.md): every user-facing
  operation has a written time budget enforced by a release-mode test
  (`tests/timing.rs`, failing at about ten times the measured time) and a
  scaling check; hot code uses standard structures (the kernel's BVH,
  grids, bulk loading, incremental updates), never loops over all pairs;
  work that grows with the input has a limit that refuses with a message;
  profile (callgrind) before recording numbers.
- Same input, same output, on every run and machine: no result may depend
  on a hash seed or on hash-table iteration order. The kernel uses
  `FxHashMap`/`FxHashSet` (`crates/kernel/clippy.toml` forbids the std
  tables there). A dependency with such a defect is used as a patched copy
  in `vendor/` (ADR 0022, `vendor/README.md`).
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
