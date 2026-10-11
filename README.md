# OpenShape

A native Linux, pure-Rust, history-based parametric 3D CAD system.
License: GPL-3.0-or-later.

OpenShape is early: Phase 0 (foundations), Phase 1 (kernel core) and the
first Phase 2 milestone (the sketcher: constraint solver, sketch mode,
dimensions, region extrude and revolve) are in place. See [STATUS.md](STATUS.md) for what works today and
[docs/ROADMAP.md](docs/ROADMAP.md) for the plan.

## Build and run (Linux)

Requirements: a stable Rust toolchain (rustup), and Vulkan drivers (Mesa
RADV on AMD, `mesa-vulkan-drivers` on Debian/Ubuntu, `vulkan-radeon` on
Arch, `mesa-vulkan-drivers` on Fedora).

```sh
cargo run --release -p openshape            # the desktop app
cargo run --release -p openshape -- --headless-demo out/   # offscreen demo
cargo run --release -p openshape -- --bench out/            # frame times
cargo test --workspace
```

[docs/RUNNING.md](docs/RUNNING.md) has the full guide: prerequisites per
distribution, the desktop app, developer checks, the offscreen demo and
benchmark, and `scripts/gpu-check.sh` for checking a GPU machine (including a
headless server) end to end.

In the app: Insert menu to add features; drag with the left mouse button to
orbit, right or middle button to pan, scroll to zoom, F to fit; Ctrl+Z and
Ctrl+Shift+Z to undo and redo; right-click a feature for edit, suppress,
reorder and delete; File menu to save/open documents and export STL.

## Layout

See [CLAUDE.md](CLAUDE.md) for the crate map and contributor rules, and
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the design.
