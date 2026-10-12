# OpenDatum

A native Linux, pure-Rust, history-based parametric 3D CAD system.
License: GPL-3.0-or-later.

OpenDatum is early: Phase 0 (foundations), Phase 1 (kernel core) and the
first Phase 2 milestones (the sketcher: constraint solver, sketch mode,
dimensions, region extrude and revolve; and the application shell: feature
list with rollback bar, toolbars, view cube, floating dialogs, tabs, command
search) are in place. See [STATUS.md](STATUS.md) for what works today and
[docs/ROADMAP.md](docs/ROADMAP.md) for the plan.

## Build and run (Linux)

Requirements: a stable Rust toolchain (rustup), and Vulkan drivers (Mesa
RADV on AMD, `mesa-vulkan-drivers` on Debian/Ubuntu, `vulkan-radeon` on
Arch, `mesa-vulkan-drivers` on Fedora).

```sh
cargo run --release -p opendatum            # the desktop app
cargo run --release -p opendatum -- --headless-demo out/   # offscreen demo
cargo run --release -p opendatum -- --bench out/            # frame times
cargo test --workspace
```

[docs/RUNNING.md](docs/RUNNING.md) has the full guide: prerequisites per
distribution, the desktop app, developer checks, the offscreen demo and
benchmark, and `scripts/gpu-check.sh` for checking a GPU machine (including a
headless server) end to end.

In the app: the toolbar adds sketches and features (S sketch, E extrude,
V revolve) and Ctrl+K searches every command; right-drag orbits, middle-drag
pans, the wheel zooms, F fits; Ctrl+Z and Ctrl+Shift+Z undo and redo;
right-click a feature to edit, rename, suppress, move or delete it, and drag
the rollback bar to roll the history back; the menu at the top left saves,
opens and exports STL.

## Layout

See [CLAUDE.md](CLAUDE.md) for the crate map and contributor rules, and
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the design.
