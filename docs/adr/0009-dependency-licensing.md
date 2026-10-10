# ADR 0009 — Dependency licensing policy

Status: accepted (Phase 0)

OpenShape is GPL-3.0-or-later. Every dependency must be compatible with
GPL-3.0. `deny.toml` allows: MIT, Apache-2.0, Apache-2.0 WITH LLVM-exception,
BSD-2-Clause, BSD-3-Clause, ISC, Zlib, 0BSD, Unicode-3.0, Unicode-DFS-2016,
BSL-1.0, CC0-1.0, MPL-2.0, LGPL-2.1-or-later, LGPL-3.0-or-later,
GPL-3.0-or-later and GPL-3.0-only. Anything else fails CI (`cargo deny check`).
Unknown registries and git sources are denied. Advisories are checked.
Geometry, solver and file-format code must be pure Rust; FFI is allowed only
for GPU drivers (wgpu/ash), windowing (winit and its system libraries) and
`spacenavd`.
