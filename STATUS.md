# OpenShape status

Updated 2026-10-11, Phase 2 (sketcher): core sketcher done; editing tools
and the new application shell in progress.

## What works

- **Op log** (`openshape-oplog`): Lamport/replica op IDs with causal deps,
  deterministic replay, undo/redo stored as ops (redo = undo of undo),
  history at any set of heads, sync primitives (`missing_for`, `integrate`
  with a pending buffer for out-of-order ops), CBOR container. Property tests
  check undo/redo against a reference model and replica convergence.
- **Document model** (`openshape-core-model`): documents with element tabs, a
  Part Studio feature list as data (insert, edit, rename, suppress, move,
  delete), a reducer over the op log, units parsing of single quantities
  (`"25.4 mm"`, `"1 in"`, `"90 deg"`; bare numbers use the default unit).
  Expressions and variables come with Variable Studios.
- **NURBS** (`openshape-nurbs`, kernel-private): rational curves and surfaces
  after Piegl & Tiller — evaluation, derivatives, knot insertion, split,
  closest point, exact arcs, interpolation, extrusion/revolution surfaces.
- **Kernel** (`openshape-kernel`, behind the `Kernel` trait):
  - B-rep arena (solid/shell/face/loop/coedge/edge/vertex) with validation:
    closed loops, manifold edge use, edges on their surfaces, Euler–Poincaré,
    outward orientation. Every operation validates its output.
  - Box, extrude (blind, symmetric; multi-region profiles with holes; lines,
    arcs, circles, splines) and revolve (full and partial, profiles touching
    the axis, offset axes). Periodic faces split at ≤ π (ADR 0008).
  - Surface classification: plane, cylinder, cone, sphere, torus, extrusion,
    revolution, NURBS.
  - Watertight tessellation with chord and angle tolerances (shared edge
    discretisation, constrained Delaunay per face, singular points handled).
  - Exact mass properties (volume, area, centroid, inertia) by boundary
    integration, accurate to 1e-9 relative on the golden corpus.
- **Features** (`openshape-features`): regenerates a Part Studio through the
  kernel with per-feature status (ok, warning, suppressed, failed with
  message). A failing feature never blocks the others. Sketch features are
  solved and their regions computed during regeneration; extrude and revolve
  take their profile from picked sketch regions, which are named
  topologically and re-found after sketch edits (design note 0006 §5).
- **Sketch solver** (`openshape-sketch-solver`, on the deterministic
  `openshape-linalg`): points, lines, circles, arcs, through-point and
  control-point splines under 15 geometric constraint kinds and 5
  dimension kinds. Weighted
  minimum-norm Gauss–Newton with Levenberg–Marquardt fallback and
  continuation (design note 0005, §9 lists what changed while building it).
  Reports degrees of freedom, free entities, redundant and conflicting
  constraints, and never reports success on a result that fails an
  independent final check. 50-entity sketches solve in under 1 ms and a
  300-entity chain in 11.6 ms (release).
- **Sketch regions** (`openshape-kernel`, `arrangement`): the planar
  arrangement of arbitrary sketch curves (lines, arcs, circles, rational
  splines; crossings, tangencies, overlaps, duplicates, self-intersections)
  with exact areas, holes, and region profiles ready for extrude and
  revolve.
- **Rendering** (`openshape-render`): wgpu renderer for the window and an
  offscreen renderer; shaded parts with edges, perspective/orthographic
  camera, orbit/pan/zoom, view presets, fit. Vulkan first, GL fallback.
- **App** (`openshape` binary, egui/eframe): File (new, open/save `.osdoc`,
  export STL), Edit (undo/redo, Ctrl+Z / Ctrl+Shift+Z / Ctrl+Y), Insert
  (sketch on the top, front or right plane; extrude and revolve from picked
  sketch regions; templates), View (presets, fit with F, perspective/ortho,
  show sketches). Feature tree with status icons; right-click to edit,
  rename, suppress, move up/down, delete. Properties panel shows mass
  properties. Dialogs take values with units.
- **Sketch mode**: line (chained), corner rectangle, circle, 3-point arc,
  centre arc, point, through-point spline, construction toggle; 12
  constraint tools and dimensions (length, distance aligned/horizontal/
  vertical, angle, radius, diameter) typed with units; snapping to points,
  curves and the origin with horizontal/vertical inference; dragging with
  live solving; box selection and delete; colour-coded state (free, fully
  constrained, redundant, conflicting) with constraint glyphs and dimension
  labels. Every edit is one undoable op. The editor is UI-toolkit
  independent and tested end to end by driving it with screen clicks.
- **UI kit** (`openshape-ui-kit`): the application's look in one place
  (ADR 0012): theme tokens, 60 original icons drawn in code, a view cube,
  and an offscreen egui renderer used for screenshots and UI tests. Not yet
  wired into the application shell.
- **I/O** (`openshape-io`): STL export, ASCII and binary.
- **Headless demo**: `cargo run -p openshape -- --headless-demo out/` builds
  a box through the op log, orbits, adds three more parts, undoes and redoes,
  saves and reopens, exports STL, and checks each step pixel-exactly. All 12
  checks pass on Mesa lavapipe.
- **Benchmark**: `openshape --bench <dir> [frames]` times offscreen frames of
  the demo scene (36,656 triangles, 1920x1080, 4x MSAA) to GPU completion.
  Lavapipe in the cloud container: about 67 ms mean (15 fps).
  `scripts/gpu-check.sh` runs build, tests, demo and benchmark on a real GPU
  and collects a report (docs/RUNNING.md).

## Verified so far

- `cargo fmt`, `cargo clippy -D warnings`, `cargo test --workspace` and
  `cargo deny check` pass.
- Real GPU, 2026-10-10, AMD Radeon 8060S (RADV, Mesa 26.0.6, Fedora 44),
  commit f844d3b: the headless demo passed all 12 checks; 156 of 158 tests
  passed. The 2 failures were render tests sampling edge lines: hardware
  1-pixel lines under MSAA came out faint on RADV (they were solid on
  lavapipe). Edges are now drawn as fixed-width screen-space quads (1.5 px
  by default, scaled by DPI in the app), so line weight no longer depends on
  the driver, and the two tests locate lines from geometry instead of one
  pixel or a colour threshold.
- Re-run on the same machine, commit d62aee9: build, demo (12/12),
  benchmark and render timing pass; 157 of 158 tests passed. 1080p
  benchmark: 0.120 ms mean, 0.230 ms p99 (8,350 fps) for the 4-part demo
  scene; 1,000 instances: 0.5 ms per frame steady, including readback.
  The one failure was a test assumption, not a rendering bug: it required
  one pixel to be bit-identical between the GPU's sRGB encoder and the
  shader's, which Vulkan does not guarantee (lavapipe also differs by one
  code on 9,184 pixels; RADV happened to differ at the sampled one). The
  test now bounds the difference at one code on every pixel away from
  antialiased edges and across full-range and dark gradients, which is
  stricter than before, and fails if the encoder is replaced by a plain
  gamma 2.2 curve. Confirming run on RADV pending.
- The windowed app has not yet been run on real hardware.

## Not yet / known limits

- Sketch editing tools are not done: trim, extend, offset, sketch fillet,
  mirror, patterns, slots, polygons, ellipses, conics, tangent arcs, centre
  and aligned rectangles, text, use/project. Parallel, perpendicular and
  tangent inference while drawing is not done.
- The application shell still has the Phase 1 layout; the new layout
  (ADR 0012: feature list with rollback bar, feature toolbar, view cube,
  element tabs, command search) is next.
- Sketches can only be placed on the three standard planes until Phase 3
  brings faces and construction planes.
- Extrude rejects some valid but extreme sketch regions with an error (never
  with wrong geometry): in the random arrangement stress tests, 92 region
  extrudes over about 1,900 random sketches failed in the kernel's
  tessellation or profile checks (slivers and narrow notches between curves,
  tolerance-scale cusps), because tessellation validates with a chord
  tolerance of 1% of the body size and the profile builder checks
  self-intersection on 48-sample polylines. Realistic cases probed all
  extrude exactly (volume within 1.2e-9 relative): walls down to 0.01 mm next to a
  hole, cusp regions next to circles tangent to two sides, 0.2 mm slots,
  0.01 mm slivers, V notches. Making both checks adapt to local feature
  size is on the kernel list.
- Solver conflict detection is local: from a start far from any solution it
  can call a sketch conflicting although another branch would satisfy it.
- No booleans: each feature makes its own part (Phase 3).
- No persistent naming or incremental regeneration: every edit regenerates
  the whole studio (Phase 3). Features don't reference model geometry yet, so
  nothing can break on upstream edits today.
- Non-manifold topology deferred to Phase 3.
- File dialogs are a text path prompt, not a native picker.
- No STEP or other exchange formats beyond STL (Phase 4; see docs/FORMATS.md).

## Next steps

1. Phase 2: design note for the constraint solver, then the sketch data model
   in core-model, the solver crate (Newton/LM, DOF analysis, over/under-
   constraint diagnostics, property tests) and region detection feeding
   extrude/revolve.
2. Re-run `scripts/gpu-check.sh` on RADV to confirm the encoder-tolerance
   change (docs/huginn-tasks.md has the prompt), and run the windowed app on
   a desktop.
3. Phase 3 groundwork: design notes for curve/surface intersection, booleans
   and persistent naming, so the sketcher's references are built on them.
