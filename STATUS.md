# OpenShape status

Updated 2026-10-12, Phase 2 (sketcher): core sketcher and the new
application shell done; sketch editing tools next.

## What works

- **Op log** (`openshape-oplog`): Lamport/replica op IDs with causal deps,
  deterministic replay, undo/redo stored as ops (redo = undo of undo),
  history at any set of heads, sync primitives (`missing_for`, `integrate`
  with a pending buffer for out-of-order ops), CBOR container. Property tests
  check undo/redo against a reference model and replica convergence.
- **Document model** (`openshape-core-model`): documents with element tabs, a
  Part Studio feature list as data (insert, edit, rename, suppress, move,
  delete), a rollback bar stored in the op log (ADR 0013; its rules are
  property-tested against an independent model), element tabs that can be
  added, renamed and deleted, a reducer over the op log, units parsing of single quantities
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
    It marks seam edges, where two pieces of one surface meet (the splits of
    ADR 0008), so the display draws only creases; the comparison is
    geometric and errs towards showing an edge.
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
  camera, orbit/pan/zoom, view presets, fit, translucent meshes, and an
  edge colour per mesh (used for selection, sketches and plane outlines).
  Vulkan first, GL fallback.
- **App shell** (`openshape` binary, egui/eframe; ADR 0012), in the familiar
  parametric-CAD layout with OpenShape's own icons, palette and name:
  - top bar: application menu (new, open, save, export STL, quit),
    document name, undo/redo, command search (Ctrl+K, fuzzy matching over
    every command with its shortcut), display modes, perspective/ortho,
    fit;
  - feature toolbar (sketch with a plane menu, extrude, revolve, box,
    templates, roll to end); in sketch mode the sketch toolbar (tools,
    construction, dimension, 12 constraints, Finish);
  - feature list: origin and the three default planes (show/hide, hover
    highlight, double-click to sketch), the features with status icons and
    tooltips, the draggable rollback bar, the parts with show/hide, and a
    properties card (volume, area, centroid, topology) for the selected
    part; inline rename; context menus (edit, rename, suppress, roll history
    here, move up/down across the bar, delete);
  - graphics area: the default planes as translucent outlined squares sized
    to the model, part and plane picking (click to select, double-click to
    edit), sketches drawn in 3D and hidden by the parts in front of them,
    the view cube, a context menu (views, display modes, show/hide
    sketches and planes), and messages;
  - floating feature dialogs with accept and cancel in their header; Enter
    accepts, Escape cancels; editing a feature shows the history up to it
    and greys the features after it (ADR 0013);
  - element tabs: add, switch, rename (double-click), delete;
  - right-drag orbits, middle-drag or Ctrl+right-drag pans, the wheel
    zooms, in both Part Studio and sketch mode.
  `openshape --screenshots <dir>` writes review screenshots of the whole
  app rendered offscreen.
- **Sketch mode**: line (chained), corner, centre-point and aligned
  rectangles, circle, 3-point, centre and tangent arcs, slot, regular
  polygon (3 to 64 sides), point, through-point spline, construction
  toggle (design note 0007 §4; each shape's degrees of freedom and the
  absence of redundant constraints are tested); 12
  constraint tools and dimensions (length, distance aligned/horizontal/
  vertical, angle, radius, diameter) typed with units; snapping to points,
  curves and the origin; while drawing a line, inference of horizontal,
  vertical, tangent (leaving an arc or spline) and perpendicular or
  parallel (to the line it continues or one the pointer passed); dragging with
  live solving; box selection and delete; colour-coded state (free, fully
  constrained, redundant, conflicting) with constraint glyphs and dimension
  labels placed outside the sketch. Trim (a piece between crossings; an
  uncrossed curve is deleted; a circle cut twice becomes an arc), extend (a
  line along its ray, an arc around its circle, to the next curve) and split
  (design note 0007 §7): cut ends are put on the crossing curve, the far end
  of a trimmed middle piece moves to a new curve with every constraint it
  had (tangencies included), curves that met at a moved end stay joined,
  and only constraints that no longer mean anything (a moved end's, a
  segment's length) are removed and counted; each edit is refused rather
  than committed if the sketch would stop solving or gain a redundant
  constraint. Fuzzed on random lines, arcs, circles, rectangles, slots,
  polygons and joined chains (12,000 sketches, about 15,000 edits, every
  one solved with nothing unplanned moving). Every edit is one undoable op; Cancel
  undoes the whole sketch session. The editor is UI-toolkit independent and
  tested end to end by driving it with screen clicks.
- **UI kit** (`openshape-ui-kit`): the application's look in one place
  (ADR 0012): theme tokens with pinned contrast ratios, 73 original icons
  drawn in code, a view cube, the command palette and search field, and an
  offscreen egui renderer used for screenshots and UI tests.
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
- The windowed app has not yet been run on real hardware. Building the
  offscreen UI renderer exposed a Phase 1 defect there: egui samples native
  textures as gamma-encoded bytes, but the app handed it the 3D view's sRGB
  view, so the window showed the 3D view darker than the offscreen renders.
  The view egui reads is now a plain `Rgba8Unorm` view of the same texture
  (or, on GL, a non-sRGB texture), and `crates/ui/tests/viewport_colour.rs`
  checks that what egui shows matches the renderer's output within one code
  on Vulkan and GL.
- Phase 2, 2026-10-11, on Mesa lavapipe: fmt, clippy with `-D warnings`, 378
  tests and `cargo deny` pass. Stress runs beyond the default test counts,
  all with zero failures:
  - solver: about 400,000 constructive cases at scales 1e-4 to 1e3 and
    540,000 garbage-input fuzz cases (no panic, no NaN, every success
    verified independently, no success with a collapsed curve);
  - regions: 8 full runs of the arrangement property tests with random
    seeds against independent oracles (polygon union, exact grid cells,
    raster flood fill), every successful extrude matching area × height.
- Integration found and fixed a solver defect before it shipped: a line
  constrained both horizontal and vertical collapsed to a point and the
  solver reported success. Any result with a zero-length line or a
  zero-radius circle or arc now fails, and the constraint that caused the
  collapse is reported as conflicting.
- App shell, 2026-10-12, on Mesa lavapipe: fmt, clippy with `-D warnings`,
  455 tests and `cargo deny` pass. `crates/ui/tests/shell.rs` drives the
  real application offscreen with synthetic mouse and keyboard input (13
  end-to-end tests: toolbar, dialogs, sketch from a plane row, cancel,
  edit-in-place, dragging the rollback bar, moving features across it,
  context menus, tabs, command search, picking, view cube, dimension
  placement). Writing them found two real defects, both fixed: Enter never
  reached feature dialogs (the global shortcuts consumed it), and the screen
  shots showed sketches drawn over the parts that hide them (they are now
  rendered in 3D with depth).

## Not yet / known limits

- Sketch editing tools are not done: offset, sketch fillet, mirror,
  patterns (in progress, design note 0007), and ellipses, conics, text,
  use/project. Trim cuts one piece per click (no drag across several
  curves yet), and splines can be trimmed only whole: their pieces are not
  exact splines of the sketch's two kinds (§7.5).
- Shell gaps: features reorder with Move up/down but not by dragging rows;
  dimension labels are placed automatically (outside the sketch) and cannot
  be dragged; feature dialogs do not preview the result live, and there is
  no toggle to see the final model while editing; parts highlight when
  selected but not on hover; there is no box selection of parts; shortcuts
  are not remappable; only the light theme exists; File > Open and Save use
  a text path prompt, not a native file picker.
- Silhouettes of curved faces are not drawn as lines (only real edges), so
  a cylinder's sides show by shading alone.
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
- No STEP or other exchange formats beyond STL (Phase 4; see docs/FORMATS.md).

## Next steps

1. Sketch editing tools (design note 0007): sketch fillet, mirror, offset,
   linear and circular patterns; draggable dimension labels.
2. Live preview in feature dialogs.
3. Re-run `scripts/gpu-check.sh` on RADV (docs/huginn-tasks.md prompt 1) and
   run the windowed app on a desktop.
4. Phase 3 groundwork: design notes for curve/surface intersection, booleans
   and persistent naming; make the kernel's tessellation and profile checks
   adapt to local feature size.
