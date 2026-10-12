# OpenShape status

Updated 2026-10-12, Phase 2 (sketcher): core sketcher, the application
shell, the editing tools, ellipses, conics, transform, curvature
continuity, text, smooth closed splines, trimming of ellipses and conics,
and curvature combs done; region detection now absorbs faces thinner than
the tolerance instead of refusing them. Phase 3 groundwork next.

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
    arcs, circles, splines, ellipses and conics as exact rational
    quadratics) and revolve (full and partial, profiles touching
    the axis, offset axes). Periodic faces split at ≤ π (ADR 0008).
  - Surface classification: plane, cylinder, cone, sphere, torus, extrusion,
    revolution, NURBS.
  - Watertight tessellation with chord and angle tolerances (shared edge
    discretisation, constrained Delaunay per face, singular points handled).
    It marks seam edges, where two pieces of one surface meet (the splits of
    ADR 0008) or where a sweep's faces from two profile curves meet at a
    point where the curves have equal first derivatives (a closed spline
    cut in two), so the display draws only creases; the comparison is
    geometric and errs towards showing an edge (a tangent joint between
    different curves stays visible). A NURBS split a rounding step from a
    knot splits at the knot, so no piece ends on a sliver knot span.
  - Exact mass properties (volume, area, centroid, inertia) by boundary
    integration with adaptive Gauss–Legendre quadrature (design note
    0003): accurate to the relative accuracy the caller asks (down to
    1e-13), measured against the boundary integral of the integrand's
    magnitude; a computation that cannot reach it fails with an error
    instead of returning an unvouched value. The golden corpus checks
    1e-9 relative and better.
- **Features** (`openshape-features`): regenerates a Part Studio through the
  kernel with per-feature status (ok, warning, suppressed, failed with
  message). A failing feature never blocks the others. Sketch features are
  solved and their regions computed during regeneration; extrude and revolve
  take their profile from picked sketch regions, which are named
  topologically and re-found after sketch edits (design note 0006 §5).
- **Sketch solver** (`openshape-sketch-solver`, on the deterministic
  `openshape-linalg`): points, lines, circles, arcs, through-point and
  control-point splines, open or smooth closed (ADR 0016: the periodic
  spline held as its exact clamped equivalent, a fixed linear map of the
  sketch's points), ellipses and elliptical arcs (design note 0008),
  conic arcs with their ρ (design note 0009) under 16 geometric constraint
  kinds (curvature continuity the latest, design note 0010) and 6
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
    and greys the features after it (ADR 0013); the dialog's values are
    previewed live (the feature regenerated as if applied, highlighted, a
    regeneration failure shown in the dialog), with nothing recorded until
    it is accepted, and STL export always writes the document, never an
    editor's view of it;
  - element tabs: add, switch, rename (double-click), delete;
  - right-drag orbits, middle-drag or Ctrl+right-drag pans, the wheel
    zooms, in both Part Studio and sketch mode.
  `openshape --screenshots <dir>` writes review screenshots of the whole
  app rendered offscreen.
- **Sketch mode**: line (chained), corner, centre-point and aligned
  rectangles, circle, 3-point, centre and tangent arcs, slot, regular
  polygon (3 to 64 sides), point, through-point and control-point splines
  (the control polygon shown dashed; clicking the first point again
  closes the curve smoothly, previewed while the pointer is there, ADR
  0016; a smooth closed spline has no ends, so tangency and curvature
  continuity at it are refused with the reason), ellipse and elliptical arc (design note
  0008: centre, axis end, then the other axis's extent, with
  construction semi-axes that can be dimensioned like any line; an arc's
  ends are clicked anywhere around it and it runs the way the pointer
  turned; both handednesses; a line can be made tangent to an ellipse and
  curves tangent at an elliptical arc's end; mirror and both patterns
  copy them; the kernel receives them exactly, so an ellipse's region is
  πab and its extrusion πab·h to rounding), conic (design note 0009: its
  two ends, the corner where its end tangents meet, then its shoulder,
  which the curve reaches up to; its control polygon shown dashed; its ρ
  dimensioned as a plain number, "ρ 0.42", or set by dragging the curve;
  lines and tangent arcs continue tangentially from its ends; mirror and
  both patterns copy it with its ρ, which a dimension on the source drives
  exactly; the kernel receives it exactly, so a segment's region is the
  closed form r(w) times its control triangle), construction
  toggle (design note 0007 §4; each shape's degrees of freedom and the
  absence of redundant constraints are tested); 12
  constraint tools and dimensions (length, distance aligned/horizontal/
  vertical, angle, radius, diameter) typed with units; snapping to points,
  curves and the origin; while drawing a line, inference of horizontal,
  vertical, tangent (leaving an arc or spline) and perpendicular or
  parallel (to the line it continues or one the pointer passed); dragging with
  live solving; box selection and delete; colour-coded state (free, fully
  constrained, redundant, conflicting) with constraint glyphs and dimension
  labels placed outside the sketch automatically, or dragged where the user
  wants them (stored per dimension as an offset from what it measures, so
  they move with the geometry; one undo step per drag; a command returns
  them to automatic placement; design note 0007 §12). Trim (a piece between crossings; an
  uncrossed curve is deleted; a circle or ellipse cut twice becomes an arc
  or elliptical arc on the same points), extend (a line along its ray, an
  arc around its circle, an elliptical arc around its ellipse, to the next
  curve, a conic along the same conic) and split (design note 0007 §7;
  ellipses 0008 §8.1, where an ellipse's own semi-axes do not cut it;
  conics 0009 §8.1, every piece an exact conic with its own corner and ρ,
  split pieces held tangent): cut ends are put on the crossing curve, the far end
  of a trimmed middle piece moves to a new curve with every constraint it
  had (tangencies included), curves that met at a moved end stay joined,
  and only constraints that no longer mean anything (a moved end's, a
  segment's length) are removed and counted; each edit is refused rather
  than committed if the sketch would stop solving or gain a redundant
  constraint. Fuzzed on random lines, arcs, circles, ellipses, elliptical
  arcs, conics, rectangles, slots, polygons and joined chains (12,000
  sketches, about 15,000 edits, every one solved with nothing unplanned
  moving; a 2,000-case run with ellipses committed 1,469 ellipse edits,
  and a 5,000-case run 2,477 conic edits, every piece on its conic). Sketch fillet rounds the
  corner of two lines or arcs with a tangent arc of the radius typed in the
  sketch dialog (line–line, line–arc, arc–arc; the centre from the curves'
  offsets); constraints on the corner, a filleted line's length and
  construction curves ending there (a centre-point rectangle's diagonals)
  move to a virtual sharp, so a dimensioned shape keeps its size and the
  degrees of freedom do not change; fuzzed over about 4,000 fillets of
  rectangles, polygons, chains and line/arc corners. Mirror (design note
  0007 §9.1) reflects a selection, or curves and points clicked one at a
  time, about a line; one derived solver constraint per line keeps every
  image the exact reflection of its source with no added degrees of
  freedom or redundancy, images share the sources' points on the line (so
  a mirrored half profile closes) and replace the point-on-line
  constraints there, and deleting a source or an image drops just that
  pair; fuzzed over 3,000 random batches of rectangles, polygons, slots,
  circles, arcs and points about random lines. Offset (design note 0007
  §10.6) moves a chain of lines and arcs joined end to end, or a circle,
  to the side the pointer is on; the distance is then typed and stays an
  editable dimension. One derived constraint holds every image (joined at
  corners where the offsets meet, through smooth joints on the normal,
  free ends square to the source) with no added freedom; changing the
  distance moves only the offset; after every solve the solver checks the
  offset is still valid (no piece turned back, no crossing, smooth joints
  still smooth) and fails naming it rather than returning a broken one.
  Fuzzed over 3,000 random rectangles, polygons, slots, circles and line
  chains: 2,738 offsets committed exact (no source moved, no crossing by
  the kernel's own intersections), the rest refused with the reason.
  Linear and circular patterns (design note 0007 §11.3) copy a selection
  2 to 100 times along a step between two clicks (a construction line
  whose length is then typed) or about a centre (copies spread at
  360°/count, the angle then typed and kept as a dimension); one derived
  constraint per use holds every copy exactly with no added freedom,
  spokes share their hub, a new step or angle moves only the copies, and
  deleting a copy drops just that copy.
  Copies that nothing else holds are computed from their sources after
  the solve rather than solved for, so 100 copies of a hexagon solve in
  7 ms (release). Fuzzed over 3,000 random sets of rectangles, polygons,
  slots, circles, arcs and points: 5,577 patterns committed (2,850 linear,
  2,727 circular), every one exact with no source moved and only the
  anchor's freedom added, and exact again after a new step or angle; 428
  were refused with the reason (nothing to pattern, or all of it on the
  centre).
  Move, rotate and scale (design note 0007 §13, ADR 0014) apply a
  similarity to a selection, by clicked points with a live preview or by
  typed offsets, angle or factor, optionally as a copy. The selection
  lands exactly (a translation adds the offset; a turn within 1e-9 rad of
  a quarter turn is snapped to it): its dimensions scale exactly, angles
  and ρ stay, fixes follow their point, quarter turns swap horizontal and
  vertical, and other turns keep lines that were horizontal or vertical
  parallel or perpendicular to one another, freeing only the rotation.
  Relations to the rest of the sketch are kept (the rest follows, as in a
  drag) or, with Keep constraints off, removed; dimensions to the rest are
  measured anew; copies get the constraints among the selection only.
  Every removed constraint is counted in the status message, and a
  transform the rest cannot follow is refused with the reason. Fuzzed
  over 3,000 random sketches (rectangles, circles, links, a mirror, a
  circular pattern, distances) and random moves, turns and scales with and
  without Copy and Keep: 2,579 committed, every moved point within
  6.4e-14 m of where the similarity puts it, no new redundancy, detached
  geometry unmoved, degrees of freedom changed only when the message says
  so; 422 refused with the reason (243 held by kept constraints, 176
  identities or copies onto the original, 3 that would collapse a line
  joined to the selection).
  Curvature continuity (G2, design note 0010) joins a spline to a line,
  arc, elliptical arc, conic or another spline at a shared end: one
  tangent line, the away tangents opposite (a cusp is never accepted), and
  equal curvature through the joint. The spline adapts and the other curve
  keeps its shape: the spline's points weigh 1e-3 in the solver's step
  metric, under a trust region of a quarter of the sketch size per step;
  only if that finds no solution does every curve give. Corners from 10°
  to 179° solve. Measured over 11,000 joints of every kind with every
  point moved by up to 15 % of the sketch size: 30 fail (0.27 %), all
  reported, none wrong; the other curve's curvature at the joint ends
  within 25 % of where it would without the constraint in 98.3 % of
  joints. Before the spline-first rule the same joints failed 1.0 % of
  the time and the solver flattened an ellipse a hundredfold in 10 % of
  elliptical joints. The tool replaces a tangency it implies, the joint
  shows "C", and trim and split hand the constraint on to the piece at
  the joint.
  Curvature combs show how curves bend: select arcs, circles, ellipses,
  conics or splines and choose Curvature comb (command search, or
  right-click a curve in the sketch, which selects it and offers "Show
  curvature comb" or "Hide curvature comb"). Each curve gets 96 teeth
  evenly spaced by arc length, standing on the outside of the turn (away
  from the centre of curvature) and as long as the curvature there, with
  their tips joined; one scale serves all combs shown, the tightest bend
  drawing 72 points long. A spline's comb shows both sides of a curvature
  jump at a knot. Across a curvature-continuous joint the teeth carry on;
  across a merely tangent one they jump. The samples come from the
  kernel (`Kernel::curvature_samples`), from exact first and second
  derivatives; combs are a view setting of the editing session, kept
  through edits but not saved with the sketch.
  Text (ADR 0015) places one line of characters in a box: type it, pick
  one of six bundled faces (Liberation Sans, Serif and Mono, regular and
  bold), click a corner, then click its height; the box grows towards the
  pointer and reads left to right. The box is four construction lines with
  a horizontal bottom, so the text has three freedoms (where it is, how
  tall) and is constrained and dimensioned through its corners and box
  lines like any geometry; a width dimension drives the height. The
  letters are the font's own quadratic and cubic outlines, shaped with
  kerning by harfrust (HarfBuzz's Rust shaper) and stored in the document
  (so a font or shaper update never changes an existing model); their
  regions extrude to their exact shapes, with holes, and keep their names
  when the text moves, turns, grows or gets more characters. Selecting a
  text opens its characters and font for editing (the box keeps its
  corner, height and direction); mirror images and pattern copies take the
  new text. Characters a font cannot draw are refused by name before
  anything is placed. Letters take no snaps, constraints, trim, extend,
  split or fillet (each refused with the reason); deleting a text deletes
  its box.
  Every edit is one undoable op; Cancel
  undoes the whole sketch session. The editor is UI-toolkit independent and
  tested end to end by driving it with screen clicks.
- **UI kit** (`openshape-ui-kit`): the application's look in one place
  (ADR 0012): theme tokens with pinned contrast ratios, 82 original icons
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
- Re-run on the same machine, 2026-10-12, commit 4de454a: build, demo
  (12/12), benchmark (0.153 ms mean, 0.224 ms p99 at 1080p, one run) and
  render timing (1,000 instances 0.3 ms per frame) pass; the test step
  stopped at its first failing binary after 307 passes, so the render
  tests (and the sRGB confirmation) did not run. The failure was in a
  test's reference formula, not the product: huginn builds with Fedora's
  Rust 1.95, whose `f64::acosh` computes ln(w + √(w² − 1)) and is 25 units
  in the last place off at w = 1 + 1e-5 (reproduced here; Rust 1.97's is
  correctly rounded), and the conic segment-area reference cancels that
  into a 4e-10 error against a 1e-10 check. The conic tests now compute
  arcosh with `ln_1p`; no product code uses `acosh`. The whole suite then
  passed on Rust 1.95 here (780 of 780, release), and
  `scripts/gpu-check.sh` now runs every test binary and lists each
  failing test.
- Second machine, 2026-10-12, commit 4de454a: AMD Radeon RX 7600 (RADV
  NAVI33, Mesa 26.2.4, Fedora 44, rustup Rust 1.99): every step passed,
  780 of 780 tests, so the sRGB test's bound holds on RADV (on this card;
  and on the 8060S, below). Demo 12/12; 1080p benchmark 0.215
  ms mean, 0.301 ms p99 (4,660 fps); 1,000 instances 0.9 ms per frame
  steady including readback.
- 8060S re-run, 2026-10-12, commit 9cbfe98, on Fedora's Rust 1.95:
  every step passed, 780 of 780 tests (so the sRGB bound holds on RADV
  STRIX_HALO too), demo 12/12; 1080p benchmark 0.154 ms mean, 0.184 ms
  p99 (6,470 fps); 1,000 instances 0.3 ms per frame steady including
  readback. Both real GPUs are now fully green.
- The windowed app ran on a desktop for the first time on 2026-10-12, on
  an AMD RX 7600 (RADV NAVI33) under GNOME on Wayland, with the 3D view's
  brightness right. Its title bar was winit's bare fallback frame (GNOME
  leaves decorations to applications); it is now winit's GNOME-style bar
  (ADR 0012 addendum). Building the
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
- The mirror fuzz found and fixed a second solver defect, in diagnosis
  only: a constraint implied exactly by others that the solver eliminates
  (an angle of 90° between a horizontal and a vertical line) could be the
  only equation of its component, and its Jacobian row, rounding noise of
  about 1e-32, counted as full rank against itself. The constraint was
  reported redundant when the input was nudged and not when it was placed
  exactly. Rows that cancel to within 1e-10 of their own terms are now
  zero (design note 0005 §9), with a regression test on both inputs.
- The offset tests found a behaviour defect in the solver: changing an
  offset's distance also moved the curves it follows (by a third of the
  change in a rectangle test), because every step was the smallest change
  of all points together. Images of derived constraints now move first:
  when moving them alone solves a step, nothing else moves (design note
  0005 §3), and the test asserts the sources stay exactly where they were.
- The pattern tests found two more. Typing a linear pattern's step resized
  the square being copied: moving the images alone cannot absorb a change
  of the direction line, so the old rule fell back to the ordinary step.
  Images now weigh nothing in every step, computed exactly as the limit
  of the weighted step (design note 0005 §5), with a regression test that
  fails on the old rule. And the first editor fuzz took 28 minutes: a
  pattern was one dense component of every copy's unknowns (100 hexagons,
  2 s a solve). Copies held only by their pattern are now peeled from the
  system and computed afterwards, exactly, with their free state from the
  same combination of null-space rows: 7 ms, and the fuzz takes seconds.
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

- Text, 2026-10-12, on Mesa lavapipe: fmt, clippy with `-D warnings`,
  756 tests and `cargo deny` pass. The font tests pin each bundled
  font's hash and licence, lay out closed contours in every face, compare
  the stored pieces of nine awkward glyphs (O, S, g, 8, @, e, Q, B, &)
  with the font's own outline read independently, and check holes,
  kerning, monospace widths and determinism; the solver
  tests check the box's freedoms, redundancy and conflicts, and letters
  mapped exactly into turned and mirrored boxes; the features tests
  extrude words letter by letter (each volume within 1e-6 of its
  region's exact area times the depth) and follow them through resizes, turns and longer text; nine
  end-to-end tests drive the tool, the edit dialog, mirror, patterns,
  deletion and every refusal. Extruding a word exposed the cost of exact
  mass properties on curved faces: edges were split to the requested
  accuracy before quadrature, so a 15-letter word at 1e-9 took 128 s and
  one extruded ellipse 20 s in release. Integration is now adaptive
  (Gander & Gautschi), the three text feature tests take 1.8 s in a
  debug build (557 s before), and
  the kernel tests that ask 1e-9 or finer were checked against the
  closed forms they already used. Building it found two defects before
  they shipped: a draft that stopped refining when bisection stalled
  returned a strongly rational conic's volume 21 % small, and a spline
  split one rounding step from a knot left a 5.6e-17 knot span on which
  quadrature never converged; both are regression tests.
- Smooth closed splines, 2026-10-12, on Mesa lavapipe: fmt, clippy with
  `-D warnings`, 780 tests and `cargo deny` pass. The solver tests check
  the closed control spline against the uniform periodic B-spline on the
  wrapped polygon and the closed through-point spline against the
  periodic cubic interpolant in second-derivative form (de Boor), both
  evaluated in the tests independently of the solver, plus freedoms, a
  point sliding across the seam, refusals, mirror images and old
  documents reading unchanged; eight end-to-end tests close both tools
  on their first point and compare the shape, the region area and the
  extruded volume with the same independent constructions, check the
  preview is the curve the click makes, the refusals, mirror, pattern,
  copy and undo. The screenshots showed two lines down the smooth side
  of an extruded closed spline: where the profile loop cuts the curve in
  two the faces meet smoothly, but the display could only hide seams
  inside one surface. The kernel now hides an edge where swept pieces
  meet with equal first derivatives (golden tests: smooth, tangent-only
  and cornered joints, extruded and revolved). The through-point side
  still showed one line: the cut at the parameter midpoint landed a
  rounding step from a knot (0.5 against 0.49999999999999994), leaving a
  sliver span whose end derivative was meaningless. Such splits now snap
  to the knot. The closure's end derivatives, which matched to 7e-14
  from the solve, are now equal by construction, so the display test
  does not depend on how well the interpolation system is conditioned.
- Trimming conics, 2026-10-12, on Mesa lavapipe: fmt, clippy with
  `-D warnings`, 795 tests and `cargo deny` pass. The solver tests
  check pieces against the conic's implicit equation inside and past the
  arc, with their end tangents, and parameter round trips; five end-to-end
  tests drive trim, split and extend on elliptic, parabolic and
  hyperbolic conics (pieces on the old conic, segment areas against the
  closed form, constraints handed on, freedoms, refusals); a wrong weight
  normalisation fails four of them. Adding conics to the trim fuzz, with
  a line held tangent into a conic's start, found two defects, both
  fixed with regression tests. The kernel reported that one tangent touch
  as three contacts up to 3.7 µm apart (its refinement stalls at
  tangencies, where the curves stay within the tolerance for about
  0.4 mm on a gently curved conic), so extend moved the line's end 2 µm
  and dropped its tangency; tangencies are now finished by Newton's
  method on a system whose root is simple there (design note 0007 §6).
  And a 2 µm piece of a conic, a straight segment to the tolerance, was
  flattened by the solver; such pieces are now refused.
  CI had failed once on main (9cbfe98) in the random region test
  (lines, circles and arcs on a lattice, perturbed by up to three
  tolerances), unrelated to this work. Running it 20,000 times instead of
  48 found four rare failures that are documented refusals, not wrong
  results, which the test now allows for perturbed input only (counted
  and printed, seeds kept as regressions): a triangle of 1e-13 m² or a
  few microns of boundary that the profile builder refuses below its
  resolution, a 0.3 µm sliver whose mass properties cannot show 1e-10
  (computed to 1e-6), and a sliver whose two sides the validator's coarse
  tessellation turns into one chord each, reading as inside-out (the
  documented extrude limitation). A fifth is a real kernel defect, about
  once in 2,500 perturbed cases: a line tangent to a circle at its end,
  with a second line starting 0.11 µm from that end, made the
  arrangement inconsistent (Euler's formula fails) and it refused the
  sketch (see the known issues).
- Curvature combs and region robustness, 2026-10-12, on Mesa lavapipe:
  fmt, clippy with `-D warnings`, 819 tests and `cargo deny` pass. The
  comb's samples are checked against closed forms (circle, ellipse,
  parabola, equal arc lengths by an independent quadrature) and against
  exact one-sided differences at spline knots; the editor tests check a
  G2 joint's two sides against the solved curves' end derivatives, and
  the right-click menu showing and hiding a comb in the rendered image.
  The region arrangement was stressed with the lattice test at 40,000
  cases instead of 48: it failed after about 2,500 cases at the start of
  the day (an internal inconsistency), and each failure it found since
  is fixed and kept as a regression test with the exact curves'
  regions: half-edges are ordered beyond a vertex's merged points, faces
  thinner than the tolerance that are artefacts of it are absorbed, and
  the profile builder refuses micron-sized loops by name. Successive
  runs passed 7,533, 13,908, 2,868, 8,360 and 14,673 cases before the
  next defect.
- Curvature continuity, 2026-10-12, on Mesa lavapipe: fmt, clippy with
  `-D warnings`, 717 tests and `cargo deny` pass. The solver tests
  check every joint against curvatures computed in the test from closed
  forms (Piegl & Tiller's end derivatives for splines, the quotient rule
  for conics, |A×B|/‖E′‖³ for ellipses), independent of the solver's
  equations; a deterministic statistics test bounds the failure rate and
  how often the other curve changes shape. Measuring found two defects
  before they shipped, both fixed: with every point weighing the same the
  solver met a wrongly bent spline by flattening the curve it joins (an
  arc's radius tripled; ellipses squashed), and an elliptical arc drawn
  with its minor axis on the other side, the same curve, failed 255 times
  in 2,000 against 13 (now 9 against 9). Two rare random-input test
  failures on main (an ellipse gradient checked with too coarse a
  difference step, and a rectangle drawn thinner than the trim's cut gap)
  were fixed separately.
- Transform, 2026-10-12, on Mesa lavapipe: fmt, clippy with
  `-D warnings`, 703 tests and `cargo deny` pass. Besides the fuzz above,
  14 end-to-end tests drive the three tools by clicks and typed values
  (snapping, previews, quarter turns, other turns, scaling areas by s²,
  copies, labels, every curve kind).
- Conics, 2026-10-12, on Mesa lavapipe: fmt, clippy with `-D warnings`,
  686 tests and `cargo deny` pass. The solver's conic points are checked
  against Lee's implicit equation of the control triangle, which is
  independent of the solver's parametric form; the fuzz generator gained
  a conic motif (20,000 cases per property, no failure), and 100,000
  perturbed re-solves per size fail no more often with conics than
  without (0, 2 and 19 at 1%, 2% and 5%, against 0, 1 and 19), always
  reported, never wrong. The kernel tests check segment areas against
  the closed form r(w) to 1e-12 (the closed form itself against an
  independent quadrature to 1e-10), a quarter ellipse drawn as a conic
  (πab/4), a paraboloid and a one-sheet hyperboloid revolved from conic
  arcs (volume and area to 1e-9), and 2,000 random extruded segments in
  release (valid, watertight, exact). Writing them found two defects,
  both fixed: the closest-point search stopped a rounding error short of
  an arc's end when that end was the answer, and mirror images could
  differ from their source's ρ in the last bit (now written back
  exactly). The test's own closed form needed `(1 − w)(1 + w)` for
  `1 − w²`: near w = 1 the rounding of w² cost five digits.
- Ellipses, 2026-10-12, on Mesa lavapipe: fmt, clippy with `-D warnings`,
  646 tests and `cargo deny` pass. The kernel tests check every curve
  point against the ellipse's implicit equation, exact areas of full
  ellipses, sectors and segments (1e-12 relative), elliptic cylinders,
  prolate and oblate spheroids and elliptic tori against closed forms;
  the editor tests were checked against deliberate defects (swapping the
  handedness rule, dropping the exact snap) and fail on each. Writing
  them found a kernel defect that predates ellipses: an extruded closed
  spline (an ellipse, or a closed control-point spline) could fail to
  mesh, because boundary samples at the curve's closure point landed on
  the wrong end of the periodic parameter range. They are now put on the
  face's side; three of the new tests fail without the fix.

## Not yet / known limits

- Sketch tools not done: use/project (needs Phase 3 faces).
  Text is one line, in six bundled faces only (no system fonts), with no
  right-to-left reordering of mixed scripts and no text on a curve; its
  letters cannot be trimmed or constrained individually. A word extruded
  in one feature makes one part with a solid per letter (each feature
  makes its own part until booleans).
  Curvature continuity fails on 0.27 % of heavily perturbed joints
  (mostly a spline meeting an interpolated spline), reported, never
  wrong; where the spline alone cannot reach the joint's curvature the
  other curve gives too (a conic may flatten) unless it is dimensioned.
  Transform has no on-screen manipulator (drag handles): the selection
  moves by clicked points or typed values.
  Conics cannot be made tangent to a curve away from their ends, so a
  sketch fillet at a conic's end is refused; a conic cannot be offset (its
  offset is not a conic). The two conics left by trimming the middle of
  one are independent afterwards (nothing says they are one conic), split
  pieces are only held tangent, and a piece whose corner would be within
  0.1 µm of its chord is refused as too short to keep its shape.
  The ρ label reads "ρ 0.42"; in the app's font at 13 points ρ looks much
  like a p.
  Ellipses cannot be made tangent to a circle, an arc or another ellipse
  except at a shared end point (no closed form; design note 0008 §4.5),
  so a sketch fillet at an elliptical arc's end is refused too. An ellipse
  cannot be offset (its offset is not a conic). Radius, diameter, equal
  and concentric do not apply to ellipses; their semi-axes are lines to
  dimension instead.
  The ellipse solver fails more often than the rest of the solver from
  far-off starts (§9 of the note: about 8 in 100,000 random sketches
  perturbed by up to 5% of their size, against 0 without ellipses),
  always reporting the failure, never a wrong result.
  Splines cannot be offset (their offsets are not splines; an
  approximation with a guaranteed tolerance is future work).
  A fillet removes (and reports) equal-length and midpoint constraints on
  the filleted lines, so a filleted regular polygon is no longer held
  regular; equal distances between virtual sharps are not expressible yet. Trim cuts one piece per click (no drag across several
  curves yet), and splines can be trimmed only whole: their pieces are not
  exact splines of the sketch's two kinds (§7.5).
- Shell gaps: features reorder with Move up/down but not by dragging rows;
  there is no toggle to see the final model while editing, and a
  new feature's preview is not listed in the feature list; parts highlight when
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
  size is on the kernel list. Mass properties have the same blind spot:
  their rounding floor scales with the body's size, not with the
  coordinates the curves are evaluated in, so a sliver 0.3 µm wide on an
  arc of 3 m radius cannot be integrated to 1e-10 and the kernel says
  so (it converges to 1e-6; the random region test computes such slivers
  there).
- Region detection refuses a sketch with a long sliver thinner than the
  tolerance that has three or more sides (design note 0006 §4.2 absorbs
  the other kinds), naming where; and the profile builder, which merges
  end points within ten tolerances, refuses a region that small as a whole
  (a loop within twenty tolerances of a point). Both are errors, never
  wrong regions; in the random lattice test they occur only for curves
  moved by about a tolerance or more.
- Solver conflict detection is local: from a start far from any solution it
  can call a sketch conflicting although another branch would satisfy it.
- No booleans: each feature makes its own part (Phase 3).
- No persistent naming or incremental regeneration: every edit regenerates
  the whole studio (Phase 3). Features don't reference model geometry yet, so
  nothing can break on upstream edits today.
- Non-manifold topology deferred to Phase 3.
- No STEP or other exchange formats beyond STL (Phase 4; see docs/FORMATS.md).

## Next steps

1. Check the new title bar on GNOME (on Pat's desktop).
2. Phase 3 groundwork: design notes for curve/surface intersection, booleans
   and persistent naming; make the kernel's tessellation and profile checks
   adapt to local feature size, and scale the mass-property rounding floor
   to the coordinates curves are evaluated in. The rest of the Phase 2
   sketcher (project/convert, intersection curves, pierce) needs Phase 3
   faces.
