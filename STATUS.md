# OpenDatum status

Updated 2026-10-16. Phase 2 (sketcher) is done but for the tools that
need model faces (project/convert, intersection curves, pierce). Phase 3
has begun: the design notes for intersections (0011), booleans (0012)
and persistent naming (0013) are written; every region of the random
sketch tests extrudes, slivers included, with mass properties exact to
rounding; and persistent naming is built through its first references:
the kernel's origin tables (ADR 0018), the `opendatum-naming` crate,
which names the faces, edges and vertices of boxes, extrudes and
revolves as the user sees them and resolves stored references, and a
sketch on a picked planar face and an extrude up to a picked parallel
face, from the viewport through the model (ADR 0019). Booleans combine
solids bounded by planes, cylinders, spheres and cones where their faces
meet in lines, circles and conics (design note 0012 §9, §11, §13 and
§15, ADR 0020): union, intersection and difference, with planes across
cones at any angle (ellipses, parabolas, hyperbolas, and the rulings of
a plane through a frustum's apex), checked by the volume identities and
point membership on tens of thousands of random pairs. Bodies with
those conic edges on cylinders and cones are operands of further
booleans (a second countersink beside a tilted one, a pin flattened
twice, chains of cuts; §16, §17, ADR 0023). Bodies whose edges and
corners rest on flat faces without crossing them are decided (§18, §19,
ADR 0024; the faces there all flat): the contact leaves no trace, or
the boolean is refused as non-manifold where the result would pinch; a fuzz of such contacts
swept across the tolerance found 24,571 faults on main in 264,000
operations (results touching themselves, wrong origins, kernel errors)
and none in this build. A fin under a tolerance thick, which main
returns where an edge lies beside the line where its face crosses a
face of the other body, is refused (§20). Next: the rest of the
contacts (curved faces, surfaces touching), then the features that use
booleans.
The product is now named OpenDatum (ADR 0017): packages `opendatum-*`,
binary `opendatum`, documents `.odatum`; documents saved under the old
name open unchanged.

## What works

- **Op log** (`opendatum-oplog`): Lamport/replica op IDs with causal deps,
  deterministic replay, undo/redo stored as ops (redo = undo of undo),
  history at any set of heads, sync primitives (`missing_for`, `integrate`
  with a pending buffer for out-of-order ops), CBOR container. Property tests
  check undo/redo against a reference model and replica convergence.
- **Document model** (`opendatum-core-model`): documents with element tabs, a
  Part Studio feature list as data (insert, edit, rename, suppress, move,
  delete), a rollback bar stored in the op log (ADR 0013; its rules are
  property-tested against an independent model), element tabs that can be
  added, renamed and deleted, a reducer over the op log, units parsing of single quantities
  (`"25.4 mm"`, `"1 in"`, `"90 deg"`; bare numbers use the default unit).
  Expressions and variables come with Variable Studios.
- **NURBS** (`opendatum-nurbs`, kernel-private): rational curves and surfaces
  after Piegl & Tiller — evaluation, derivatives, knot insertion, split,
  closest point, exact arcs, interpolation, extrusion/revolution surfaces.
- **Kernel** (`opendatum-kernel`, behind the `Kernel` trait):
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
    Every triangle of a curved face is refined until it faces the way
    the surface does at its corners and centroid (unless its corners lie
    in a line within the model tolerance) and is within the chord
    tolerance of the surface as estimated, not proven: measured at
    thirteen points of it and, between them, by the quartic through them,
    bounded by its Bernstein-Bézier coefficients (design note 0002 §4).
    Unit tests measure the meshes against exact distances on a grid of
    24ths of every triangle (a ball, polar caps of it, a torus, solids
    turned from S-shaped splines: 0.96 to 0.999 of the tolerance) and
    against the
    surfaces' closest points at fifteen points (cones to the apex, balls
    to the poles, faces cut across their rulings, thin wedges, a revolved
    cubic bent nearly to a crease). A triangle too small to split that
    folds over is mended by an edge flip, or else split at its centroid; a
    sliver too thin for either, refinement that cannot finish, or a face
    past two million vertices is an error. On four benchmark bodies (a
    ball, a ring, a turned S-spline, a creased revolve), meshing takes 1.6
    to 2.4 times as long a triangle as before this was built (docs/
    PERFORMANCE.md, fifth pass).
    It marks seam edges, where two pieces of one surface meet (the splits of
    ADR 0008) or where a sweep's faces from two profile curves meet at a
    point where the curves have equal first derivatives (a closed spline
    cut in two), so the display draws only creases; the comparison is
    geometric and errs towards showing an edge (a tangent joint between
    different curves stays visible). A NURBS split a rounding step from a
    knot splits at the knot, so no piece ends on a sliver knot span.
    The edge samples adapt to each face's boundary (design note 0002 §5):
    pieces of boundary that cross, or may come closer than their chords
    show, are bisected on the exact curves until the polygon has the
    face's topology, and where two curves leave a vertex tangentially (a
    cusp) their samples are matched at equal distances along the common
    direction, so their polygons cannot cross however close the curves
    run.
  - Profiles are checked for curves that touch or cross with the exact
    contact finder of the sketch arrangement (design note 0001), not on
    sampled polylines, and measured with its exact areas and winding
    numbers; the refusal names the place.
  - Exact mass properties (volume, area, centroid, inertia) by boundary
    integration with adaptive Gauss–Legendre quadrature (design note
    0003): accurate to the relative accuracy the caller asks (down to
    1e-13), measured against the boundary integral of the integrand's
    magnitude; a computation that cannot reach it fails with an error
    instead of returning an unvouched value. The golden corpus checks
    1e-9 relative and better. Edges that lie on their two faces only
    within the tolerance (merged end points) get the flux of the strip
    between the faces added, so such bodies are exact to rounding too.
  - Origin tables (ADR 0018): every face, edge and vertex of a box,
    extrude or revolve records how it arose (the profile piece and the
    sweep positions it spans, or the region a cap covers), addressed by
    per-body handles in mesh order; `Kernel::adjacency` gives each edge's
    faces, each vertex's edges and the seams, and the mesh each vertex's
    position. Tests check the tables complete, unique, consistent with
    the topology and right by an independent geometric check, on tilted
    planes, symmetric extrudes, and partial and full revolves either way
    round, and that they follow the curves, not their order.
  - Booleans (`Kernel::boolean`, design note 0012 §9, §11, §13, §15,
    §17 and §19, ADRs 0020 and 0024): union, intersection and difference of solids bounded by
    planes and circular cylinders (boxes, extrudes of profiles with
    lines, arcs and circles, holes included, at any orientation), with
    coincident faces, edges and vertices, voids and several solids; and
    with spheres and cones where every meeting is a circle (a ball with
    planes at any angle or with another ball anywhere; cones, cylinders
    and balls about one axis, with planes square to it: dimples,
    countersinks, pointed and rounded ends); and planes across cones at
    any angle, in whole and half ellipses, parabolas and hyperbolas
    (tilted countersinks, drafted pins flattened along their axes,
    sections through a countersink's axis), and in two rulings where a
    plane passes through the apex of a cone whose faces stop short of
    it. Bodies such cuts left are operands again: their conic edges,
    stored as exact rational quadratic NURBS curves, are followed as
    stored in the charts of cylinders and cones, with a band bounded
    from the curve's own derivatives and their rounding (ADR 0023,
    design note 0012 §16). Where a plane
    or a cylinder meets another surface the curve is found in closed form
    (lines, circles, ellipses, conics in vertex form with rational
    quadratic edges), once per pair of surfaces, and faces are
    split in charts in which lengths are lengths in space. Numerical
    decisions are made once (a sign per vertex and plane, clusters of
    points within the tolerance, one incidence decision per point and
    edge) and each piece of a face is classified by exact ray parity
    against the other body (Shewchuk's orientation predicate) at a point
    well inside it. Spatial indexes (a bounding-volume hierarchy) find the
    faces that meet and the edges and points near each other, and a work
    limit refuses a boolean too large to finish in reasonable time
    (`TooLarge`). Coplanar neighbours and
    collinear edges are merged; origins record which operand faces,
    edges and vertices each result entity comes from, and where faces
    meet. What the tolerance cannot decide is refused with a message and
    a place (non-manifold results, points or edges closer than the
    tolerance, slivers, faces running within the tolerance of the other
    surface, a plane through a pointed cone's apex where both bodies
    reach it). Where an edge or a corner of one body rests on a flat face of
    the other without crossing it, from outside or inside, the faces there
    all flat, the contact leaves no trace in the result, or the boolean is
    refused as non-manifold at the contact where the result would pinch
    there (a cube on its corner united with a plate, a wedge's edge in a
    groove); what the tolerance cannot place is refused with its own
    message, and every result is checked for a place within the tolerance of
    a flat face bounded by straight edges that it is not on (§19). Spheres
    and cones with cylinders or cones off their axes, cylinders at an angle
    to each other, and contacts that do not cross where a face at the
    contact is curved (an edge or a corner resting on a curved face, a
    curved edge on a flat one, faces touching faces) fail as not supported
    yet.
- **Features** (`opendatum-features`): regenerates a Part Studio through the
  kernel with per-feature status (ok, warning, suppressed, failed with
  message). A failing feature never blocks the others. Sketch features are
  solved and their regions computed during regeneration; extrude and revolve
  take their profile from picked sketch regions, which are named
  topologically and re-found after sketch edits (design note 0006 §5).
  Every part carries the names of its faces, edges and vertices, and
  features refer to earlier parts through them (ADR 0019): a sketch on a
  planar face gets a frame from the plane alone (seen from outside, the
  standard planes' frames on faces parallel to them, upright on faces
  that face sideways) and follows the face through upstream edits; an
  extrude ends on a planar face parallel to its sketch. A reference to a
  part that is gone, suppressed or later, to a face that no longer
  exists, is curved, or is not parallel, fails the feature with a
  message saying which; a body that cannot be named gets a warning.
- **Persistent names** (`opendatum-naming`, design note 0013 §11): names
  for the faces, edges and vertices of box, extrude and revolve bodies as
  the user sees them, made from the origin tables and the sketch keys of
  the profile: a side face by its sketch edge, a cap by the selected
  regions it covers (by the keys the feature stores), an edge by its two
  faces, a vertex by the faces around it, with a sketch-level tie-break
  only where those collide (a lens's two lateral edges, the two points of
  a half disc on its revolve axis). Kernel faces joined across seams are
  one face (π-split pieces, a circle's halves, a sketch line split by a
  vertex), named by all their sources. Resolution (§6): the exact name;
  else the entities descending from it (a sketch edge split, extended or
  re-keyed; regions merged into one cap or parted), all of them for a
  carrier or set reference and the one the hint lies on for a single
  one; else "no longer exists", never another entity of the same feature
  and role. Tests check every name against an independent geometric
  oracle (which faces an edge or vertex touches, which curve a side lies
  on, which end a cap is at and which region it covers), names unchanged
  by dimensions, curve order and revolve angles across π, and each
  resolution rule. Sketches on faces and extrudes up to faces resolve
  through them (see Features).
- **Sketch solver** (`opendatum-sketch-solver`, on the deterministic
  `opendatum-linalg`): points, lines, circles, arcs, through-point and
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
- **Sketch regions** (`opendatum-kernel`, `arrangement`): the planar
  arrangement of arbitrary sketch curves (lines, arcs, circles, rational
  splines; crossings, tangencies, overlaps, duplicates, self-intersections)
  with exact areas, holes, and region profiles ready for extrude and
  revolve.
- **Rendering** (`opendatum-render`): wgpu renderer for the window and an
  offscreen renderer; shaded parts with edges, perspective/orthographic
  camera, orbit/pan/zoom, view presets, fit, translucent meshes, and an
  edge colour per mesh (used for selection, sketches and plane outlines).
  Vulkan first, GL fallback.
- **App shell** (`opendatum` binary, egui/eframe; ADR 0012), in the familiar
  parametric-CAD layout with OpenDatum's own icons, palette and name:
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
    to the model; picking of the parts' faces, edges and vertices (ADR
    0019: a vertex within 7 points of the pointer wins, then an edge
    within 5, then the face under it; hidden ones are skipped) with the
    whole entity highlighted while hovered and when selected, and of
    planes (click to select, double-click to edit); sketches drawn in 3D
    and hidden by the parts in front of them, the view cube, a context
    menu (views, display modes, show/hide sketches and planes), and
    messages;
  - Sketch with a planar face selected starts a sketch on that face (a
    curved face is refused with a message); the extrude dialog ends at
    a blind depth or "Up to face", picked in the view and highlighted
    while the dialog is open;
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
  `opendatum --screenshots <dir>` writes review screenshots of the whole
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
- **UI kit** (`opendatum-ui-kit`): the application's look in one place
  (ADR 0012): theme tokens with pinned contrast ratios, 82 original icons
  drawn in code, a view cube, the command palette and search field, and an
  offscreen egui renderer used for screenshots and UI tests.
- **I/O** (`opendatum-io`): STL export, ASCII and binary.
- **Headless demo**: `cargo run -p opendatum -- --headless-demo out/` builds
  a box through the op log, orbits, adds three more parts, undoes and redoes,
  saves and reopens, exports STL, and checks each step pixel-exactly. All 12
  checks pass on Mesa lavapipe.
- **Benchmark**: `opendatum --bench <dir> [frames]` times offscreen frames of
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
- Rename and mass rounding floor, 2026-10-12, on Mesa lavapipe: fmt,
  clippy with `-D warnings`, 833 tests and `cargo deny` pass (832 for
  the rename to OpenDatum, ADR 0017, with a new op-log test that
  documents saved under the old name still open). The lattice test's
  40,000-case rerun with the vertex-order fix passed; it refused, as
  documented, 265 perturbed inputs unresolvable at the tolerance, 19
  profile loops below the builder's resolution and 46 pinched profiles,
  and 97 slivers fell back to 1e-6 for their mass. Mass properties now
  scale their rounding floor with the coordinates a body's points are
  computed from (design note 0003): the lens 0.3 µm wide between arcs of
  3 m radius that was refused at 1e-10 converges, within 1e-9 of its
  exact volume (regression test, fails without the change; it allows
  1e-8, what rounding of the points can cost), and the
  random region tests no longer let any sliver fall back to 1e-6. With
  the change, all 12 random arrangement tests passed at raised counts
  (2,000 polygons, 2,000 grids, 1,000 scales, 500 splines, 20,000
  garbage inputs and the rest), and a fresh 40,000-case lattice run
  passed with every region's mass to 1e-10, no fallback; it refused 313
  unresolvable perturbed inputs, 18 loops below resolution and 44
  pinched profiles (the cases are random per run).
- Booleans of planes and cylinders (design note 0012 §11), 2026-10-13,
  on Mesa lavapipe: fmt, clippy with `-D warnings`, 925 tests and
  `cargo deny` pass, and the release time guards pass as CI runs them.
  Golden cases exact to 1e-9 in volume and area, with topology and
  origins: a plate with a round hole, a boss through and on a plate, a
  D-cut shaft, two parallel cylinders (all three operations), a
  counterbore, voids with flat and curved faces, and a slab cut across
  a cylinder (two solids). 10,000 random pairs of boxes at any angle and
  prisms with arcs (40,000 operations): every result valid with its
  origins on the operands' surfaces, every volume identity holding where
  the results exist, 13 refusals in 1,000, nearly all tangent contacts.
  The run found one fault, fixed and kept as a golden case: a void
  inside a curved body was placed by a corner of the outer shell, so a
  body cut from inside another was refused. The polyhedral fuzz at 100
  times its default size refuses exactly the same operations as before.
  Drilling a round hole into a plate with 99 fell from 211 ms to 40 ms
  (docs/PERFORMANCE.md, third pass); CI now runs the curved fuzz too.
- A missed contact in the planes-and-cylinders booleans, found and fixed
  the same day (2026-10-13): a face on a cylinder whose boundary stayed
  outside another cylinder was taken to be apart from it, so a rod
  dipping into a pipe's wall between the pipe's seams (no edge of either
  crossing the other's surface) came back as two bodies that do not
  meet: a union of two overlapping solids, a difference that removed
  nothing, an empty intersection. The validator passed it and the volume
  identities hold for it. The rule now keeps a face apart only where the
  convex hull of its boundary does, or where the surfaces meet in closed
  form (lines along the axis, circles about it); the rod case is
  refused as a quartic (`a_rod_dipping_into_a_pipe_wall_is_refused`).
  The curved fuzz now also checks point membership on meshes that share
  no code with the booleans (24 points a pair where the operands' boxes
  meet, each inside a result exactly when the operation's rule says so),
  which catches the old rule, and has a second run of prisms along two
  different axes: 2,000 such pairs, every result correct, 383 refusals
  in 1,000 operations, all quartics. The parallel fuzz at 10 times its
  default size refuses exactly the same operations as before the fix.
- Wrong volumes for turned solids in the planes-and-cylinders booleans,
  found and fixed the same day (2026-10-13), while building cones and
  spheres. A revolve's flat cap is a plane recorded as a surface of
  revolution, in two faces of half a turn each (ADR 0008). Where a
  boolean kept both halves, it joined them into one face around the
  axis on that record, whose parameters cannot hold it: a turned rod cut
  square by a box came back valid with the wrong volume (π/2 kept as
  1.309; 3π/2 as 3.927). A joined face on such a plane is now recorded
  as a plane, and the validator refuses any face of revolution whose
  boundary leaves its frame's half turn, which alone refuses every case
  the fault made. Golden cases with exact volumes: turned rods and tubes
  cut by boxes across the top and down the side. The rod golden test
  the entry above names was missing from its commit (a scratch test went
  in instead); it is added now.
  A fuzz run with turned rods, tubes and stepped rods among the random
  solids (4,000 pairs along one axis and 2,000 at an angle, at 20 times
  the default size) found more, which also comes before the cap fix and
  is not yet fixed, so turned solids are not in the CI fuzz yet. Seven
  pairs along one axis break the volume identities by 3e-9 to 2e-6 (at
  most 3e-7 of the volumes, against 1e-9 allowed), with no point found
  on the wrong side: whether the results or their mass properties are
  off is not yet known. Loudly: the mass properties of 102 results do
  not converge; 8 operations are refused as "a curved segment could not
  be resolved into chords" and 4 as a piece the triangulation finds too
  small; and the tangent plane of a turned rod's cap at its centre (on
  the axis) is missing, which the fuzz's own check of origins needs (218
  operations). All fixed in the entry below.
- What turned solids met in the booleans, fixed the same day
  (2026-10-13; design note 0012 §11, "Turned solids"). Three faults, each
  kept as a test that fails without its fix: (1) a turned solid of half
  a turn has its two flat sides on one plane, joined into one face, and
  each cap's two radii were then joined into a diameter through the
  axis, where the cap's surface of revolution is singular; the boundary
  integrand of the mass properties has a kink there, so they did not
  converge, or (7 pairs) converged to volumes 3e-7 off. A point on the
  axis of a face's surface of revolution now stays a vertex, and the
  validator refuses an edge through such an axis between its ends. (2)
  A curve's end moved to its cluster by more than a sixteenth of the
  tolerance (clusters take in points within it) left its end chord
  halved without end, as the move is part of the chord's band that
  halving cannot narrow; chords are now halved by their sag alone (also
  in the search for a piece's point and the face outlines' point test,
  where the same move could have looped forever). (3) A cap's centre came
  out of its chart at 8e-49, which the triangulation refuses as too
  small; coordinates under 1e-30 of a polygon's size are taken as zero.
  Also: the tangent plane of a planar face is now exact everywhere,
  centres of caps included (it was missing there). Turned solids (a
  fifth of the random solids) are now in the curved fuzz and in CI. At
  20 times the default size: 4,000 pairs along one axis and 2,000 at an
  angle, every result valid and correct by the identities and by point
  membership; 20 refusals in 1,000 along one axis (312 tangent
  contacts, 5 non-manifold of 16,000) and, at an angle, all quartics
  but one non-manifold result and one elliptical edge. The same day a
  property test of elliptical arcs failed the gate by chance (it summed
  2,000 rounded steps against the rounding of one); it now counts whole
  turns, and the failing case is kept.
- Booleans of spheres and coaxial cones (design note 0012 §13, item
  3a), 2026-10-13, on Mesa lavapipe: fmt, clippy with `-D warnings`, 948
  tests and `cargo deny` pass, and the release time guards pass as CI
  runs them. Golden cases with all four operations, exact to 1e-9 in
  volume, with topology and origins: a ball resting in a plate, a ball
  cut by a box, a spherical pocket, two balls side by side and on one
  axis, a countersink, a cone on a rod, a pointed cone cut square, and a
  ball crossed twice by a cone. Planes across a cone off its axis are
  refused as conics (tested: a hyperbola, an ellipse, the apex). A fuzz
  of balls, cones, frustums, rods with round ends, tubes and stepped
  rods about one axis, with each other, with boxes at any angle, with
  balls anywhere and with slabs, at 50 times its default size (5,000
  pairs, 20,000 operations): every result valid, correct by the volume
  identities and by point membership, with its origins on its operands'
  surfaces; 56 refusals in 1,000 (1,096 tangent contacts, mostly on the
  grid, and 26 non-manifold). Its first runs found 5 correct results in
  40,000 operations refused as inside out: a thin ring a ball keeps
  outside a cone meshes coarsely to no volume, and validation took the
  shell's orientation from that mesh alone. It now decides as assembly
  does, by the shell's mass properties where the mesh cannot (kept as a
  golden case and a unit test). The other fuzzes at scale on the final
  code refuse exactly what they refused before: polyhedra at 100 times
  (47 of 120,000 operations, every result correct), planes and
  cylinders at 50 times (827 of 40,000) and at an angle (9,598 of
  20,000, quartics). New time guards (docs/PERFORMANCE.md, fourth pass)
  for dimpled and countersunk plates and a ball's mesh found two things,
  both fixed before this commit: the validation change at first made
  every thin flat body measured exactly (the 1024-sided prism's
  validation 10 ms to 87 ms), and meshing a curved face searched √n
  triangles for each point it added (a unit ball's mesh now 220–249 ms,
  from 364–411 ms). The curved fuzz's own checks cost 2.4 times its
  booleans, from about 8 (shared meshes; ray parity in a tree of boxes,
  checked against the winding number).
- Booleans of planes across cones (design note 0012 §15, item 3b),
  2026-10-15, on Mesa lavapipe: fmt, clippy with `-D warnings`, 980
  tests and `cargo deny` pass, and the release time guards pass as CI
  runs them. Golden cases with all four operations, exact to 1e-9 in
  volume against the closed-form areas of the elliptical, parabolic and
  hyperbolic segments cut from each cone, with topology and origins: a
  pointed cone under a tilted roof (a whole ellipse), a countersink
  leaning 15° into a plate, a cone cut along a ruling (a parabola) and a
  thousandth of a radian either side, a cone cut parallel to its axis (a
  hyperbola), a slight taper cut along its axis, a frustum cut through
  its apex (two rulings), and the cases review and the fuzzes found
  (0012 §15). A fuzz of tapered solids with boxes and slabs across them,
  through their apexes and nearly tangent along a ruling, at 50 times
  its size (20,000 operations): every result valid, correct by the
  volume identities and by point membership; 5,672 refused, every one a
  plane through or near a pointed cone's apex or nearly tangent along a
  ruling (220 of those on frustums). The other fuzzes at scale refuse
  what they refused before, except 63 operations of 60,000 (52 of
  40,000 among planes and cylinders, 11 of 20,000 at an angle) where an
  edge touches a face without crossing it: the new rule for touches
  refuses them, where the code before it had done them. Three
  adversarial reviews before this commit found what the tests had not,
  all fixed with tests: planes nearly tangent through a frustum's apex
  that gave two solids a nanometre apart (474 of an 8,000-case sweep)
  or inconsistent cuts (582); meshes of curved faces far past their
  chord tolerance (16 tolerances by a pointed cone's apex, 5% on a
  ball); a fuzz arm that never ran; and curve sampling that took
  tolerances that were not positive. Meshing was rebuilt to hold its
  tolerance, as estimated between the points it measures (design note
  0002 §4): on a grid of every triangle against
  closest points, a ball, a ring, S-splines, a nearly closed torus and a
  revolved cubic bent nearly to a crease are all within it (worst
  0.99997), where the old mesher strayed up to 2.4 times the tolerance
  (2.37, a nearly closed torus). On four benchmark bodies it costs 1.6
  to 2.4 times as long a triangle; profiling took the first version
  from 3.6 times main's time on a ball to 2.1 times (docs/PERFORMANCE.md,
  fifth pass), with new time guards for the tilted countersinks, the
  flattened pins and the creased revolve.
- Operands with conic edges on curved faces (design note 0012 §16 and
  §17, ADR 0023), 2026-10-15, on Mesa lavapipe: fmt, clippy with `-D
  warnings`, 1,011 tests and `cargo deny` pass, and the release time
  guards pass as CI runs them. Golden second cuts with all four
  operations, exact in volume, with topology and origins: a conic-edged
  body against a box far away, a rod cut by a tilted plane and again by
  a plane across its ellipse, a pointed cone under a ridge, the tilted
  countersink plate halved by its mirror plane and cut beside its hole,
  a flattened pin halved through its hyperbola's vertex and flattened
  again on its other side, cuts repeated in the same plane and with the
  same cone, and a cut leaving a void beside an ellipse. Far from the
  origin: a rod cut twice at the origin, 1 km and 30 km from it, within
  a millionth in volume, and refused as too far at 300 km and 3,000 km;
  a cone's hyperbola a millimetre off its axis cut again up to 30 km
  from the origin, within 2·10⁻⁹ in volume of the same cuts at the
  origin, and refused at 100 km; and cones 1 m to 400 km tall cut by a
  plane parallel to their axis (first cuts), within 10⁻¹² in volume
  over size cubed, one 600 km tall refused (the hyperbola's vertex too
  far away to compute within the tolerance, 0012 §14). A fuzz of second
  cuts at 50 times its size: on rods cut by a slab and again by a box,
  slab or rod,
  all 25,020 operations done, valid and correct by the volume identities
  and point membership; on cut cones 137 in 1,000 refused, nearly all
  at a pointed cone's apex or tangent contacts; in chains of three or
  four cuts 38 in 1,000, likewise. The fuzzes main ran refuse what they
  refused on main, but for two unions now done and correct. The fuzz
  found a fault in a first cut (a slab through a frustum's apex nearly
  along its seam failed with an internal error where a refusal
  belonged), now refused (0012 §17). Five adversarial reviews (0012
  §17), each finding confirmed by at least two of three skeptics or
  checked against the code: the first corrected the rules of §16; the
  second, third, fourth and fifth had six findings each, every code fix
  with a test that fails without it, nine of them false claims in the
  documentation and one an unfilled test count, corrected. The third review's rounding findings were settled by
  measurement against double-double references: a conic track's
  allowance is now 32 ε of
  its reach (was 10⁻¹³, which refused first cuts of cones 20 km tall
  that main did), a NURBS course on a plane carries a bound derived from
  its degree, size and weights (its floor in proportion to the weights'
  ratio times its whole size refused second cuts across uneven
  hyperbolas from 1 km), and every course's positional allowance counts
  the chart's size; tests against those references fail without each
  allowance, floor or term (ten mutations checked). The fourth found
  halving a conic track's chords stopped at four times its allowance,
  as other bounds' do, though its bound settles on the allowance itself:
  that could refuse first cuts of cones from about 235 km tall that
  more halving resolves (it refused one 400 km tall), and it now halves
  to within a thousandth of the allowance. The fifth found that meet
  accepted conics whose allowance came within that thousandth of a
  sixteenth of the tolerance, where the halving can stop over it: needle
  cones 879 km to 879.6 km long were refused as too far from the origin
  where meet should refuse them first; meet now refuses with two
  thousandths to spare, from about 878 km of reach. The validator's check of edges against a turned face's
  axis now takes
  closed forms first: validating the tilted countersunk plate fell from
  14.3 ms to 4.2 ms and cutting its holes in one boolean from 41.1 ms to
  24.3 ms, medians against main (docs/PERFORMANCE.md, sixth pass); a cut
  one at a time into conic-edged cones takes about twice as long as one
  into straight countersinks (16 take 197 ms against 102 ms). The scan
  that keeps the kernel's sources from reading the environment now
  leaves out modules built only for tests (the rounding tests read
  their sample sizes from it) and also catches variables read when the
  kernel is built (`crates/kernel/tests/sources.rs`).
- Contacts that do not cross on flat faces (design note 0012 §18 and
  §19, ADR 0024), 2026-10-16, on Mesa lavapipe: fmt, clippy with `-D
  warnings`, 1,072 tests and `cargo deny` pass, and the release time
  guards pass as CI runs them. 31 goldens (`tests/touching.rs`, golden
  11 in `tests/booleans.rs`) assert all six operations (union,
  intersection and difference, each in both orders) with counts, volumes, areas, origins, and refusal texts and places: a
  wedge on its edge and a cube on its corner on a plate from outside and
  inside, a groove's root, a notch, the slit a knife left on main,
  contacts on a face's boundary, many contacts in one face, millimetre
  and kilometre scales, and thirteen found by review and the fuzz. A
  contact fuzz (a random polyhedron with a box on its corner or a prism
  on its edge resting on one of its faces, from outside or inside, moved
  across the tolerance band in quarter steps): at a hundred times its
  run, all 264,000 operations answered (by the volume identities, an
  origins check and an oracle that finds results touching themselves)
  or refused, none failing; on main the same run, with the final
  oracle, fails 24,571 (20,457
  results that validate and touch themselves, 2,148 naming origins that
  do not bound their entities, 1,926 shell-volume kernel errors, 40
  others). The flat fuzz and the motions at a hundred times refuse what
  main refuses; the curved fuzzes at twenty times likewise, but for two
  pairs where a box's corner rests on a prism's flat face from inside,
  now answered but for the pinch, and pair 3493, where main returned a
  result touching itself, now refused. Every result is checked for a
  place within the tolerance of a flat face bounded by straight edges
  that it is not on; disabled one
  at a time, each of its three rules for sharp corners, and its measure
  of a point's distance, fails a test that pins it (a fourth rule,
  which none needed, was dropped). The slice costs 5.6% more
  instructions on a union of two 256-sided star prisms and 8.1% on 64
  pins cut one at a time; a union with a plate of 400 flush pins, which
  a first version of the check took quadratic time on, is guarded
  (95–113 ms, limit 1,000 ms) with a scaling check (docs/PERFORMANCE.md,
  seventh pass).
- Geometry for booleans of curved solids (design note 0012 §10.7, item
  1), 2026-10-13, on Mesa lavapipe: fmt, clippy with `-D warnings`,
  915 tests and `cargo deny` pass. Tested on its own; the boolean stages
  do not use it yet, so curved bodies are still refused. Planes,
  cylinders, cones and spheres are recognised exactly from the surfaces
  extrudes and revolves make (tori, NURBS and elliptic cylinders are
  not), with signed distances checked at 2,000 random points at known
  offsets (1e-13). An edge (line, arc, elliptical arc, conic, cubic
  spline) is split where a surface's implicit function turns, found as
  the roots, in Bernstein form, of its derivative along the edge's
  rational Bézier pieces (the NURBS crate now extracts them by
  blossoming), and at the pieces' joints. Against 3,000-point sampling of
  3,000 random edge–surface pairs: no crossing is missed, the implicit
  function is monotone between splits and turns at each, and every
  crossing is on the surface to 1e-12; edges lying on a surface (circles
  of cylinders, cones and spheres, an ellipse on a cylinder) have no
  spurious turns. Where two surfaces meet, in closed form: two planes, a
  plane and a cylinder (two lines or an ellipse), cylinders with
  parallel axes, and any two surfaces of revolution about one axis
  (circles where their meridians meet: a plane square to the axis,
  cylinders, cones, spheres centred on it). Over 4,800 random pairs every
  curve is on both surfaces to 1e-12 and every crossing found by
  sampling one surface is on a curve; touching, coincident and separate
  pairs are decided within the tolerance, and configurations within a
  sixteenth of the tolerance of a special one snap to it and stay within
  that sixteenth of both surfaces. The tests were checked against ten
  deliberately injected faults (dropped joints, no rounding floor, a
  dropped derivative term, a cone's slope off by 10%, a lost shift of a
  cone's reference point, a reversed bisection, a stretched ellipse,
  wrong circle heights, a shifted line, a wrong separation band); each
  fails a test. The cone's slope first passed, so the tests were
  strengthened: the implicit function is now checked directly, not only
  through the distance. A plane across a cone off its
  axis (conics) and pairs meeting in quartics are refused as unsupported.
- Near-coincidences in booleans, 2026-10-13, on Mesa lavapipe: fmt,
  clippy with `-D warnings`, 901 tests and `cargo deny` pass. Crossings
  within the tolerance of each other along an edge (an edge meeting a
  plane at a shallow angle) are now one point: a copy of a body moved by
  a quarter or a half of the tolerance is never refused (was 4 and 24 of
  800 operations), and at 0.9 tolerances 8 are (was 108); the test now
  requires none at half a tolerance or less. A run at 100 times the
  default (30,000 random pairs and 7,000 motions, 148,000 operations)
  found two older faults, both fixed and kept as tests that fail on the
  old code. Two sides of a prism meeting at a slight crease were taken
  for one plane, because the test meant to accept only rounding compared
  the cosine of the angle (1 − cos θ ≤ 1e-10, so angles up to 1.4e-5 rad)
  instead of its sine; a boolean then merged them into one face with a
  corner 4 µm off its plane, which the validator allowed, since it held
  planes to the 100 tolerances curved surfaces get (planes now 10; a
  result's vertices are at most 1.06 tolerances off their planes in the
  fuzz). The same mistake hid such a crease as a seam in the display and
  called a line revolved 1e-7 rad off its axis a cylinder; both fixed.
  And a strip 2.4 tolerances wide between the holes of a body and of its
  copy moved by three tolerances was classified on the other's face,
  because its point was within the tolerance of a corner of that face;
  that is now decided by parity alone, and the strip is refused. At 100
  times every result is now valid and consistent; 47 of 120,000
  operations on random pairs are refused, all non-manifold results.
- Second efficiency pass (docs/PERFORMANCE.md), 2026-10-13, on Mesa
  lavapipe: fmt, clippy with `-D warnings`, 897 tests and `cargo deny`
  pass, and the release time guards pass as CI runs them. Faces clear of
  the other body skip the point-inside and parity test; convex pieces
  get their point without a triangulation; convex flat faces are meshed
  without a Delaunay search. Measured by wall clock: piece points over the
  64 cuts 45 ms to 25 ms, a 256-sided union 48 ms to 43 ms, meshing a
  1024-sided prism 9 ms to 8 ms. The boolean fuzz at 20 times the
  default refuses exactly the same 2,191 operations for the same reasons
  as before. A patched incircle predicate (no buffer zeroing) matched the
  published crate on 36 million calls but saved only 3–6%, so it is not
  carried; the earlier note that the zeroing took 60% of meshing came
  from an instruction count and was wrong.
- Time guards and the first efficiency pass (docs/PERFORMANCE.md),
  2026-10-13, on Mesa lavapipe: fmt, clippy with `-D warnings`,
  893 tests and `cargo deny` pass. Release-mode time guards for
  extrude, mesh, validate, mass properties, booleans, sketch regions and
  regenerating a document (each failing at about ten times its measured
  time) and three scaling checks (eight times the size; growth bounded at
  24) pass; CI now runs them in release, and every job and long step has
  a timeout. The boolean fuzz at 20 times the default gives exactly the
  outcomes it gave before the indexes (6,000 pairs, 11 refusals; 1,400
  motions, 2,180 refusals). Meshes and boolean results are the same on
  every run: the triangulation library is used as a copy patched to
  insert deferred points in order (ADR 0022), and a new golden test
  repeats meshes and booleans with fresh hash seeds (it fails with the
  published library). Union of two 1024-sided prisms
  5.4 s to under 0.6 s; extruding an 8192-gon 315 ms to 92 ms; a plate
  with 1,600 holes 1.35 s to 0.21 s. The full run also caught a flaw in
  a NURBS property test (not in the code): mirroring a parameter as
  a + b − u rounds by an ulp, which a steep random curve turns into
  1e-11 of distance; the test now allows its speed times that rounding.
- Booleans of polyhedral solids (design note 0012 §9), 2026-10-13, on
  Mesa lavapipe: fmt, clippy with `-D warnings`, 881 tests and
  `cargo deny` pass. Golden cases check volume, area, topology counts,
  genus and origins (overlapping cubes, merged and flush boxes, a T, a
  slot and a through hole, a void, disjoint, identical and touching
  bodies, at 1 mm and 1 km too). Every result is validated and its
  origins checked against the operands' geometry (planes, edge lines,
  vertex positions). Random pairs check the three volume identities to
  1e-9 relative: at 20 times the default run, 6,000 pairs
  with every result valid and 11 of 24,000 operations refused, each a
  non-manifold result. A copy of a random body moved by 0.25 to 10
  tolerances (1,400 pairs) gives no invalid or inconsistent result;
  refusals per 800 operations are 4 and 24 at a quarter and a half
  tolerance, 108 to 800 between 0.9 and 3, and 156 at 10. Writing the
  tolerance test found five inconsistencies (a point near a corner on
  two edges, points close in the face's plane but not in space, pieces
  classified within the tolerance of the other surface); each is now a
  named refusal. CI runs both random tests at ten times the default size
  and prints the rates.
- Picking and references in the application (ADR 0019), 2026-10-13, on
  Mesa lavapipe: fmt, clippy with `-D warnings`, 867 tests and `cargo
  deny` pass. The new UI tests (4) drive the real app offscreen: the
  pointer picks a corner vertex before its edges and an edge before its
  faces, each named as expected, with the hovered face lit; a click on
  the plate's top face and Sketch makes a sketch on that face (its frame
  at the face's height, cancelling removes it), while the cylinder's
  side is refused; an extrude of the plate's region up to the cylinder's
  picked top refuses the curved side with a message and ends at 30 mm
  with the exact volume; nothing is picked where there are no parts.
  New features are numbered one past the highest name of their type
  (deleting one no longer leads to a duplicate name). Screenshots 33 to
  37 show hovering, an edge and vertex, a sketch on a face and an
  extrude up to a face.
- References from features (ADR 0019), model side, 2026-10-13, on Mesa
  lavapipe: fmt, clippy with `-D warnings`, 863 tests and `cargo deny`
  pass. Every arrangement test region's profile now also checks the
  sketch edge reported for each curve (each curve spans its edge, runs
  with the selection on its left, every boundary edge reported once);
  every origin test body checks its face planes against its
  outward-wound mesh; the reference tests (10) cover box faces, caps and
  sides of sketch extrudes and inline cylinders, a face in two pieces,
  edits that move, redraw or remove the referenced face, suppression,
  reordering and deletion, and a document saved before the change. The
  text and sketch tests, which now name every part they make, all named
  without a warning.
- Origin tables and persistent names, 2026-10-12, on Mesa lavapipe:
  fmt, clippy with `-D warnings`, 851 tests and `cargo deny` pass. The
  origin tables pass
  their completeness, topology and geometry checks on every box, extrude
  and revolve test body; the naming tests (14) check names against
  geometry and every resolution rule. Writing them found that one sketch
  line split by a vertex gave two faces of one plane, which the display
  showed as one face (the edge between them is a seam): logical faces now
  follow the seams, so selection matches what is seen.
- Feature-size-adaptive tessellation and exact profile checks,
  2026-10-12, on Mesa lavapipe: fmt, clippy with `-D warnings`, 830 tests
  and `cargo deny` pass. The random sketch tests no longer exempt any
  region from extrusion (they had excused the kernel's refusals of
  slivers, narrow notches and cusps). Making that pass found and fixed:
  mass properties wrong by up to 1.6e-5 (a sliver) and about 1e-7
  (ordinary bodies) wherever merged end points leave edges off their
  faces by up to five tolerances (now exact to rounding: design note
  0003, "Tolerant edges"); a thin lens tessellated as an inside-out
  triangle; and a crescent between tangent circles whose cusp refinement
  chased to rounding (now matched samples, 0002 §5). Drafts of the
  refinement also exhausted memory on an ellipse (a closure-point
  projection) and ran away at two sharp corners (166 s and 74 s for one
  extrude); the final version takes milliseconds there. Each case is a
  regression test, checked to fail without its fix.
  Stress runs at raised case counts: polygons and circles 2,000, grids
  2,000, scale invariance 1,000, splines 500, garbage input 20,000, all
  passed. The lattice test's 40,000-case run failed after 24,130 cases:
  two arcs tangent to a line at its end, all three leaving the vertex
  towards −x, made Euler's formula fail at scale 10³ only. The vertex's
  half-edges formed one run of nearly equal directions with a single
  gap, and the run's middle direction was taken half a turn past its
  first angle, as for a run with no gap at all, so it pointed away from
  the curves and their angles seen from it fell on both sides of ±π.
  Fixed (a run closed by a gap spans its own angles) and kept as a
  regression test, checked to fail without the fix (the rerun is in the
  entry above). Documented refusals in the failed run: lattice inputs
  perturbed by about a tolerance or more refused as unresolvable (198),
  profile loops below the builder's resolution (14), slivers whose mass
  cannot show 1e-10 (57, computed to 1e-6) and pinched profiles (30);
  none was a wrong result.
- Curvature combs and region robustness, 2026-10-12, on Mesa lavapipe:
  fmt, clippy with `-D warnings`, 822 tests and `cargo deny` pass. The
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
  thinner than the tolerance that are artefacts of it are absorbed
  (before the face checks, including faces whose area is rounding noise),
  points of a curve inside a merged vertex count as at it, the crossing
  of two nearly parallel lines far from the origin is no longer lost to
  rounding, and the profile builder refuses micron-sized loops by name.
  Successive runs passed 7,533, 13,908, 2,868, 8,360, 14,673, 39,619 and
  21,360 cases before the next defect.
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
- Tessellation (and so extrude's validation) still fails with an error,
  never wrong geometry, where two boundary curves of a face cross or
  touch within the tolerance in its parameter space (refinement stops at
  segments a tolerance long) or where a face's boundary would need more
  than eight times its initial samples plus 100,000. The chord tolerance
  is estimated, not bounded: where a surface bends faster than a quartic
  between the points measured, the mesh can stray past it, by an amount
  nothing bounds (no
  test case does; a rigorous bound from the surface's derivatives is not
  built). A ball meshes to about 2.6 times the triangles an even mesh
  would need at the same tolerance (126,656 at the defaults): 4.5 times
  near the poles, where its parameters stretch the triangulation, which
  is Delaunay in them, and 1.8 times at the equator, from the triangles'
  shape and refinement stopping short of the tolerance. In the random sketch
  tests every region now extrudes, cusps and slivers included (design
  note 0002 §5). Mass properties are exact to the rounding of the
  coordinates a body's points are computed from: a body far from the
  origin, or a sliver on an arc of large radius, gets its integrals to
  about 1e-13 of those coordinates times its size to one dimension less,
  not 1e-13 of its own size (a lens 0.3 µm wide on arcs of 3 m radius:
  its volume to 1e-9).
- Region detection refuses a sketch with a long sliver thinner than the
  tolerance that has three or more sides (design note 0006 §4.2 absorbs
  the other kinds), naming where; and the profile builder, which merges
  end points within ten tolerances, refuses a region that small as a whole
  (a loop within twenty tolerances of a point). Both are errors, never
  wrong regions; in the random lattice test they occur only for curves
  moved by about a tolerance or more.
- Solver conflict detection is local: from a start far from any solution it
  can call a sketch conflicting although another branch would satisfy it.
- Booleans handle solids bounded by planes and circular cylinders,
  spheres and cones where every meeting is a circle (design note 0012
  §13), and planes across cones (§15). They refuse as not supported
  yet: spheres and cones with cylinders or cones off their axes, and two
  cylinders at an angle (quartics); tangent contacts where a face at
  the contact is curved, faces touching without crossing (a box's face
  resting on a cylinder's side, a slot tangent to a hole) or an edge or
  a corner resting on a face: 324 of 16,000 random operations in the
  fuzz of planes and cylinders at 20 times (292 faces touching, 32 an
  edge or a corner), all its refusals but five non-manifold results,
  and 4 of 8,000 at an angle (an edge or a corner resting on a flat
  face of such a body, the faces there all flat, is decided since §19);
  planes nearly tangent along a cone's
  ruling through or near its apex, where the cone stays within the
  tolerance of the plane out to the faces (the cone fuzz tilts them
  10⁻⁴ to 10⁻⁷ rad from tangent: 220 of its 20,000 operations on
  frustums). An operand's NURBS edge on a sphere, or of degree above
  five on any curved face, is refused (no operation makes one yet:
  planes meet spheres in circles, and the conics cuts leave are
  quadratic). The fuzz of second cuts at 50 times refuses none of
  25,020 operations on cut rods and 137 in 1,000 on cut cones, four
  fifths of those through or near a pointed cone's apex and most of the
  rest tangent contacts, and 38 in 1,000 in chains of three or four
  cuts, likewise (design note 0012 §17).
- Contacts that do not cross on flat faces (design note 0012 §19) are
  decided, but much of the tolerance band is refused: the contact fuzz
  at a hundred times its run refuses 80,627 of 264,000 operations (5,573
  of 24,000 with the contact exact, 12,095 with the body moved a
  tolerance against the face's normal), where main refused 53,168 and
  failed 24,571. A fin
  under a tolerance thick, which main returns where an edge of one body
  lies on a face's plane beside the line where its own face crosses
  that face, past reflex corners, and the face holds the strip between
  them (§19's known fault, with more cases found by review: the face's
  edge on that line, a hole's edge, a notch's floor), is refused (§20,
  goldens 32 and 33); where the face does not hold the strip (B resting
  on its edge, or over a hole) the contact is answered as before. The
  check does not itself decide a corner of the face within the
  tolerance beside that line under the strip: probes of that were
  refused by other checks, but no proof is claimed. Within three quarters of a tolerance of the exact contact, 704
  of 83,065 results differ in their counts from the exact contact's
  (mostly an edge left split), 673 of them as on main: results are not
  yet canonical under motions within the tolerance. The oracle that
  finds contacts the booleans missed does not run on the curved fuzzes
  yet: its test of a point in a curved face is not exact.
- Far from the origin: a boolean is refused as too far from the origin
  for the tolerance where the rounding of the bound on a curved edge's
  course keeps a chord that could hide a contact wider than a sixteenth
  of the tolerance (design note 0012 §16.2). At the default tolerance:
  from about 36 km from the origin for the conics cuts leave on
  cylinders and cones (a rod cut twice is cut at 30 km and refused at
  300 km; a cone's hyperbola cut again is cut at 30 km and refused at
  100 km); for the conic where a plane cuts a cone, which applies to
  first cuts too, never: its vertex is refused first as too far away to
  compute within the tolerance (design note 0012 §14: from about 878 km
  of the vertex's coordinates, its distance from where the faces meet
  and that region's size; a cone 400 km tall at the origin is cut by a
  plane parallel to its axis, one 600 km tall refused; needle cones
  whose hyperbolas' vertices are at the origin are cut to 877.9 km long
  and refused from there); and, by the formula for its allowance, from about 160 km
  (degree five) to 310 km (a conic) for a small NURBS edge on a plane
  whose chart is anchored near the world's origin (farther when it is
  anchored near the edge), less for a large edge whose weights are very
  uneven.
- At a cone's apex: a plane through a pointed cone's apex, or whose
  curve comes within 128 tolerances of the axis there, is refused where
  the apex lies near where faces of both bodies can be, which takes in
  nearly every section of a pointed cone through its axis (a half
  section of a drill point: the result's edges would meet at the cone's
  singular point). Planes through the apex of a frustum, whose faces stop
  short of it, work, its axis leaning or not. A box whose face runs
  within 100 tolerances (10 µm) of a turned solid's axis, but not on
  it, is refused where an edge of the result would pass the centre of a
  flat end within that ("a turned face's centre"; tested from 0.5 to
  5 µm). A ball halved through its axis at the angle of its faces'
  seams is refused as edges or points meeting within the tolerance of
  each other at its poles; at other angles the halves are exact (both
  in `a_ball_halved_through_its_axis`).
  A corner of one body resting exactly where a face's test point falls
  (the centre of a cylinder's end face) is refused, though the result
  exists (once in 10,000 random pairs).
  No feature uses booleans yet: each feature still makes its own part.
  Near-coincidences the tolerance cannot decide are refused, not
  resolved: a copy of a body moved by half a tolerance or less is always
  combined, but 0.5–1% of operations are refused at 0.9 tolerances (mostly a
  point within the tolerance of two edges near a sharp corner), and most
  at about a tolerance or more, where the bodies leave slivers thinner
  than it (design note 0012 §9 gives the rates). The result's vertices sit at the clusters' points,
  so a vertex can be up to twice the tolerance from where an exact
  construction would put it; the surface-area identity is not tested
  yet.
- References to model geometry (ADR 0019): up to a face works only for
  a planar face parallel to the sketch; others fail with a message
  until surface intersection exists. Picked edges and vertices select
  and highlight but no feature uses them yet (fillets, chamfers and
  sketch projection will). Picking runs on the display mesh, so a
  point picked on a curved edge or face lies within the display chord
  of it (vertices, straight edges and planar faces are exact). A sketch on a face that tilts through 45° between facing up
  and facing sideways turns with the frame rule (ADR 0019 §2). Hints are
  checked against a mesh within 0.1 mm, so ties between pieces of an
  entity finer than that fail as ambiguous. A cap reference survives
  sketch edits while the feature's stored region key finds its region;
  once a region is selected anew with other boundary curves (a region
  cut in two, say), references to its caps fail rather than guess.
  Splits and merges by later operations come with booleans. No
  incremental regeneration: every edit regenerates the whole studio.
- Known slow spots (docs/PERFORMANCE.md, measured by wall clock): flat
  faces are triangulated by constrained Delaunay, whose exact incircle
  test is slow on points sharing a circle (round holes, rectangles); the
  flat faces take two thirds of meshing a plate with round holes, mostly
  in that test. Validation meshes every result; small pieces that are
  not convex get their interior point from a triangulation; curved
  segments carry about eleven extra chords at each end, so splitting a
  face with 100 round holes takes 9 ms; every boolean reads both whole
  bodies, so its time grows with the bodies rather than with where they
  meet (drilling a hole into a plate with 99 takes 40 ms); and each
  boolean groups a body's faces by surface with a linear search over
  the groups found so far, faces times distinct surfaces (4.8% of the
  instructions of a union with a plate of 900 flush pins, growing 14.5
  times from 100 to 900 pins). All are within their budgets.
- Non-manifold topology deferred to Phase 3.
- No STEP or other exchange formats beyond STL (Phase 4; see docs/FORMATS.md).

## Next steps

1. Check the new title bar on GNOME (on Pat's desktop).
2. Phase 3, in the order the design notes give (0012 §12.4 item 3:
   balls and coaxial cones, item 3a, §13; planes across cones, item 3b,
   §15, operands with conic edges on cylinders and cones, §17, and
   contacts that do not cross on flat faces, §19 and §20, are done).
   Next: the rest of the contacts' slice 1
   (§18.6 steps 5 to 8: the cut-pair filter, edges and corners on curved
   faces, 76 of the curved fuzzes' refusals at twenty times, the Apart
   band, the contact time guards and mutation checks), then slice 2,
   surfaces that touch (§18.3), 1,304 of those refusals; then the apex
   left refused (a plane through a pointed cone's apex where both
   bodies reach it: edges meeting at the cone's singular point), naming
   through booleans (design note 0013 §10 step 4), and extrude and
   revolve that add, remove and intersect. Alongside, absorb the boolean
   refusals left within the tolerance (0.5–1% of operations at 0.9
   tolerances, mostly a point within the tolerance of two edges near a
   sharp corner). The rest of the Phase 2 sketcher (project/convert,
   intersection curves, pierce) needs Phase 3 faces.
3. Efficiency, alongside Phase 3 (docs/PERFORMANCE.md, known slow spots):
   flat faces triangulated from orientation tests alone (no incircle
   test), which speeds meshing and validation; a cheaper exact stage of
   the incircle test, which curved faces' meshes reach on points of one
   circle (18% of a ball's mesh); faces grouped by surface through an
   index of the surfaces, not a linear search (on main too, §19's
   "Cost"); booleans that work near
   where the bodies meet (validation and classification of untouched
   faces carried over); then lower the budgets to match.
