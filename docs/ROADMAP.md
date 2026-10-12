# OpenShape roadmap

Items are ticked only when implemented **and tested**. Each phase ends
building, tested and runnable, with this file and `STATUS.md` updated.
Feature-area checklists below the phases list the full functional scope; each
item is delivered in the phase noted.

## Phase 0 — Foundations
- [x] ADRs: kernel (0001), UI toolkit (0002), persistent naming (0003), op-log & sync (0004), native format (0005), units & tolerance (0006), rendering (0007), periodic face splitting (0008), dependency licensing (0009); later: solver linear algebra (0010), sketch document model (0011), UI look (0012), rollback bar (0013)
- [x] Cargo workspace, CLAUDE.md, ARCHITECTURE, ROADMAP, STATUS
- [x] CI: fmt, clippy `-D warnings`, tests, `cargo-deny` (licenses/advisories/bans/sources)
- [x] Op-log skeleton: stable IDs, causal ordering, deterministic replay, undo/redo as ops, history-at-heads, CBOR container, convergence property test
- [x] Document model with a Part Studio feature list and a reducer over the op log
- [x] wgpu viewport (Vulkan/RADV primary, GL fallback), verified offscreen on lavapipe
- [x] Kernel-made box on screen, built through the op log, with undo/redo
- [x] Orbit/zoom/pan camera

## Phase 1 — Kernel core
- [x] NURBS curves: evaluation, derivatives, knot insertion, split, closest point, exact arcs, interpolation
- [x] NURBS surfaces: evaluation, derivatives, normals, closest point, revolution/extrusion construction
- [x] Analytic curves (line, circle) and surfaces (plane, extrusion, revolution with cylinder/cone/sphere/torus classification)
- [x] B-rep topology (solid/shell/face/loop/coedge/edge/vertex) with validation (closed loops, manifold edge use, edges on surfaces, Euler–Poincaré, outward orientation)
- [x] Tessellation with chord/angle tolerance, watertight across faces, singular points handled
- [x] Extrude (blind, symmetric; profiles with holes and multiple regions; lines, arcs, circles, splines)
- [x] Revolve (full and partial; profiles touching the axis; cones, spheres, tori)
- [x] Mass properties (volume, area, centroid, inertia) by exact boundary integration (1e-9 relative on the golden corpus)
- [x] Golden corpus: volume, area, topology counts, validity, watertightness
- [x] STL export (ASCII and binary)
- [ ] Non-manifold topology (deferred to Phase 3 with booleans)

## Phase 2 — Sketcher
- [x] Sketch data model in core-model; sketch mode in the UI (ADR 0011)
- [x] Constraint solver (Newton/LM with DOF analysis, over/under-constraint diagnostics; design note 0005)
- [x] Property tests over random constraint systems; fuzzing (constructive and garbage-input proptests)
- [x] Region detection from arbitrary sketch curves (intersections) feeding extrude/revolve (design note 0006)
- [ ] Sketch tools and constraints (see §2 below; core set, editing tools, ellipses, conics, transform, curvature and text done; project/convert, intersection curves and pierce need Phase 3 faces)
- [x] Application shell in the common parametric-CAD layout with our own look (ADR 0012), driven end to end by offscreen UI tests

## Phase 3 — Booleans and core features
- [ ] Curve/curve, curve/surface, surface/surface intersection (design note first)
- [ ] Robust booleans; fuzzing; `|A∪B| = |A|+|B|−|A∩B|` property test; failure-rate metric in CI
- [ ] Incremental regeneration engine (`regen`): dependency graph, caching, cancellation, error propagation
- [ ] Persistent naming (`naming`) with reorder/suppress/pattern-count tests
- [ ] Hole, patterns, mirror, construction planes
- [ ] Milestone: a real bracket survives upstream edits with no broken references

## Phase 4 — Fillets, chamfers, shell, draft, offsets, STEP
- [ ] Constant/variable/conic fillets, vertex blends; chamfers
- [ ] Offsets, shell, thicken, draft
- [ ] STEP AP203/214/242 import/export; healing and sewing

## Phase 5 — Documents and history
- [ ] Native format file I/O and open spec (docs/NATIVE_FORMAT.md)
- [ ] Versions, branches, history restore UI, compare
- [ ] Autosave and crash recovery from the op log; solo-app polish

## Phase 6 — Assemblies
- [ ] Mate solver, mates and relations, limits; standard content; interference

## Phase 7 — Drawings and 2D formats
- [ ] Drawing engine, views, annotations, tables; PDF/SVG/DXF/DWG

## Phase 8 — Advanced modeling
- [ ] Sweep, loft, surfacing, sheet metal, frames, variables, configurations

## Phase 9 — Custom features, API, plugins
- [ ] `script` language with stdlib over the kernel; API; WASM plugins

## Phase 10 — Server/client
- [ ] Sync, real-time collaboration, merge UI, PDM/release

## Phase 11 — Simulation
## Phase 12 — CAM and ECAD
## Phase 13 — Render studio and remaining formats

---

## Full functional scope

### 1. Documents and data management
- [ ] Tabs: Part Studios ✔ (add, rename, delete, switch), Assemblies, Drawings, Variable Studios, scripts, attached files (PDF, images, other)
- [ ] Folders, labels, search, recents, thumbnails, public/private (server)
- [x] Undo/Redo (op-log based)
- [ ] Full edit history with restore-to-any-point (log support ✔, UI pending)
- [ ] Named immutable Versions
- [ ] Branches
- [ ] Feature-level merge with conflict UI
- [ ] Visual and feature-tree compare
- [ ] Linked documents pinned to versions, update notifications
- [ ] PDM: part numbers, revisions, release workflow with approvers, custom properties, where-used, obsolescence

### 2. Part Studio
Sketch tools
- [ ] Lines ✔, arcs (3-point ✔, tangent ✔, center ✔), circles ✔, ellipses and elliptical arcs ✔ (trim/extend/split of them pending), splines (through-point ✔, control-point ✔, smooth closed), conics ✔ (trim/extend/split of them pending), rectangles (corner ✔, center ✔, aligned ✔), polygons ✔, slots ✔, points ✔
- [ ] Construction geometry ✔, text ✔ (one line, six bundled faces), image tracing
- [ ] Use/project/convert, offset ✔, trim ✔, extend ✔, split ✔, sketch fillet ✔, mirror ✔, linear/circular sketch patterns ✔, transform (move/rotate/scale, copy) ✔, intersection curves

Sketch constraints
- [ ] Coincident ✔, concentric ✔, parallel ✔, perpendicular ✔, horizontal ✔, vertical ✔, tangent ✔, equal ✔, midpoint ✔, normal ✔, pierce, symmetric ✔, fix ✔, curvature ✔
- [x] Dimensions: linear (aligned, horizontal, vertical), angular, radial, diameter, a conic's ρ; typed with units
- [x] Inference while drawing (horizontal/vertical, parallel, perpendicular, tangent; point, curve and origin snaps); colour-coded DOF state

Features
- [ ] Extrude: blind ✔, symmetric ✔, from picked sketch regions ✔, up-to, through-all, draft, second direction, thin
- [ ] Revolve (full/partial ✔, about a sketch line ✔; up-to and two-direction pending)
- [ ] Sweep: guides, twist, profile control
- [ ] Loft: guides, end conditions, connections
- [ ] Thicken, enclose
- [ ] Fillet: constant, variable, conic, full-round
- [ ] Chamfer
- [ ] Hole: simple, counterbore, countersink, tapped; ANSI/ISO tables
- [ ] Draft, shell, rib
- [ ] Thread, cosmetic thread
- [ ] Booleans: union, subtract, intersect
- [ ] Split, transform
- [ ] Direct edits: delete face, move face, replace face
- [ ] Offset surface
- [ ] Patterns: linear, circular, curve (features, faces, parts)
- [ ] Mirror
- [ ] Construction: planes, mate connectors, helix
- [ ] Curves: 3D fit spline, bridging, projected, composite, isocline

Surfacing
- [ ] Fill, boundary surface, ruled surface, extend, trim/split
- [ ] Analysis: zebra, curvature combs, draft analysis, deviation

Sheet metal
- [ ] From sketch/part/thicken; flange, hem, tab, reliefs, corner break, form tools, make joint
- [ ] Live flat pattern, bend table, K-factor and bend radius

Frames / weldments
- [ ] Frames from sketch/path with ANSI/ISO profile library; trim/cap, gussets, cut list

Variables and configurations
- [ ] Variables, Variable Studio, expressions with units (units parsing ✔)
- [ ] Configurations: input tables, checkboxes, lists

Inspection and appearance
- [ ] Measure: distance, angle, area, volume
- [ ] Mass properties with materials library (kernel mass properties ✔)
- [ ] Section view, interference check, properties panel
- [ ] Per-part and per-face appearance

### 3. Assemblies
- [ ] Instances from same/linked documents and standard content
- [ ] Mates: fastened, revolute, slider, planar, cylindrical, pin-slot, ball, parallel, tangent, width
- [ ] Relations: gear, rack-and-pinion, screw, linear
- [ ] Mate connectors, groups, limits; replicate; linear/circular/feature patterns
- [ ] DOF animation, named positions, exploded views
- [ ] Section view, interference detection, mass properties
- [ ] Rigid/flexible sub-assemblies; in-context design
- [ ] Parametric fastener library (ANSI/ISO), auto-sizing to holes

### 4. Drawings
- [ ] Sheets, templates, title blocks (ANSI, ISO)
- [ ] Projected, section (full/half/aligned/offset), detail, auxiliary, broken, broken-out, isometric views
- [ ] Associative dimensions, GD&T, datums, surface finish, weld symbols, notes, centerlines, center marks, hatching
- [ ] BOM tables + balloons, hole tables, bend tables, revision tables
- [ ] Export PDF, SVG, DXF, DWG

### 5. Rendering
- [ ] Shaded ✔, wireframe ✔, hidden-line display
- [ ] Silhouette edges, transparency, section caps
- [ ] Named views (presets ✔), perspective ✔ / orthographic ✔, view cube ✔
- [ ] GPU picking with selection filters; box/lasso select
- [ ] LOD and instancing; 60 fps at 1,000+ parts on RADV
- [ ] Render studio: PBR, HDRI, decals, progressive path tracer (wgpu compute)

### 6. Simulation
- [ ] Linear static, modal, thermal steady-state; pure-Rust tet meshing; loads/constraints; mates → contact
- [ ] Stress/displacement/safety-factor plots; animated modes; textbook benchmark validation

### 7. CAM
- [ ] Facing, pocket, contour, drilling, adaptive; 3-axis parallel/scallop; 4th-axis indexing
- [ ] Tool library, stock, WCS; material-removal simulation with gouge/collision checks
- [ ] Post-processors: Fanuc, Haas, LinuxCNC, Grbl, Mach3/4

### 8. ECAD
- [ ] IDF and IDX import, board 3D model, push-back of mechanical changes

### 9. Collaboration
- [ ] Real-time multi-user editing on the op log (convergence ✔ in oplog tests); presence; follow mode
- [ ] Comments and @mentions pinned to geometry
- [ ] Sharing permissions, teams, activity feed, audit log, admin tools
- [ ] Docker and systemd packaging, TLS; offline merge on reconnect

### 10. Import/export
- [ ] STEP AP203/214/242 (+PMI where feasible)
- [ ] IGES
- [ ] STL ASCII/binary export ✔ — import pending
- [ ] 3MF, OBJ, glTF/GLB, PLY, VRML, AMF
- [ ] Parasolid text `.x_t` (only if clean-room from public docs is feasible — see FORMATS.md)
- [ ] DXF/DWG import/export; SVG/PDF export; images
- [ ] Native format spec (draft ✔ in NATIVE_FORMAT.md; file I/O pending)
- [x] Proprietary formats documented in FORMATS.md with legal routes
- [ ] Fuzz target + real-world corpus for every importer

### 11. Extensibility
- [ ] Custom-feature language editor (highlighting, diagnostics, debug output), shareable library
- [ ] Automation API (queries, feature create/edit, exports, webhooks)
- [ ] WASM plugins: features, panels, importers

### 12. UX baseline
- [ ] Feature tree: rollback bar ✔, suppress ✔, reorder (move up/down ✔, drag pending), edit-in-place ✔, failure indicators ✔
- [ ] Context menus ✔, remappable shortcuts, command palette ✔
- [ ] 3D mouse via spacenavd
- [ ] Light/dark themes, HiDPI, multi-window
- [ ] Autosave and crash recovery from the op log
