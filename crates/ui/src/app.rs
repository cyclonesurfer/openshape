//! The egui desktop front end: the application shell of ADR 0012.
//!
//! [`OpenDatumApp`] holds the view state around a [`Model`]. Each frame
//! [`OpenDatumApp::draw`] lays out the top bar, the feature toolbar, the
//! element tabs, the feature list and the graphics area
//! ([`crate::shell`], [`crate::feature_list`], [`crate::graphics`]), then the
//! floating dialogs ([`crate::floating`]). Widgets do not change the model
//! while they are drawn; they queue [`Action`]s that are applied at the end
//! of the frame. The same `draw` runs in the window (eframe) and offscreen
//! ([`crate::screenshots`]), so tests and screenshots exercise the real
//! application.

use std::collections::{BTreeSet, HashMap};

use eframe::egui;
use egui::{Key, Modifiers, Pos2, Rect};
use opendatum_core_model::naming::Use;
use opendatum_core_model::{DVec3, ElementId, Feature, FeatureId, InsertAt, ModelOp, PlaneRef};
use opendatum_naming::Pick;
use opendatum_render::{Camera, DisplayMode, Projection};
use opendatum_ui_kit::palette::CommandPalette;
use opendatum_ui_kit::view_cube::CubeTarget;
use opendatum_ui_kit::Theme;

use crate::commands::{Command, Mode};
use crate::dialogs::FeatureDialog;
use crate::entity_pick::EntityHit;
use crate::graphics::Viewport;
use crate::planes::StdPlane;
use crate::profile_pick::{ProfilePick, SweepKind};
use crate::sketch_mode::{normal_view, SketchDisplay, SketchEditor, Tool};
use crate::{FeaturePreview, Model};

/// An axis-aligned box `(min, max)` in model space.
pub(crate) type Bounds = ([f32; 3], [f32; 3]);

/// The smallest box holding both (either may be absent).
pub(crate) fn union(a: Option<Bounds>, b: Option<Bounds>) -> Option<Bounds> {
    match (a, b) {
        (Some((a0, a1)), Some((b0, b1))) => Some((
            std::array::from_fn(|i| a0[i].min(b0[i])),
            std::array::from_fn(|i| a1[i].max(b1[i])),
        )),
        (a, None) => a,
        (None, b) => b,
    }
}

/// What is selected outside sketch mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Selection {
    /// A feature in the feature list.
    Feature(FeatureId),
    /// A part, by the feature that made it.
    Part(FeatureId),
    /// A face, edge or vertex of a part, picked in the view (ADR 0019).
    Entity(EntityHit),
    /// A default plane.
    Plane(StdPlane),
    /// The origin.
    Origin,
}

/// What the file prompt is for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FileAction {
    /// Open a document.
    Open,
    /// Save the document.
    Save,
    /// Export STL.
    ExportStl,
}

/// An inline rename in progress.
#[derive(Clone, Debug, PartialEq)]
pub enum Renaming {
    /// A feature, with the text typed so far.
    Feature(FeatureId, String),
    /// An element (tab), with the text typed so far.
    Element(ElementId, String),
}

/// A change requested by a widget while drawing; applied after the frame's
/// widgets are drawn.
#[derive(Clone, Debug)]
pub enum Action {
    /// Runs a command.
    Run(Command),
    /// Changes the selection.
    Select(Option<Selection>),
    /// Opens a feature's editor.
    EditFeature(FeatureId),
    /// Applies a model edit.
    Apply(ModelOp),
    /// Moves the rollback bar.
    SetRollback(InsertAt),
    /// Starts an inline rename.
    StartRename(Renaming),
    /// Commits the inline rename.
    CommitRename,
    /// Abandons the inline rename.
    CancelRename,
    /// Switches to an element (tab).
    Activate(ElementId),
    /// Deletes an element (tab).
    DeleteElement(ElementId),
    /// Shows or hides a default plane.
    TogglePlane(StdPlane),
    /// Shows or hides a part.
    TogglePart(FeatureId),
    /// Shows every part.
    ShowAllParts,
    /// Fits the view to one part.
    ZoomTo(FeatureId),
    /// Turns the view to a view-cube target.
    SetView(CubeTarget),
}

/// Screen rectangles of the shell's controls in the last frame, for tests
/// and automation that drive the application with pointer events.
#[derive(Clone, Debug)]
pub struct ShellLayout {
    /// Toolbar and top-bar buttons.
    pub commands: Vec<(Command, Rect)>,
    /// Feature-list rows of features.
    pub feature_rows: Vec<(FeatureId, Rect)>,
    /// Feature-list rows of the default planes.
    pub plane_rows: Vec<(StdPlane, Rect)>,
    /// Rows of the parts list.
    pub part_rows: Vec<(FeatureId, Rect)>,
    /// The rollback bar.
    pub rollback_bar: Option<Rect>,
    /// Element tabs.
    pub tabs: Vec<(ElementId, Rect)>,
    /// The add-element button.
    pub add_tab: Option<Rect>,
    /// The command search field.
    pub search: Option<Rect>,
    /// The graphics area.
    pub graphics: Rect,
    /// The view cube.
    pub view_cube: Option<Rect>,
    /// The open floating dialog.
    pub dialog: Option<Rect>,
}

impl Default for ShellLayout {
    fn default() -> Self {
        ShellLayout {
            commands: Vec::new(),
            feature_rows: Vec::new(),
            plane_rows: Vec::new(),
            part_rows: Vec::new(),
            rollback_bar: None,
            tabs: Vec::new(),
            add_tab: None,
            search: None,
            graphics: Rect::NOTHING,
            view_cube: None,
            dialog: None,
        }
    }
}

impl ShellLayout {
    /// The button of a command.
    pub fn command(&self, c: Command) -> Option<Rect> {
        self.commands.iter().find(|(k, _)| *k == c).map(|(_, r)| *r)
    }

    /// The row of a feature.
    pub fn feature_row(&self, f: FeatureId) -> Option<Rect> {
        self.feature_rows
            .iter()
            .find(|(k, _)| *k == f)
            .map(|(_, r)| *r)
    }

    /// The tab of an element.
    pub fn tab(&self, e: ElementId) -> Option<Rect> {
        self.tabs.iter().find(|(k, _)| *k == e).map(|(_, r)| *r)
    }
}

/// A part's measured properties, cached per regeneration.
#[derive(Clone, Debug)]
pub(crate) struct PartProperties {
    pub generation: u64,
    pub feature: FeatureId,
    pub lines: Vec<(String, String)>,
}

/// The application.
pub struct OpenDatumApp {
    pub(crate) model: Model,
    pub(crate) theme: Theme,
    pub(crate) camera: Camera,
    pub(crate) viewport: Option<Viewport>,
    pub(crate) selection: Option<Selection>,
    pub(crate) dialog: Option<FeatureDialog>,
    pub(crate) pick: Option<ProfilePick>,
    pub(crate) sketch: Option<SketchEditor>,
    /// The undo depth when sketch mode was entered (cancel returns to it).
    pub(crate) sketch_undo_depth: usize,
    pub(crate) renaming: Option<Renaming>,
    pub(crate) file_prompt: Option<(FileAction, String)>,
    pub(crate) palette: Option<CommandPalette>,
    /// The plane menu of the Sketch command is open.
    pub(crate) plane_menu: bool,
    /// The last tool used from each toolbar group, by the group's first
    /// tool.
    pub(crate) tool_choice: HashMap<Tool, Tool>,
    pub(crate) message: Option<(String, f64)>,
    pub(crate) adapter: String,
    pub(crate) display_mode: DisplayMode,
    pub(crate) fit_pending: bool,
    pub(crate) show_sketches: bool,
    pub(crate) planes_visible: [bool; 3],
    pub(crate) hidden_parts: BTreeSet<FeatureId>,
    /// A plane row under the pointer (highlighted in the view).
    pub(crate) hover_plane: Option<StdPlane>,
    /// The face, edge or vertex under the pointer (highlighted in the
    /// view).
    pub(crate) hover_entity: Option<EntityHit>,
    /// The model generation a selected face, edge or vertex was picked in:
    /// a regeneration may renumber the part's entities, so it is cleared
    /// then.
    pub(crate) entity_generation: u64,
    /// [`OpenDatumApp::sketch_bounds`] for a model generation.
    pub(crate) sketch_bounds_cache: std::cell::Cell<Option<(u64, Option<Bounds>)>>,
    pub(crate) properties: Option<PartProperties>,
    pub(crate) layout: ShellLayout,
    pub(crate) doc_name: String,
    pub(crate) actions: Vec<Action>,
    /// Seconds since start, from egui's input.
    pub(crate) time: f64,
    /// The last pointer position over the graphics area.
    pub(crate) pointer: Option<Pos2>,
}

impl OpenDatumApp {
    /// Creates the app in an eframe window.
    pub fn new(cc: &eframe::CreationContext<'_>, model: Model) -> Self {
        Self::with_render_state(&cc.egui_ctx, cc.wgpu_render_state.as_ref(), model)
    }

    /// Creates the app on `ctx` with the wgpu device in `render_state` (the
    /// 3D view is unavailable without one). Installs the theme.
    pub fn with_render_state(
        ctx: &egui::Context,
        render_state: Option<&egui_wgpu::RenderState>,
        model: Model,
    ) -> Self {
        let theme = Theme::light();
        theme.apply(ctx);
        let (viewport, adapter) = match render_state {
            Some(rs) => {
                let info = rs.adapter.get_info();
                log::info!(
                    "opendatum: window uses {} ({:?}, {:?}, driver {} {})",
                    info.name,
                    info.backend,
                    info.device_type,
                    info.driver,
                    info.driver_info
                );
                (
                    Some(Viewport::new(rs, &theme)),
                    format!("{} ({:?})", info.name, info.backend),
                )
            }
            None => (None, "no GPU device".to_string()),
        };
        Self {
            model,
            theme,
            camera: Camera::isometric(),
            viewport,
            selection: None,
            dialog: None,
            pick: None,
            sketch: None,
            sketch_undo_depth: 0,
            renaming: None,
            file_prompt: None,
            palette: None,
            plane_menu: false,
            tool_choice: HashMap::new(),
            message: None,
            adapter,
            display_mode: DisplayMode::ShadedWithEdges,
            fit_pending: true,
            show_sketches: true,
            planes_visible: [true; 3],
            hidden_parts: BTreeSet::new(),
            hover_plane: None,
            hover_entity: None,
            entity_generation: 0,
            sketch_bounds_cache: std::cell::Cell::new(None),
            properties: None,
            layout: ShellLayout::default(),
            doc_name: "Untitled document".into(),
            actions: Vec::new(),
            time: 0.0,
            pointer: None,
        }
    }

    /// The model.
    pub fn model(&self) -> &Model {
        &self.model
    }

    /// The model, for scripted edits (tests, screenshots).
    pub fn model_mut(&mut self) -> &mut Model {
        &mut self.model
    }

    /// The camera.
    pub fn camera_mut(&mut self) -> &mut Camera {
        &mut self.camera
    }

    /// Where the shell's controls were drawn in the last frame.
    pub fn layout(&self) -> &ShellLayout {
        &self.layout
    }

    /// The current selection.
    pub fn selection(&self) -> Option<Selection> {
        self.selection
    }

    /// Selects something (or nothing).
    pub fn select(&mut self, s: Option<Selection>) {
        self.selection = s;
        self.entity_generation = self.model.generation;
    }

    /// The face, edge or vertex under the pointer in the view.
    pub fn hovered_entity(&self) -> Option<EntityHit> {
        self.hover_entity
    }

    /// The sketch editor, in sketch mode.
    pub fn sketch_editor(&self) -> Option<&SketchEditor> {
        self.sketch.as_ref()
    }

    /// The sketch editor, in sketch mode, for scripted edits.
    pub fn sketch_editor_mut(&mut self) -> Option<&mut SketchEditor> {
        self.sketch.as_mut()
    }

    /// The open profile picker (extrude or revolve), if any.
    pub fn profile_pick_mut(&mut self) -> Option<&mut ProfilePick> {
        self.pick.as_mut()
    }

    /// The open template dialog, if any.
    pub fn feature_dialog(&self) -> Option<&FeatureDialog> {
        self.dialog.as_ref()
    }

    /// The open template dialog, for scripted edits.
    pub fn feature_dialog_mut(&mut self) -> Option<&mut FeatureDialog> {
        self.dialog.as_mut()
    }

    /// True while the command palette is open.
    pub fn palette_open(&self) -> bool {
        self.palette.is_some()
    }

    /// The command palette, while open.
    pub fn palette_mut(&mut self) -> Option<&mut CommandPalette> {
        self.palette.as_mut()
    }

    /// The message shown in the graphics area, if any.
    pub fn message(&self) -> Option<&str> {
        self.message.as_ref().map(|(m, _)| m.as_str())
    }

    /// The document name shown in the top bar.
    pub fn set_document_name(&mut self, name: &str) {
        self.doc_name = name.to_string();
    }

    /// Shows a short message in the graphics area.
    pub fn notify(&mut self, text: impl Into<String>) {
        self.message = Some((text.into(), self.time));
    }

    /// Which commands apply now.
    pub fn mode(&self) -> Mode {
        if self.sketch.is_some() {
            Mode::Sketch
        } else {
            Mode::PartStudio
        }
    }

    /// True if `c` can run now.
    pub fn enabled(&self, c: Command) -> bool {
        if !c.available_in(self.mode()) {
            return false;
        }
        match c {
            Command::Undo => self.model.session().can_undo(),
            Command::Redo => self.model.session().can_redo(),
            Command::Extrude | Command::Revolve => !self.model.sketches().is_empty(),
            Command::RollToEnd => self.model.studio().rollback != InsertAt::End,
            Command::Constrain(t) => self
                .sketch
                .as_ref()
                .is_some_and(|ed| ed.can_constrain(&self.model, t)),
            Command::DeleteSelection => self
                .sketch
                .as_ref()
                .is_some_and(|ed| !ed.selection.is_empty()),
            Command::CurvatureComb => self
                .sketch
                .as_ref()
                .is_some_and(|ed| !ed.comb_candidates(&self.model).is_empty()),
            Command::ResetLabels => self
                .sketch
                .as_ref()
                .is_some_and(|ed| !ed.labels_to_reset(&self.model).is_empty()),
            _ => true,
        }
    }

    /// True if `c` is a toggled-on state (a selected tool, the display mode).
    pub fn active(&self, c: Command) -> bool {
        match c {
            Command::Tool(t) => self.sketch.as_ref().is_some_and(|ed| ed.tool == t),
            Command::Construction => self.sketch.as_ref().is_some_and(|ed| ed.construction),
            Command::CurvatureComb => self
                .sketch
                .as_ref()
                .is_some_and(|ed| ed.combs_shown(&self.model)),
            Command::Display(m) => self.display_mode == m,
            Command::Perspective => {
                matches!(self.camera.projection, Projection::Perspective { .. })
            }
            Command::Orthographic => matches!(self.camera.projection, Projection::Orthographic),
            Command::ToggleSketches => self.show_sketches,
            Command::TogglePlanes => self.planes_visible.iter().any(|v| *v),
            Command::Extrude => self
                .pick
                .as_ref()
                .is_some_and(|p| p.kind == SweepKind::Extrude),
            Command::Revolve => self
                .pick
                .as_ref()
                .is_some_and(|p| p.kind == SweepKind::Revolve),
            Command::Template(t) => self.dialog.as_ref().is_some_and(|d| d.template == t),
            _ => false,
        }
    }

    /// Runs a command (does nothing when it is disabled).
    pub fn run(&mut self, c: Command) {
        if !self.enabled(c) {
            return;
        }
        match c {
            Command::NewDocument => {
                self.close_editors();
                self.model = Model::new();
                self.selection = None;
                self.hidden_parts.clear();
                self.doc_name = "Untitled document".into();
                self.fit_pending = true;
            }
            Command::OpenDocument => {
                self.file_prompt = Some((FileAction::Open, "part.odatum".into()));
            }
            Command::SaveDocument => {
                self.file_prompt = Some((FileAction::Save, "part.odatum".into()));
            }
            Command::ExportStl => {
                self.file_prompt = Some((FileAction::ExportStl, "part.stl".into()));
            }
            Command::Undo => {
                if !self.model.undo() {
                    self.notify("Nothing to undo");
                }
            }
            Command::Redo => {
                if !self.model.redo() {
                    self.notify("Nothing to redo");
                }
            }
            Command::NewPartStudio => {
                self.close_editors();
                self.model.add_part_studio();
                self.selection = None;
                self.fit_pending = true;
            }
            Command::CommandSearch => self.palette = Some(CommandPalette::new()),
            Command::Sketch => match self.selection {
                Some(Selection::Plane(p)) => self.new_sketch(p),
                Some(Selection::Entity(hit)) if matches!(hit.pick, Pick::Face(_)) => {
                    self.new_sketch_on_face(&hit)
                }
                _ => self.plane_menu = true,
            },
            Command::SketchOn(p) => self.new_sketch(p),
            Command::Extrude | Command::Revolve => {
                let kind = if c == Command::Extrude {
                    SweepKind::Extrude
                } else {
                    SweepKind::Revolve
                };
                self.close_editors();
                self.pick = Some(ProfilePick::new(kind));
            }
            Command::Template(t) => {
                self.close_editors();
                self.dialog = Some(FeatureDialog::new(t));
            }
            Command::RollToEnd => self.model.set_rollback(InsertAt::End),
            Command::FitView => self.fit(),
            Command::View(t) => self.set_view(t),
            Command::Display(m) => self.display_mode = m,
            Command::Perspective => {
                self.camera.projection = Projection::Perspective {
                    fov_y_rad: opendatum_render::DEFAULT_FOV_Y_RAD,
                }
            }
            Command::Orthographic => self.camera.projection = Projection::Orthographic,
            Command::ToggleSketches => self.show_sketches = !self.show_sketches,
            Command::TogglePlanes => {
                let show = !self.planes_visible.iter().any(|v| *v);
                self.planes_visible = [show; 3];
            }
            Command::Tool(t) => {
                if let Some(ed) = self.sketch.as_mut() {
                    ed.set_tool(t);
                }
                self.tool_choice.insert(t.group()[0], t);
            }
            Command::Construction => {
                if let Some(ed) = self.sketch.as_mut() {
                    ed.toggle_construction(&mut self.model);
                }
            }
            Command::Constrain(t) => {
                if let Some(ed) = self.sketch.as_mut() {
                    ed.constrain(&mut self.model, t);
                }
            }
            Command::CurvatureComb => {
                if let Some(ed) = self.sketch.as_mut() {
                    ed.toggle_combs(&self.model);
                }
            }
            Command::Dimension => {
                let at = self
                    .pointer
                    .unwrap_or_else(|| self.layout.graphics.center());
                if let Some(ed) = self.sketch.as_mut() {
                    ed.dimension(&self.model, at);
                }
            }
            Command::ResetLabels => {
                if let Some(ed) = self.sketch.as_mut() {
                    ed.reset_labels(&mut self.model);
                }
            }
            Command::DeleteSelection => {
                if let Some(ed) = self.sketch.as_mut() {
                    ed.delete_selection(&mut self.model);
                }
            }
            Command::FinishSketch => self.sketch = None,
            Command::CancelSketch => {
                self.sketch = None;
                self.model.undo_to(self.sketch_undo_depth);
            }
        }
    }

    /// Closes sketch mode (keeping its edits), the profile picker and the
    /// template dialog.
    fn close_editors(&mut self) {
        self.sketch = None;
        self.pick = None;
        self.dialog = None;
        self.plane_menu = false;
        self.sync_edit_point();
    }

    /// The feature an open editor changes: the sketch in sketch mode or the
    /// feature an open dialog edits.
    pub fn editing(&self) -> Option<FeatureId> {
        self.sketch
            .as_ref()
            .map(|e| e.feature)
            .or(self.pick.as_ref().and_then(|p| p.editing))
            .or(self.dialog.as_ref().and_then(|d| d.editing))
    }

    /// What the open feature dialog would make, while its values build a
    /// feature.
    pub fn dialog_preview(&self) -> Option<FeaturePreview> {
        let (editing, built) = match (&self.dialog, &self.pick) {
            (Some(d), _) => (d.editing, d.build()),
            (None, Some(p)) => (p.editing, p.build()),
            (None, None) => return None,
        };
        let feature = built.ok()?;
        Some(match editing {
            Some(id) => FeaturePreview::Edit(id, feature),
            None => FeaturePreview::New(feature),
        })
    }

    /// Makes the model show the history up to the edited feature, with the
    /// open dialog's feature as if applied; returns whether it changed.
    fn sync_edit_point(&mut self) -> bool {
        let point = self.editing();
        let preview = self.dialog_preview();
        self.model.set_view(point, preview)
    }

    /// Enters sketch mode on a sketch feature, looking straight at its plane.
    pub fn enter_sketch(&mut self, feature: FeatureId) {
        let Some(r) = self.model.sketch_result(feature) else {
            self.notify(if self.model.studio().is_rolled_back(feature) {
                "Move the rollback bar below the sketch to edit it"
            } else {
                "Unsuppress the sketch to edit it"
            });
            return;
        };
        let frame = r.frame;
        let d = SketchDisplay::from_result(r);
        self.close_editors();
        self.selection = Some(Selection::Feature(feature));
        self.sketch = Some(SketchEditor::new(feature));
        self.sketch_undo_depth = self.model.session().undo_depth();
        self.sync_edit_point();
        let (yaw, pitch) = normal_view(&frame);
        self.camera.yaw = yaw;
        self.camera.pitch = pitch;
        let pts: Vec<DVec3> = d
            .points()
            .map(|(_, p)| frame.point2(p))
            .chain(
                d.curves
                    .iter()
                    .flat_map(|c| c.polyline.iter().map(|p| frame.point2(*p))),
            )
            .collect();
        if let Some((lo, hi)) = bounds(&pts) {
            self.camera.fit(lo, hi);
        } else if self.model.meshes().is_empty() {
            // An empty sketch in an empty document: show about 100 mm.
            self.camera.target = frame.origin;
            self.camera.distance = 0.25;
        }
    }

    /// Starts a new sketch on a default plane. Cancelling sketch mode
    /// removes it again.
    pub fn new_sketch(&mut self, plane: StdPlane) {
        self.close_editors();
        let depth = self.model.session().undo_depth();
        let id = self.model.add_sketch(plane.plane_ref());
        self.enter_sketch(id);
        self.sketch_undo_depth = depth;
    }

    /// Starts a new sketch on a planar face of a part (ADR 0019); a curved
    /// face or one that cannot be referenced is refused with a message.
    /// Cancelling sketch mode removes the sketch again.
    pub fn new_sketch_on_face(&mut self, hit: &EntityHit) {
        let Pick::Face(face) = hit.pick else {
            self.notify("Select a planar face or a plane to sketch on");
            return;
        };
        if self.model.face_plane(hit.feature, face).is_none() {
            self.notify("Sketches go on planar faces: select a flat face or a plane");
            return;
        }
        let reference = match self.model.reference(hit, Use::Carrier) {
            Ok(r) => r,
            Err(e) => {
                self.notify(e);
                return;
            }
        };
        self.close_editors();
        let depth = self.model.session().undo_depth();
        let id = self.model.add_sketch(PlaneRef::Face(Box::new(reference)));
        self.enter_sketch(id);
        self.sketch_undo_depth = depth;
    }

    /// Opens the editor for a feature: sketch mode, the profile picker or a
    /// parameter dialog.
    pub fn edit_feature(&mut self, id: FeatureId) {
        let Some(f) = self.model.feature(id).cloned() else {
            return;
        };
        if self.model.studio().is_rolled_back(id) {
            self.notify("Move the rollback bar below the feature to edit it");
            return;
        }
        if matches!(f.feature, Feature::Sketch(_)) {
            self.enter_sketch(id);
        } else if let Some(p) = ProfilePick::edit(&f) {
            self.close_editors();
            self.selection = Some(Selection::Feature(id));
            self.pick = Some(p);
        } else {
            self.close_editors();
            self.selection = Some(Selection::Feature(id));
            self.dialog = FeatureDialog::edit(&f);
            if self.dialog.is_none() {
                self.notify("This feature has no editor yet");
            }
        }
        self.sync_edit_point();
    }

    /// Fits the view to the visible parts, sketches and default planes.
    pub fn fit(&mut self) {
        let sketches = self.show_sketches.then(|| self.sketch_bounds()).flatten();
        let h = self.plane_half();
        let planes = StdPlane::ALL
            .into_iter()
            .filter(|p| self.sketch.is_none() && self.planes_visible[p.index()])
            .flat_map(|p| crate::planes::corners(p, h))
            .map(|c| {
                let c = c.as_vec3().to_array();
                (c, c)
            })
            .reduce(|a, b| union(Some(a), Some(b)).expect("both present"));
        let content = union(union(self.visible_bounds(), sketches), planes);
        let (lo, hi) = content.unwrap_or_else(|| {
            let h = crate::planes::half_size(None, None) as f32;
            ([-h, -h, -h], [h, h, h])
        });
        self.camera.fit(lo, hi);
    }

    /// Bounds of the visible parts.
    pub(crate) fn visible_bounds(&self) -> Option<Bounds> {
        self.model
            .meshes()
            .iter()
            .zip(self.model.mesh_owners())
            .filter(|(_, f)| !self.hidden_parts.contains(f))
            .filter_map(|(m, _)| m.bounds())
            .reduce(|a, b| union(Some(a), Some(b)).expect("both present"))
    }

    /// Bounds of the regenerated sketches' curves in model space.
    pub(crate) fn sketch_bounds(&self) -> Option<Bounds> {
        if let Some((generation, b)) = self.sketch_bounds_cache.get() {
            if generation == self.model.generation {
                return b;
            }
        }
        let mut b = None;
        for s in self.model.sketches() {
            let f = s.frame;
            for c in SketchDisplay::from_result(s).curves {
                for p in c.polyline {
                    let w = (f.origin + f.x * p.x + f.y * p.y).as_vec3().to_array();
                    b = union(b, Some((w, w)));
                }
            }
        }
        self.sketch_bounds_cache
            .set(Some((self.model.generation, b)));
        b
    }

    /// Turns the view to a view-cube target; the home target also fits.
    pub fn set_view(&mut self, t: CubeTarget) {
        let (yaw, pitch) = t.yaw_pitch();
        self.camera.yaw = yaw;
        self.camera.pitch = pitch;
        if t == CubeTarget::HOME {
            self.fit();
        }
    }

    /// Moves the rollback bar so the features in `features[..gap]` are
    /// active.
    pub fn roll_to_gap(&mut self, gap: usize) {
        let st = self.model.studio();
        let at = if gap >= st.features.len() {
            InsertAt::End
        } else if gap == 0 {
            InsertAt::Start
        } else {
            InsertAt::After(st.features[gap - 1].id)
        };
        self.model.set_rollback(at);
    }

    /// The ops that move a feature one row up or down in the feature list,
    /// where the rollback bar counts as a row: crossing the bar moves the
    /// bar instead of the feature (ADR 0013).
    pub fn move_feature_ops(&self, id: FeatureId, up: bool) -> Option<ModelOp> {
        let st = self.model.studio();
        let element = self.model.studio_id();
        let i = st.index_of(id)?;
        let b = st.rollback_index();
        let n = st.features.len();
        let f = |k: usize| st.features[k].id;
        let swap = |k: usize| ModelOp::MoveFeature {
            element,
            feature: f(k),
            to: InsertAt::After(f(k + 1)),
        };
        let bar_after = |k: Option<usize>| ModelOp::SetRollback {
            element,
            at: k.map_or(InsertAt::Start, |k| InsertAt::After(f(k))),
        };
        if up {
            if i == b && b < n {
                Some(bar_after(Some(i)))
            } else if i > 0 {
                Some(swap(i - 1))
            } else {
                None
            }
        } else if i + 1 == b {
            Some(bar_after(i.checked_sub(1)))
        } else if i + 1 < n {
            Some(swap(i))
        } else {
            None
        }
    }

    fn apply_actions(&mut self) {
        for action in std::mem::take(&mut self.actions) {
            match action {
                Action::Run(c) => self.run(c),
                Action::Select(s) => {
                    self.selection = s;
                    self.entity_generation = self.model.generation;
                    self.plane_menu = false;
                }
                Action::EditFeature(f) => self.edit_feature(f),
                Action::Apply(op) => {
                    if let Err(reason) = self.model.apply(op) {
                        self.notify(reason);
                    }
                }
                Action::SetRollback(at) => self.model.set_rollback(at),
                Action::StartRename(r) => self.renaming = Some(r),
                Action::CommitRename => match self.renaming.take() {
                    Some(Renaming::Feature(id, name)) if !name.trim().is_empty() => {
                        let element = self.model.studio_id();
                        if self
                            .model
                            .feature(id)
                            .is_some_and(|f| f.name != name.trim())
                        {
                            let r = self.model.apply(ModelOp::RenameFeature {
                                element,
                                feature: id,
                                name: name.trim().to_string(),
                            });
                            if let Err(reason) = r {
                                self.notify(reason);
                            }
                        }
                    }
                    Some(Renaming::Element(id, name))
                        if !name.trim().is_empty()
                            && self
                                .model
                                .elements()
                                .iter()
                                .any(|e| e.id == id && e.name != name.trim()) =>
                    {
                        if let Err(reason) = self.model.rename_element(id, name.trim()) {
                            self.notify(reason);
                        }
                    }
                    _ => {}
                },
                Action::CancelRename => self.renaming = None,
                Action::Activate(e) => {
                    if e != self.model.studio_id() {
                        self.close_editors();
                        self.selection = None;
                        self.model.set_active(e);
                        self.fit_pending = true;
                    }
                }
                Action::DeleteElement(e) => {
                    if let Err(m) = self.model.delete_element(e) {
                        self.notify(m);
                    } else {
                        self.close_editors();
                        self.selection = None;
                        self.fit_pending = true;
                    }
                }
                Action::TogglePlane(p) => {
                    self.planes_visible[p.index()] = !self.planes_visible[p.index()]
                }
                Action::TogglePart(f) => {
                    if !self.hidden_parts.remove(&f) {
                        self.hidden_parts.insert(f);
                    }
                }
                Action::ShowAllParts => self.hidden_parts.clear(),
                Action::ZoomTo(f) => {
                    let b = self
                        .model
                        .meshes()
                        .iter()
                        .zip(self.model.mesh_owners())
                        .find(|(_, o)| **o == f)
                        .and_then(|(m, _)| m.bounds());
                    if let Some((lo, hi)) = b {
                        self.camera.fit(lo, hi);
                    }
                }
                Action::SetView(t) => self.set_view(t),
            }
        }
    }

    /// Keyboard shortcuts of the commands available in this mode, plus the
    /// keys sketch mode handles itself.
    fn shortcuts(&mut self, ctx: &egui::Context) {
        if ctx.egui_wants_keyboard_input() || self.palette.is_some() {
            return;
        }
        let mut commands: Vec<Command> = Command::all()
            .into_iter()
            .chain([Command::CommandSearch])
            .filter(|c| c.available_in(self.mode()) && c.shortcut().is_some())
            .collect();
        // Shortcuts with more modifiers first: egui matches Ctrl+Z
        // "logically" against Ctrl+Shift+Z too.
        commands.sort_by_key(|c| {
            let m = c.shortcut().map_or(Modifiers::NONE, |s| s.modifiers);
            std::cmp::Reverse(m.shift as u8 + m.alt as u8 + m.command as u8)
        });
        let mut run = Vec::new();
        ctx.input_mut(|i| {
            for c in &commands {
                if let Some(s) = c.shortcut() {
                    if i.consume_shortcut(&s) {
                        run.push(*c);
                    }
                }
            }
            if i.consume_key(Modifiers::COMMAND, Key::Y) {
                run.push(Command::Redo);
            }
        });
        let sketching = self.sketch.is_some();
        // Enter finishes a sketch; elsewhere it is left for the open
        // dialog, which accepts on it. Escape with a menu open is left for
        // the menu, which closes on it, and does nothing else.
        let menu_open = egui::Popup::is_any_open(ctx);
        let (escape, enter, backspace) = ctx.input_mut(|i| {
            (
                !menu_open && i.consume_key(Modifiers::NONE, Key::Escape),
                sketching && i.consume_key(Modifiers::NONE, Key::Enter),
                i.consume_key(Modifiers::NONE, Key::Backspace),
            )
        });
        if backspace && sketching {
            run.push(Command::DeleteSelection);
        }
        for c in run {
            self.run(c);
        }
        if let Some(ed) = self.sketch.as_mut() {
            if escape {
                ed.escape();
            }
            if enter {
                ed.finish(&mut self.model);
            }
        } else if escape {
            if self.plane_menu {
                self.plane_menu = false;
            } else if self.pick.is_some() || self.dialog.is_some() {
                self.pick = None;
                self.dialog = None;
            } else {
                self.selection = None;
            }
        }
        // Delete removes the selected feature outside sketch mode.
        if !sketching && ctx.input_mut(|i| i.consume_key(Modifiers::NONE, Key::Delete)) {
            if let Some(Selection::Feature(f)) = self.selection {
                let element = self.model.studio_id();
                if let Err(reason) = self.model.apply(ModelOp::DeleteFeature {
                    element,
                    feature: f,
                }) {
                    self.notify(reason);
                }
                self.selection = None;
            }
        }
    }

    /// Draws one frame of the application into `ui` (the whole window).
    /// `render_state` is the wgpu device the 3D view renders with.
    pub fn draw(&mut self, ui: &mut egui::Ui, render_state: Option<&egui_wgpu::RenderState>) {
        let ctx = ui.ctx().clone();
        self.time = ctx.input(|i| i.time);
        self.layout = ShellLayout::default();
        crate::shell::clear_menu_rows(&ctx);
        self.sync_edit_point();
        // Leave sketch mode if undo or a delete removed, rolled back or
        // suppressed the sketch.
        if let Some(ed) = &self.sketch {
            if self.model.sketch_result(ed.feature).is_none() {
                self.sketch = None;
                self.notify("The sketch was removed, rolled back or suppressed");
            }
        }
        if let Some(Selection::Feature(f) | Selection::Part(f)) = self.selection {
            if self.model.feature(f).is_none() {
                self.selection = None;
            }
        }
        if matches!(self.selection, Some(Selection::Entity(_)))
            && self.entity_generation != self.model.generation
        {
            self.selection = None;
        }
        self.shortcuts(&ctx);
        self.top_bar(ui);
        self.toolbar(ui);
        self.tab_bar(ui);
        self.feature_list(ui);
        self.graphics(ui, render_state);
        // The palette first: it takes the keys it uses while open.
        self.command_palette(&ctx);
        self.floating_dialogs(&ctx);
        self.apply_actions();
        // The view was drawn before the dialogs changed: draw it again.
        if self.sync_edit_point() {
            ctx.request_repaint();
        }
    }
}

/// Bounds of world points as renderer arrays.
fn bounds(pts: &[DVec3]) -> Option<([f32; 3], [f32; 3])> {
    let first = *pts.first()?;
    let (lo, hi) = pts
        .iter()
        .fold((first, first), |(lo, hi), p| (lo.min(*p), hi.max(*p)));
    Some((lo.as_vec3().to_array(), hi.as_vec3().to_array()))
}

impl eframe::App for OpenDatumApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let rs = frame.wgpu_render_state().cloned();
        self.draw(ui, rs.as_ref());
    }
}
