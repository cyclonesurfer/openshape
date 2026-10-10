//! The egui desktop front end.

use eframe::egui;
use egui::{Color32, Key, KeyboardShortcut, Modifiers, PointerButton, Sense};
use openshape_core_model::sketch::DistanceKind;
use openshape_core_model::{DVec3, Feature, FeatureId, InsertAt, ModelOp, PlaneRef};
use openshape_features::FeatureStatus;
use openshape_kernel::{Frame, Kernel};
use openshape_render::{wgpu, Camera, DisplayMode, Projection, SceneRenderer, ViewPreset};

use crate::dialogs::{FeatureDialog, FieldKind, Template};
use crate::profile_pick::{ProfilePick, SweepKind};
use crate::sketch_mode::{normal_view, ConstraintTool, ScreenMap, SketchDisplay, SketchEditor, Tool};
use crate::sketch_paint;
use crate::Model;

/// Offscreen target the 3D view renders into; shown by egui as a texture.
struct Viewport {
    scene: SceneRenderer,
    texture: Option<(wgpu::Texture, wgpu::TextureView)>,
    size: (u32, u32),
    texture_id: Option<egui::TextureId>,
    uploaded_generation: u64,
}

const VIEW_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;

/// What the file prompt is for.
#[derive(Clone, Copy, PartialEq, Eq)]
enum FileAction {
    Open,
    Save,
    ExportStl,
}

/// The application.
pub struct OpenShapeApp {
    model: Model,
    camera: Camera,
    viewport: Option<Viewport>,
    selected: Option<FeatureId>,
    dialog: Option<FeatureDialog>,
    rename: Option<(FeatureId, String)>,
    file_prompt: Option<(FileAction, String)>,
    message: String,
    adapter: String,
    display_mode: DisplayMode,
    fit_pending: bool,
    /// Sketch mode, when a sketch is being edited.
    sketch: Option<SketchEditor>,
    /// An extrude or revolve picking its profile from sketch regions.
    pick: Option<ProfilePick>,
    show_sketches: bool,
    /// Where the 3D view was last drawn (screen points).
    view_rect: egui::Rect,
}

impl OpenShapeApp {
    /// Creates the app; sets up GPU resources on eframe's wgpu device.
    pub fn new(cc: &eframe::CreationContext<'_>, model: Model) -> Self {
        let (viewport, adapter) = match cc.wgpu_render_state.as_ref() {
            Some(rs) => {
                let info = rs.adapter.get_info();
                (
                    Some(Viewport {
                        scene: SceneRenderer::new(&rs.device, VIEW_FORMAT),
                        texture: None,
                        size: (0, 0),
                        texture_id: None,
                        uploaded_generation: u64::MAX,
                    }),
                    format!("{} ({:?})", info.name, info.backend),
                )
            }
            None => (None, "no wgpu device".to_string()),
        };
        Self {
            model,
            camera: Camera::isometric(),
            viewport,
            selected: None,
            dialog: None,
            rename: None,
            file_prompt: None,
            message: String::new(),
            adapter,
            display_mode: DisplayMode::ShadedWithEdges,
            fit_pending: true,
            sketch: None,
            pick: None,
            show_sketches: true,
            view_rect: egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(800.0, 600.0)),
        }
    }

    /// Enters sketch mode on a sketch feature, looking straight at its plane.
    fn enter_sketch(&mut self, feature: FeatureId) {
        let Some(r) = self.model.sketch_result(feature) else {
            self.message = "Unsuppress the sketch to edit it".into();
            return;
        };
        let frame = r.frame;
        let d = SketchDisplay::from_result(r);
        self.dialog = None;
        self.pick = None;
        self.selected = Some(feature);
        self.sketch = Some(SketchEditor::new(feature));
        let (yaw, pitch) = normal_view(&frame);
        self.camera.yaw = yaw;
        self.camera.pitch = pitch;
        let pts: Vec<DVec3> = d
            .points()
            .map(|(_, p)| frame.point2(p))
            .chain(d.curves.iter().flat_map(|c| c.polyline.iter().map(|p| frame.point2(*p))))
            .collect();
        if let Some((lo, hi)) = bounds(&pts) {
            self.camera.fit(lo, hi);
        } else if self.model.meshes().is_empty() {
            // An empty sketch in an empty document: show about 100 mm.
            self.camera.target = frame.origin;
            self.camera.distance = 0.25;
        }
    }

    /// Starts a new sketch on a standard plane.
    fn new_sketch(&mut self, plane: PlaneRef) {
        let id = self.model.add_sketch(plane);
        self.enter_sketch(id);
    }

    /// Leaves sketch mode.
    fn finish_sketch(&mut self) {
        self.sketch = None;
    }

    /// Opens the editor for a feature: sketch mode, the profile picker or a
    /// parameter dialog.
    fn edit_feature(&mut self, id: FeatureId) {
        let Some(f) = self.model.feature(id).cloned() else {
            return;
        };
        if matches!(f.feature, Feature::Sketch(_)) {
            self.enter_sketch(id);
        } else if let Some(p) = ProfilePick::edit(&f) {
            self.sketch = None;
            self.pick = Some(p);
        } else {
            self.dialog = FeatureDialog::edit(&f);
            if self.dialog.is_none() {
                self.message = "This feature has no editor yet".into();
            }
        }
    }

    fn fit(&mut self) {
        if let Some((lo, hi)) = self.model.bounds() {
            self.camera.fit(lo, hi);
        }
    }

    fn undo(&mut self) {
        if !self.model.undo() {
            self.message = "Nothing to undo".into();
        }
    }

    fn redo(&mut self) {
        if !self.model.redo() {
            self.message = "Nothing to redo".into();
        }
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        if ctx.egui_wants_keyboard_input() {
            return;
        }
        if self.sketch.is_some() {
            self.sketch_shortcuts(ctx);
        }
        let redo_a = KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::Z);
        let redo_b = KeyboardShortcut::new(Modifiers::COMMAND, Key::Y);
        let undo = KeyboardShortcut::new(Modifiers::COMMAND, Key::Z);
        let (r, u, f) = ctx.input_mut(|i| {
            let r = i.consume_shortcut(&redo_a) || i.consume_shortcut(&redo_b);
            let u = i.consume_shortcut(&undo);
            let f = i.consume_key(Modifiers::NONE, Key::F);
            (r, u, f)
        });
        if r {
            self.redo();
        } else if u {
            self.undo();
        }
        if f {
            self.fit();
        }
    }

    fn sketch_shortcuts(&mut self, ctx: &egui::Context) {
        let Some(ed) = self.sketch.as_mut() else {
            return;
        };
        let pointer = ctx
            .input(|i| i.pointer.hover_pos())
            .unwrap_or(self.view_rect.center());
        let pressed = |k: Key| ctx.input_mut(|i| i.consume_key(Modifiers::NONE, k));
        if pressed(Key::Escape) {
            ed.escape();
        }
        if pressed(Key::Enter) {
            ed.finish(&mut self.model);
        }
        if pressed(Key::Delete) || pressed(Key::Backspace) {
            ed.delete_selection(&mut self.model);
        }
        if pressed(Key::D) {
            ed.dimension(&self.model, pointer);
        }
        if pressed(Key::Q) {
            ed.toggle_construction(&mut self.model);
        }
        for t in Tool::ALL {
            if let Some(k) = t.key() {
                if pressed(k) {
                    ed.set_tool(t);
                }
            }
        }
    }

    fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("New").clicked() {
                    self.model = Model::new();
                    self.selected = None;
                    self.fit_pending = true;
                }
                if ui.button("Open…").clicked() {
                    self.file_prompt = Some((FileAction::Open, "part.osdoc".into()));
                }
                if ui.button("Save…").clicked() {
                    self.file_prompt = Some((FileAction::Save, "part.osdoc".into()));
                }
                if ui.button("Export STL…").clicked() {
                    self.file_prompt = Some((FileAction::ExportStl, "part.stl".into()));
                }
                ui.separator();
                if ui.button("Quit").clicked() {
                    ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                }
            });
            ui.menu_button("Edit", |ui| {
                if ui
                    .add_enabled(
                        self.model.session().can_undo(),
                        egui::Button::new("Undo  (Ctrl+Z)"),
                    )
                    .clicked()
                {
                    self.undo();
                }
                if ui
                    .add_enabled(
                        self.model.session().can_redo(),
                        egui::Button::new("Redo  (Ctrl+Shift+Z)"),
                    )
                    .clicked()
                {
                    self.redo();
                }
            });
            ui.menu_button("Insert", |ui| {
                ui.menu_button("Sketch", |ui| {
                    for (name, plane) in [
                        ("On Top plane", PlaneRef::Top),
                        ("On Front plane", PlaneRef::Front),
                        ("On Right plane", PlaneRef::Right),
                    ] {
                        if ui.button(name).clicked() {
                            self.new_sketch(plane);
                            ui.close();
                        }
                    }
                });
                if ui.button("Extrude…").clicked() {
                    self.sketch = None;
                    self.pick = Some(ProfilePick::new(SweepKind::Extrude));
                }
                if ui.button("Revolve…").clicked() {
                    self.sketch = None;
                    self.pick = Some(ProfilePick::new(SweepKind::Revolve));
                }
                ui.separator();
                ui.menu_button("Templates", |ui| {
                    for t in Template::ALL {
                        if ui.button(t.label()).clicked() {
                            self.dialog = Some(FeatureDialog::new(t));
                            ui.close();
                        }
                    }
                });
            });
            ui.menu_button("View", |ui| {
                for (name, preset) in [
                    ("Isometric", ViewPreset::Isometric),
                    ("Front", ViewPreset::Front),
                    ("Top", ViewPreset::Top),
                    ("Right", ViewPreset::Right),
                ] {
                    if ui.button(name).clicked() {
                        self.camera.set_view(preset);
                        self.fit();
                    }
                }
                if ui.button("Fit  (F)").clicked() {
                    self.fit();
                }
                ui.checkbox(&mut self.show_sketches, "Show sketches");
                ui.separator();
                let ortho = matches!(self.camera.projection, Projection::Orthographic);
                if ui.selectable_label(!ortho, "Perspective").clicked() {
                    self.camera.projection = Projection::Perspective {
                        fov_y_rad: openshape_render::DEFAULT_FOV_Y_RAD,
                    };
                }
                if ui.selectable_label(ortho, "Orthographic").clicked() {
                    self.camera.projection = Projection::Orthographic;
                }
                ui.separator();
                for (name, mode) in [
                    ("Shaded", DisplayMode::Shaded),
                    ("Shaded with edges", DisplayMode::ShadedWithEdges),
                    ("Wireframe", DisplayMode::Wireframe),
                ] {
                    if ui
                        .selectable_label(self.display_mode == mode, name)
                        .clicked()
                    {
                        self.display_mode = mode;
                    }
                }
            });
            ui.separator();
            if ui
                .add_enabled(self.model.session().can_undo(), egui::Button::new("⟲ Undo"))
                .clicked()
            {
                self.undo();
            }
            if ui
                .add_enabled(self.model.session().can_redo(), egui::Button::new("⟳ Redo"))
                .clicked()
            {
                self.redo();
            }
        });
    }

    fn feature_tree(&mut self, ui: &mut egui::Ui) {
        ui.heading("Features");
        ui.separator();
        let element = self.model.studio_id();
        let features = self.model.studio().features.clone();
        let mut op: Option<ModelOp> = None;
        let mut edit: Option<FeatureId> = None;
        egui::ScrollArea::vertical().show(ui, |ui| {
            for (i, f) in features.iter().enumerate() {
                let status = self.model.status(f.id).cloned();
                let (glyph, color) = match &status {
                    Some(FeatureStatus::Ok) => ("✔", Color32::from_rgb(90, 170, 90)),
                    Some(FeatureStatus::Warning(_)) => ("⚠", Color32::from_rgb(220, 150, 30)),
                    Some(FeatureStatus::Suppressed) => ("⏸", Color32::GRAY),
                    Some(FeatureStatus::Failed(_)) => ("✖", Color32::from_rgb(220, 70, 60)),
                    None => ("?", Color32::GRAY),
                };
                let r = ui
                    .horizontal(|ui| {
                        ui.colored_label(color, glyph);
                        let text = if f.suppressed {
                            egui::RichText::new(&f.name).strikethrough().weak()
                        } else {
                            egui::RichText::new(&f.name)
                        };
                        ui.selectable_label(self.selected == Some(f.id), text)
                    })
                    .inner;
                let r = match &status {
                    Some(FeatureStatus::Failed(m)) | Some(FeatureStatus::Warning(m)) => {
                        r.on_hover_text(m)
                    }
                    _ => r.on_hover_text(f.feature.type_name()),
                };
                if r.clicked() {
                    self.selected = Some(f.id);
                }
                if r.double_clicked() {
                    edit = Some(f.id);
                }
                r.context_menu(|ui| {
                    if ui.button("Edit…").clicked() {
                        edit = Some(f.id);
                        ui.close();
                    }
                    if ui.button("Rename…").clicked() {
                        self.rename = Some((f.id, f.name.clone()));
                        ui.close();
                    }
                    let label = if f.suppressed {
                        "Unsuppress"
                    } else {
                        "Suppress"
                    };
                    if ui.button(label).clicked() {
                        op = Some(ModelOp::SetSuppressed {
                            element,
                            feature: f.id,
                            suppressed: !f.suppressed,
                        });
                        ui.close();
                    }
                    if ui
                        .add_enabled(i > 0, egui::Button::new("Move up"))
                        .clicked()
                    {
                        let to = if i >= 2 {
                            InsertAt::After(features[i - 2].id)
                        } else {
                            InsertAt::Start
                        };
                        op = Some(ModelOp::MoveFeature {
                            element,
                            feature: f.id,
                            to,
                        });
                        ui.close();
                    }
                    if ui
                        .add_enabled(i + 1 < features.len(), egui::Button::new("Move down"))
                        .clicked()
                    {
                        op = Some(ModelOp::MoveFeature {
                            element,
                            feature: f.id,
                            to: InsertAt::After(features[i + 1].id),
                        });
                        ui.close();
                    }
                    if ui.button("Delete").clicked() {
                        op = Some(ModelOp::DeleteFeature {
                            element,
                            feature: f.id,
                        });
                        ui.close();
                    }
                });
            }
            if features.is_empty() {
                ui.weak("Insert a feature from the Insert menu.");
            }
        });
        if let Some(op) = op {
            self.model.apply(op);
        }
        if let Some(id) = edit {
            self.edit_feature(id);
        }
        // Properties of the selected feature's part.
        if let Some(sel) = self.selected {
            if let Some(b) = self.model.bodies().bodies.iter().find(|b| b.feature == sel) {
                ui.separator();
                ui.heading("Properties");
                let k = self.model.kernel();
                if let Ok(m) = k.mass_properties(&b.body, 1e-3) {
                    ui.label(format!("Volume: {:.3} cm³", m.volume * 1e6));
                    ui.label(format!("Area: {:.3} cm²", m.area * 1e4));
                    ui.label(format!(
                        "Centroid: ({:.3}, {:.3}, {:.3}) mm",
                        m.centroid.x * 1e3,
                        m.centroid.y * 1e3,
                        m.centroid.z * 1e3
                    ));
                }
                let t = k.topology(&b.body);
                ui.label(format!(
                    "Faces {} · Edges {} · Vertices {}",
                    t.faces, t.edges, t.vertices
                ));
            }
        }
    }

    fn sketch_toolbar(&mut self, ui: &mut egui::Ui) {
        let Some(ed) = self.sketch.as_mut() else {
            return;
        };
        let model = &mut self.model;
        let mut finish = false;
        ui.horizontal_wrapped(|ui| {
            ui.strong(
                model
                    .feature(ed.feature)
                    .map(|f| f.name.clone())
                    .unwrap_or_default(),
            );
            ui.separator();
            for t in Tool::ALL {
                let label = match t.key() {
                    Some(k) => format!("{} ({})", t.label(), k.name()),
                    None => t.label().to_string(),
                };
                if ui.selectable_label(ed.tool == t, label).clicked() {
                    ed.set_tool(t);
                }
            }
            ui.separator();
            let mut construction = ed.construction;
            if ui
                .toggle_value(&mut construction, "Construction (Q)")
                .on_hover_text("New geometry is construction; with a selection, toggles it")
                .clicked()
            {
                ed.toggle_construction(model);
            }
            ui.separator();
            for c in ConstraintTool::ALL {
                let enabled = ed.can_constrain(model, c);
                if ui
                    .add_enabled(enabled, egui::Button::new(c.label()))
                    .on_disabled_hover_text(format!("Select {}", c.needs()))
                    .clicked()
                {
                    ed.constrain(model, c);
                }
            }
            ui.separator();
            if ui
                .button("Dimension (D)")
                .on_hover_text("Dimension the selection; double-click a dimension to change it")
                .clicked()
            {
                let at = ui.ctx().input(|i| i.pointer.hover_pos()).unwrap_or_default();
                ed.dimension(model, at);
            }
            if ui
                .add_enabled(!ed.selection.is_empty(), egui::Button::new("Delete"))
                .clicked()
            {
                ed.delete_selection(model);
            }
            ui.separator();
            if ui.button("✔ Finish sketch").clicked() {
                finish = true;
            }
        });
        ui.horizontal(|ui| {
            let d = ed.display(model);
            let color = d
                .as_ref()
                .map(sketch_paint::status_color)
                .unwrap_or(Color32::GRAY);
            ui.colored_label(color, ed.status(model));
            ui.separator();
            ui.weak(ed.tool.hint(ed.clicks.len()));
            if let Some(m) = &ed.message {
                ui.separator();
                ui.label(m);
            }
        });
        if finish {
            self.finish_sketch();
        }
    }

    fn sketch_windows(&mut self, ctx: &egui::Context) {
        let Some(ed) = self.sketch.as_mut() else {
            return;
        };
        let Some(mut prompt) = ed.prompt.clone() else {
            return;
        };
        let mut commit = false;
        let mut cancel = false;
        let mut kind = None;
        egui::Window::new("Dimension")
            .title_bar(false)
            .resizable(false)
            .fixed_pos(prompt.at + egui::vec2(16.0, 16.0))
            .show(ctx, |ui| {
                if let Some(k) = prompt.distance_kind() {
                    ui.horizontal(|ui| {
                        for (name, v) in [
                            ("Aligned", DistanceKind::Aligned),
                            ("Horizontal", DistanceKind::Horizontal),
                            ("Vertical", DistanceKind::Vertical),
                        ] {
                            if ui.selectable_label(k == v, name).clicked() && k != v {
                                kind = Some(v);
                            }
                        }
                    });
                }
                ui.horizontal(|ui| {
                    let r = ui.add(
                        egui::TextEdit::singleline(&mut prompt.text).desired_width(90.0),
                    );
                    if !r.has_focus() && !r.lost_focus() {
                        r.request_focus();
                    }
                    ui.label(if prompt.is_angle() { "°" } else { "mm" });
                    if r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                        commit = true;
                    }
                    if ui.input(|i| i.key_pressed(Key::Escape)) {
                        cancel = true;
                    }
                    if ui.button("OK").clicked() {
                        commit = true;
                    }
                });
                if let Some(e) = &prompt.error {
                    ui.colored_label(sketch_paint::CONFLICT, e);
                }
            });
        if cancel {
            ed.prompt = None;
            return;
        }
        let text = prompt.text.clone();
        if let Some(p) = ed.prompt.as_mut() {
            p.text = text;
        }
        if let Some(k) = kind {
            ed.set_prompt_kind(&self.model, k);
        } else if commit {
            ed.commit_prompt(&mut self.model);
        }
    }

    fn pick_window(&mut self, ctx: &egui::Context) {
        let Some(mut p) = self.pick.take() else {
            return;
        };
        let mut open = true;
        let mut done = false;
        egui::Window::new(p.title())
            .collapsible(false)
            .resizable(false)
            .open(&mut open)
            .show(ctx, |ui| {
                ui.weak(p.hint());
                let picked = p.picked(&self.model).len();
                ui.horizontal(|ui| {
                    ui.label(format!(
                        "{} region{} picked",
                        picked,
                        if picked == 1 { "" } else { "s" }
                    ));
                    if picked != p.regions.len() {
                        ui.colored_label(
                            sketch_paint::CONFLICT,
                            format!("({} no longer exist)", p.regions.len() - picked),
                        );
                    }
                    if ui.button("Clear").clicked() {
                        p.regions.clear();
                    }
                });
                egui::Grid::new("sweep").num_columns(2).show(ui, |ui| match p.kind {
                    SweepKind::Extrude => {
                        ui.label("Distance");
                        ui.text_edit_singleline(&mut p.distance);
                        ui.end_row();
                        ui.label("Symmetric");
                        ui.checkbox(&mut p.symmetric, "");
                        ui.end_row();
                    }
                    SweepKind::Revolve => {
                        ui.label("Axis");
                        ui.horizontal(|ui| {
                            ui.label(if p.axis.is_some() { "a sketch line" } else { "none" });
                            if ui.selectable_label(p.picking_axis, "Pick…").clicked() {
                                p.picking_axis = !p.picking_axis;
                            }
                        });
                        ui.end_row();
                        ui.label("Angle");
                        ui.text_edit_singleline(&mut p.angle);
                        ui.end_row();
                    }
                });
                if let Some(e) = &p.error {
                    ui.colored_label(sketch_paint::CONFLICT, e);
                }
                ui.horizontal(|ui| {
                    if ui.button("OK").clicked() {
                        match p.build() {
                            Ok(feature) => {
                                let element = self.model.studio_id();
                                match p.editing {
                                    Some(id) => self.model.apply(ModelOp::EditFeature {
                                        element,
                                        feature: id,
                                        definition: feature,
                                    }),
                                    None => {
                                        let name = self.model.next_name(feature.type_name());
                                        let id = self.model.add_feature(&name, feature);
                                        self.selected = Some(id);
                                    }
                                }
                                done = true;
                            }
                            Err(e) => p.error = Some(e),
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        done = true;
                    }
                });
            });
        if open && !done {
            self.pick = Some(p);
        }
    }

    fn dialogs(&mut self, ctx: &egui::Context) {
        if let Some(mut d) = self.dialog.take() {
            let mut open = true;
            let mut done = false;
            let title = if d.editing.is_some() {
                "Edit feature"
            } else {
                "New feature"
            };
            egui::Window::new(title)
                .collapsible(false)
                .resizable(false)
                .open(&mut open)
                .show(ctx, |ui| {
                    ui.strong(d.template.label());
                    egui::Grid::new("fields").num_columns(2).show(ui, |ui| {
                        for field in &mut d.fields {
                            ui.label(field.label);
                            match field.kind {
                                FieldKind::Flag => {
                                    let mut b = field.text == "true";
                                    ui.checkbox(&mut b, "");
                                    field.text = b.to_string();
                                }
                                FieldKind::Length | FieldKind::Angle => {
                                    ui.text_edit_singleline(&mut field.text);
                                }
                            }
                            ui.end_row();
                        }
                    });
                    if let Some(e) = &d.error {
                        ui.colored_label(Color32::from_rgb(220, 70, 60), e);
                    }
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            match d.build() {
                                Ok(feature) => {
                                    let element = self.model.studio_id();
                                    match d.editing {
                                        Some(id) => self.model.apply(ModelOp::EditFeature {
                                            element,
                                            feature: id,
                                            definition: feature,
                                        }),
                                        None => {
                                            let name = self.model.next_name(feature.type_name());
                                            let id = self.model.add_feature(&name, feature);
                                            self.selected = Some(id);
                                            self.fit_pending = self.model.meshes().len() == 1;
                                        }
                                    }
                                    done = true;
                                }
                                Err(e) => d.error = Some(e),
                            }
                        }
                        if ui.button("Cancel").clicked() {
                            done = true;
                        }
                    });
                });
            if open && !done {
                self.dialog = Some(d);
            }
        }
        if let Some((id, mut name)) = self.rename.take() {
            let mut keep = true;
            egui::Window::new("Rename feature")
                .collapsible(false)
                .show(ctx, |ui| {
                    ui.text_edit_singleline(&mut name);
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            let element = self.model.studio_id();
                            self.model.apply(ModelOp::RenameFeature {
                                element,
                                feature: id,
                                name: name.clone(),
                            });
                            keep = false;
                        }
                        if ui.button("Cancel").clicked() {
                            keep = false;
                        }
                    });
                });
            if keep {
                self.rename = Some((id, name));
            }
        }
        if let Some((action, mut path)) = self.file_prompt.take() {
            let mut keep = true;
            let title = match action {
                FileAction::Open => "Open document",
                FileAction::Save => "Save document",
                FileAction::ExportStl => "Export STL (mm)",
            };
            egui::Window::new(title).collapsible(false).show(ctx, |ui| {
                ui.label("Path:");
                ui.text_edit_singleline(&mut path);
                ui.horizontal(|ui| {
                    if ui.button("OK").clicked() {
                        let p = std::path::PathBuf::from(&path);
                        let r = match action {
                            FileAction::Open => Model::open(&p).map(|m| {
                                self.model = m;
                                self.selected = None;
                                self.fit_pending = true;
                            }),
                            FileAction::Save => self.model.save(&p),
                            FileAction::ExportStl => self.model.export_stl(&p),
                        };
                        self.message = match r {
                            Ok(()) => format!("{title}: {path}"),
                            Err(e) => format!("{title} failed: {e}"),
                        };
                        keep = false;
                    }
                    if ui.button("Cancel").clicked() {
                        keep = false;
                    }
                });
            });
            if keep {
                self.file_prompt = Some((action, path));
            }
        }
    }

    fn viewport_ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let avail = ui.available_size();
        let (rect, response) = ui.allocate_exact_size(avail, Sense::click_and_drag());
        self.view_rect = rect;
        let ppp = ui.ctx().pixels_per_point();
        let sketching = self.sketch.is_some();
        // In sketch mode the primary button draws and drags geometry, and
        // the secondary button orbits.
        let orbit_button = if sketching {
            PointerButton::Secondary
        } else {
            PointerButton::Primary
        };
        if response.dragged_by(orbit_button) {
            let d = response.drag_delta() * ppp;
            self.camera.orbit(d.x, d.y);
        }
        if response.dragged_by(PointerButton::Middle)
            || (!sketching && response.dragged_by(PointerButton::Secondary))
        {
            let d = response.drag_delta() * ppp;
            self.camera.pan(d.x, d.y, rect.height() * ppp);
        }
        self.sketch_input(ui, &response, rect);
        self.pick_input(&response, rect);
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                self.camera.zoom((scroll as f64 * 0.002).exp());
            }
        }
        if self.fit_pending && !self.model.meshes().is_empty() {
            self.fit();
            self.fit_pending = false;
        }
        let (Some(vp), Some(rs)) = (self.viewport.as_mut(), frame.wgpu_render_state()) else {
            ui.painter().text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "No GPU device: the 3D view is unavailable",
                egui::FontId::proportional(16.0),
                Color32::LIGHT_RED,
            );
            return;
        };
        let size = (
            ((rect.width() * ppp).round() as u32).max(1),
            ((rect.height() * ppp).round() as u32).max(1),
        );
        if vp.uploaded_generation != self.model.generation {
            vp.scene.set_meshes(&rs.device, self.model.meshes());
            vp.uploaded_generation = self.model.generation;
        }
        vp.scene.set_display_mode(self.display_mode);
        // Same apparent line weight on high-DPI screens.
        vp.scene
            .set_edge_width(openshape_render::DEFAULT_EDGE_WIDTH_PX * ppp);
        if vp.texture.is_none() || vp.size != size {
            let texture = rs.device.create_texture(&wgpu::TextureDescriptor {
                label: Some("openshape viewport"),
                size: wgpu::Extent3d {
                    width: size.0,
                    height: size.1,
                    depth_or_array_layers: 1,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: VIEW_FORMAT,
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                    | wgpu::TextureUsages::TEXTURE_BINDING
                    | wgpu::TextureUsages::COPY_SRC,
                view_formats: &[],
            });
            let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
            let mut renderer = rs.renderer.write();
            match vp.texture_id {
                Some(id) => renderer.update_egui_texture_from_wgpu_texture(
                    &rs.device,
                    &view,
                    wgpu::FilterMode::Linear,
                    id,
                ),
                None => {
                    vp.texture_id = Some(renderer.register_native_texture(
                        &rs.device,
                        &view,
                        wgpu::FilterMode::Linear,
                    ))
                }
            }
            vp.texture = Some((texture, view));
            vp.size = size;
        }
        if let Some((_, view)) = &vp.texture {
            vp.scene
                .render(&rs.device, &rs.queue, &self.camera, view, size);
        }
        if let Some(id) = vp.texture_id {
            ui.painter().image(
                id,
                rect,
                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                Color32::WHITE,
            );
        }
        self.paint_overlays(&ui.painter().with_clip_rect(rect), rect);
    }

    /// Routes viewport input to sketch mode.
    fn sketch_input(&mut self, ui: &egui::Ui, response: &egui::Response, rect: egui::Rect) {
        let Some(ed) = self.sketch.as_mut() else {
            return;
        };
        let Some(frame) = ed.frame(&self.model) else {
            return;
        };
        let map = ScreenMap::new(&self.camera, frame, rect);
        let model = &mut self.model;
        let shift = ui.input(|i| i.modifiers.shift);
        if response.drag_started_by(PointerButton::Primary) {
            let from = ui
                .input(|i| i.pointer.press_origin())
                .or(response.interact_pointer_pos());
            if let Some(from) = from {
                ed.drag_start(model, &map, from);
            }
        }
        if let Some(pos) = response.interact_pointer_pos() {
            if response.dragged_by(PointerButton::Primary) {
                ed.drag_to(model, &map, pos);
            }
            if response.drag_stopped_by(PointerButton::Primary) {
                ed.drag_end(model, &map, pos);
            }
            if response.double_clicked_by(PointerButton::Primary) {
                ed.double_click(model, &map, pos);
            } else if response.clicked_by(PointerButton::Primary) {
                ed.click(model, &map, pos, shift);
            }
        }
        if response.secondary_clicked() {
            ed.escape();
        }
        if !response.dragged() {
            ed.hover(model, &map, response.hover_pos());
        }
    }

    /// Routes viewport clicks to the profile picker.
    fn pick_input(&mut self, response: &egui::Response, rect: egui::Rect) {
        let Some(p) = self.pick.as_mut() else {
            return;
        };
        let view = ScreenMap::new(&self.camera, Frame::WORLD, rect);
        if response.clicked_by(PointerButton::Primary) {
            if let Some(pos) = response.interact_pointer_pos() {
                p.click(&self.model, &view, pos);
            }
        }
        p.hover(&self.model, &view, response.hover_pos());
    }

    fn paint_overlays(&mut self, painter: &egui::Painter, rect: egui::Rect) {
        let view = ScreenMap::new(&self.camera, Frame::WORLD, rect);
        let editing = self.sketch.as_ref().map(|e| e.feature);
        if self.show_sketches || self.pick.is_some() {
            for s in self.model.sketches() {
                if Some(s.feature) != editing {
                    let d = SketchDisplay::from_result(s);
                    sketch_paint::paint_inactive(painter, &view.with_frame(s.frame), &d);
                }
            }
        }
        if let Some(p) = self.pick.as_mut() {
            let mut regions: Vec<(FeatureId, usize, Color32)> = Vec::new();
            if let Some(f) = p.sketch {
                regions.extend(
                    p.picked(&self.model)
                        .into_iter()
                        .map(|i| (f, i, sketch_paint::REGION_FILL)),
                );
            }
            if let Some((f, i)) = p.hover {
                if !regions.iter().any(|r| r.0 == f && r.1 == i) {
                    regions.push((f, i, sketch_paint::REGION_HOVER));
                }
            }
            for (f, i, color) in regions {
                let Some(frame) = self.model.sketch_result(f).map(|s| s.frame) else {
                    continue;
                };
                let outline = p.outline(&self.model, f, i).to_vec();
                sketch_paint::paint_region(painter, &view.with_frame(frame), &outline, color);
            }
            if let (Some(f), Some(axis)) = (p.sketch, p.axis) {
                if let Some(s) = self.model.sketch_result(f) {
                    let d = SketchDisplay::from_result(s);
                    if let Some(c) = d.curve(axis) {
                        for run in view.with_frame(s.frame).polyline(&c.polyline) {
                            painter.add(egui::Shape::line(
                                run,
                                egui::Stroke::new(2.5, sketch_paint::SELECTED),
                            ));
                        }
                    }
                }
            }
        }
        if let Some(ed) = &self.sketch {
            if let (Some(d), Some(frame)) = (ed.display(&self.model), ed.frame(&self.model)) {
                sketch_paint::paint_editor(painter, &view.with_frame(frame), &d, ed);
            }
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

impl eframe::App for OpenShapeApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        // Leave sketch mode if undo or a delete removed or suppressed the sketch.
        if let Some(ed) = &self.sketch {
            if self.model.sketch_result(ed.feature).is_none() {
                self.sketch = None;
                self.message = "The sketch was removed or suppressed".into();
            }
        }
        self.shortcuts(&ctx);
        egui::Panel::top("menu").show(ui, |ui| self.menu_bar(ui));
        if self.sketch.is_some() {
            egui::Panel::top("sketch tools").show(ui, |ui| self.sketch_toolbar(ui));
        }
        egui::Panel::bottom("status").show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(format!(
                    "{} parts · {} ops in log · {}",
                    self.model.meshes().len(),
                    self.model.session().log().len(),
                    self.adapter
                ));
                if !self.message.is_empty() {
                    ui.separator();
                    ui.label(&self.message);
                }
            });
        });
        egui::Panel::left("tree")
            .resizable(true)
            .default_size(240.0)
            .show(ui, |ui| self.feature_tree(ui));
        egui::CentralPanel::no_frame().show(ui, |ui| self.viewport_ui(ui, frame));
        self.dialogs(&ctx);
        self.sketch_windows(&ctx);
        self.pick_window(&ctx);
    }
}
