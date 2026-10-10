//! The egui desktop front end.

use eframe::egui;
use egui::{Color32, Key, KeyboardShortcut, Modifiers, PointerButton, Sense};
use openshape_core_model::{FeatureId, InsertAt, ModelOp};
use openshape_features::FeatureStatus;
use openshape_kernel::Kernel;
use openshape_render::{wgpu, Camera, DisplayMode, Projection, SceneRenderer, ViewPreset};

use crate::dialogs::{FeatureDialog, FieldKind, Template};
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
                for t in Template::ALL {
                    if ui.button(t.label()).clicked() {
                        self.dialog = Some(FeatureDialog::new(t));
                    }
                }
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
        egui::ScrollArea::vertical().show(ui, |ui| {
            for (i, f) in features.iter().enumerate() {
                let status = self.model.status(f.id).cloned();
                let (glyph, color) = match &status {
                    Some(FeatureStatus::Ok) => ("✔", Color32::from_rgb(90, 170, 90)),
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
                    Some(FeatureStatus::Failed(m)) => r.on_hover_text(m),
                    _ => r.on_hover_text(f.feature.type_name()),
                };
                if r.clicked() {
                    self.selected = Some(f.id);
                }
                if r.double_clicked() {
                    self.dialog = FeatureDialog::edit(f);
                }
                r.context_menu(|ui| {
                    if ui.button("Edit…").clicked() {
                        self.dialog = FeatureDialog::edit(f);
                        if self.dialog.is_none() {
                            self.message = "This feature has no editor yet".into();
                        }
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
        let ppp = ui.ctx().pixels_per_point();
        if response.dragged_by(PointerButton::Primary) {
            let d = response.drag_delta() * ppp;
            self.camera.orbit(d.x, d.y);
        }
        if response.dragged_by(PointerButton::Secondary)
            || response.dragged_by(PointerButton::Middle)
        {
            let d = response.drag_delta() * ppp;
            self.camera.pan(d.x, d.y, rect.height() * ppp);
        }
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
    }
}

impl eframe::App for OpenShapeApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.shortcuts(&ctx);
        egui::Panel::top("menu").show(ui, |ui| self.menu_bar(ui));
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
    }
}
