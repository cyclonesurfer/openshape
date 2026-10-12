//! Floating dialogs over the graphics area (ADR 0012): feature dialogs with
//! accept and cancel in their header, sketch mode's dialog, the dimension
//! prompt, the file prompt and the command palette.

use eframe::egui;
use egui::{
    Align, Area, Id, Key, Layout, Order, Pos2, Rect, Response, RichText, Sense, Stroke, Ui, Vec2,
};
use openshape_core_model::sketch::DistanceKind;
use openshape_core_model::ModelOp;
use openshape_features::FeatureStatus;
use openshape_text::Font;
use openshape_ui_kit::palette::{PaletteEntry, PaletteOutcome};
use openshape_ui_kit::{icon_button, paint_icon, Icon, Theme};

use crate::app::{Action, FileAction, OpenShapeApp, Selection};
use crate::commands::Command;
use crate::dialogs::FieldKind;
use crate::profile_pick::SweepKind;
use crate::sketch_mode::{CurveEdit, Tool};
use crate::sketch_ops::shapes::{MAX_SIDES, MIN_SIDES};
use crate::sketch_paint;
use crate::Model;

/// Width of a feature dialog, in points.
pub const DIALOG_WIDTH: f32 = 264.0;

/// What the user did with a dialog this frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DialogAction {
    /// Nothing yet.
    None,
    /// Accept (the check mark, or Enter where a dialog allows it).
    Accept,
    /// Cancel (the cross, or Escape).
    Cancel,
}

/// The accept button of a dialog header: a check mark on the accent.
fn accept_button(ui: &mut Ui, theme: &Theme, tooltip: &str) -> Response {
    let (rect, r) = ui.allocate_exact_size(Vec2::new(30.0, 24.0), Sense::click());
    let fill = if r.hovered() {
        theme.accent_hover
    } else {
        theme.accent
    };
    ui.painter().rect_filled(rect, theme.rounding, fill);
    paint_icon(
        ui.painter(),
        Rect::from_center_size(rect.center(), Vec2::splat(18.0)),
        Icon::Check,
        theme.on_accent,
    );
    r.on_hover_text(tooltip)
}

/// A dialog anchored at `at` with an icon and title, accept and cancel in
/// its header, and `body` below. Returns the action and the dialog's rect.
pub fn floating_dialog(
    ctx: &egui::Context,
    theme: &Theme,
    id: &str,
    at: Pos2,
    icon: Icon,
    title: &str,
    body: impl FnOnce(&mut Ui),
) -> (DialogAction, Rect) {
    let mut action = DialogAction::None;
    let r = Area::new(Id::new(("floating dialog", id)))
        .order(Order::Middle)
        .fixed_pos(at)
        .show(ctx, |ui| {
            egui::Frame::new()
                .fill(theme.panel)
                .stroke(Stroke::new(1.0, theme.border))
                .corner_radius(theme.rounding + 1.0)
                .shadow(ui.style().visuals.popup_shadow)
                .inner_margin(10)
                .show(ui, |ui| {
                    ui.set_width(DIALOG_WIDTH);
                    ui.horizontal(|ui| {
                        let (ir, _) = ui.allocate_exact_size(Vec2::splat(20.0), Sense::hover());
                        paint_icon(ui.painter(), ir, icon, theme.accent);
                        ui.label(
                            RichText::new(title)
                                .strong()
                                .size(theme.heading_size)
                                .color(theme.text),
                        );
                        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                            ui.spacing_mut().item_spacing.x = 4.0;
                            if accept_button(ui, theme, "Accept").clicked() {
                                action = DialogAction::Accept;
                            }
                            let cancel = icon_button(ui, theme, Icon::Close, "Cancel", false, true);
                            if cancel.clicked() {
                                action = DialogAction::Cancel;
                            }
                        });
                    });
                    ui.add_space(4.0);
                    let (line, _) = ui
                        .allocate_exact_size(Vec2::new(ui.available_width(), 1.0), Sense::hover());
                    ui.painter().rect_filled(line, 0.0, theme.border);
                    ui.add_space(6.0);
                    body(ui);
                });
        });
    (action, r.response.rect)
}

/// A labelled text field with a unit after it.
/// Width of the sketch dialog's text fields (points).
const TEXT_FIELD_WIDTH: f32 = 150.0;

/// A drop-down of the bundled fonts.
fn font_choice(ui: &mut Ui, id: &str, font: &mut Font) {
    egui::ComboBox::from_id_salt(id)
        .selected_text(font.name())
        .show_ui(ui, |ui| {
            for f in Font::ALL {
                ui.selectable_value(font, f, f.name());
            }
        });
}

fn field(ui: &mut Ui, theme: &Theme, label: &str, text: &mut String, unit: &str) {
    ui.label(RichText::new(label).color(theme.text_muted));
    ui.horizontal(|ui| {
        ui.add(egui::TextEdit::singleline(text).desired_width(110.0));
        if !unit.is_empty() {
            ui.label(RichText::new(unit).color(theme.text_muted));
        }
    });
    ui.end_row();
}

impl OpenShapeApp {
    /// Why the open dialog's previewed feature does not regenerate.
    fn preview_failure(&self) -> Option<String> {
        match self.model.preview_status()? {
            FeatureStatus::Failed(e) => Some(e.clone()),
            _ => None,
        }
    }

    fn dialog_anchor(&self) -> Pos2 {
        self.layout.graphics.min + Vec2::new(12.0, 12.0)
    }

    pub(crate) fn floating_dialogs(&mut self, ctx: &egui::Context) {
        self.template_dialog(ctx);
        self.pick_dialog(ctx);
        self.sketch_dialog(ctx);
        self.dimension_prompt(ctx);
        self.file_prompt(ctx);
    }

    fn template_dialog(&mut self, ctx: &egui::Context) {
        let Some(mut d) = self.dialog.take() else {
            return;
        };
        let theme = self.theme.clone();
        let title = match d.editing.and_then(|id| self.model.feature(id)) {
            Some(f) => f.name.clone(),
            None => d.template.label().to_string(),
        };
        let icon = Command::Template(d.template).icon().unwrap_or(Icon::Box);
        let live = match d.build() {
            Ok(_) => self.preview_failure(),
            Err(e) => Some(e),
        };
        let (mut action, rect) = floating_dialog(
            ctx,
            &theme,
            "template",
            self.dialog_anchor(),
            icon,
            &title,
            |ui| {
                egui::Grid::new("fields")
                    .num_columns(2)
                    .spacing([10.0, 6.0])
                    .show(ui, |ui| {
                        for f in &mut d.fields {
                            match f.kind {
                                FieldKind::Flag => {
                                    ui.label(RichText::new(f.label).color(theme.text_muted));
                                    let mut b = f.text == "true";
                                    ui.checkbox(&mut b, "");
                                    f.text = b.to_string();
                                    ui.end_row();
                                }
                                FieldKind::Length | FieldKind::Angle => {
                                    field(ui, &theme, f.label, &mut f.text, "");
                                }
                            }
                        }
                    });
                if let Some(e) = live.as_ref().or(d.error.as_ref()) {
                    ui.add_space(4.0);
                    ui.colored_label(theme.danger, e);
                }
            },
        );
        self.layout.dialog = Some(rect);
        if ctx.input(|i| i.key_pressed(Key::Enter)) {
            action = DialogAction::Accept;
        }
        match action {
            DialogAction::Accept => match d.build() {
                Ok(feature) => {
                    self.model.clear_preview();
                    let element = self.model.studio_id();
                    match d.editing {
                        Some(id) => self.actions.push(Action::Apply(ModelOp::EditFeature {
                            element,
                            feature: id,
                            definition: feature,
                        })),
                        None => {
                            let name = self.model.next_name(feature.type_name());
                            let first = self
                                .model
                                .mesh_owners()
                                .iter()
                                .all(|f| *f == crate::PREVIEW_FEATURE);
                            let id = self.model.add_feature(&name, feature);
                            self.selection = Some(Selection::Feature(id));
                            self.fit_pending = first;
                        }
                    }
                }
                Err(e) => {
                    d.error = Some(e);
                    self.dialog = Some(d);
                }
            },
            DialogAction::Cancel => {}
            DialogAction::None => self.dialog = Some(d),
        }
    }

    fn pick_dialog(&mut self, ctx: &egui::Context) {
        let Some(mut p) = self.pick.take() else {
            return;
        };
        let theme = self.theme.clone();
        let title = match p.editing.and_then(|id| self.model.feature(id)) {
            Some(f) => f.name.clone(),
            None => p.title().to_string(),
        };
        let icon = match p.kind {
            SweepKind::Extrude => Icon::Extrude,
            SweepKind::Revolve => Icon::Revolve,
        };
        let picked = p.picked(&self.model).len();
        let live = match p.build() {
            Ok(_) => self.preview_failure(),
            Err(e) if p.ready() => Some(e),
            Err(_) => None,
        };
        let (action, rect) = floating_dialog(
            ctx,
            &theme,
            "pick",
            self.dialog_anchor(),
            icon,
            &title,
            |ui| {
                ui.label(
                    RichText::new(p.hint())
                        .color(theme.text_muted)
                        .size(theme.small_size),
                );
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    let (ir, _) = ui.allocate_exact_size(Vec2::splat(16.0), Sense::hover());
                    paint_icon(ui.painter(), ir, Icon::Sketch, theme.text);
                    ui.label(format!(
                        "{picked} region{} picked",
                        if picked == 1 { "" } else { "s" }
                    ));
                    if picked != p.regions.len() {
                        ui.colored_label(
                            theme.danger,
                            format!("({} missing)", p.regions.len() - picked),
                        );
                    }
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui.small_button("Clear").clicked() {
                            p.regions.clear();
                        }
                    });
                });
                ui.add_space(6.0);
                egui::Grid::new("sweep")
                    .num_columns(2)
                    .spacing([10.0, 6.0])
                    .show(ui, |ui| match p.kind {
                        SweepKind::Extrude => {
                            field(ui, &theme, "Depth", &mut p.distance, "mm");
                            ui.label(RichText::new("Symmetric").color(theme.text_muted));
                            ui.checkbox(&mut p.symmetric, "");
                            ui.end_row();
                        }
                        SweepKind::Revolve => {
                            ui.label(RichText::new("Axis").color(theme.text_muted));
                            ui.horizontal(|ui| {
                                ui.label(if p.axis.is_some() {
                                    "Sketch line"
                                } else {
                                    "None"
                                });
                                if ui.selectable_label(p.picking_axis, "Pick…").clicked() {
                                    p.picking_axis = !p.picking_axis;
                                }
                            });
                            ui.end_row();
                            field(ui, &theme, "Angle", &mut p.angle, "°");
                        }
                    });
                if let Some(e) = live.as_ref().or(p.error.as_ref()) {
                    ui.add_space(4.0);
                    ui.colored_label(theme.danger, e);
                }
            },
        );
        self.layout.dialog = Some(rect);
        let action = if ctx.input(|i| i.key_pressed(Key::Enter)) {
            DialogAction::Accept
        } else {
            action
        };
        match action {
            DialogAction::Accept => match p.build() {
                Ok(feature) => {
                    self.model.clear_preview();
                    let element = self.model.studio_id();
                    match p.editing {
                        Some(id) => self.actions.push(Action::Apply(ModelOp::EditFeature {
                            element,
                            feature: id,
                            definition: feature,
                        })),
                        None => {
                            let name = self.model.next_name(feature.type_name());
                            let first = self
                                .model
                                .mesh_owners()
                                .iter()
                                .all(|f| *f == crate::PREVIEW_FEATURE);
                            let id = self.model.add_feature(&name, feature);
                            self.selection = Some(Selection::Feature(id));
                            self.fit_pending = first;
                        }
                    }
                }
                Err(e) => {
                    p.error = Some(e);
                    self.pick = Some(p);
                }
            },
            DialogAction::Cancel => {}
            DialogAction::None => self.pick = Some(p),
        }
    }

    fn sketch_dialog(&mut self, ctx: &egui::Context) {
        if let Some(ed) = self.sketch.as_mut() {
            ed.sync_text_edit(&self.model);
        }
        let Some(ed) = self.sketch.as_ref() else {
            return;
        };
        let theme = self.theme.clone();
        let model: &Model = &self.model;
        let name = model
            .feature(ed.feature)
            .map(|f| f.name.clone())
            .unwrap_or_default();
        let plane = model
            .sketch_def(ed.feature)
            .map(|d| {
                crate::planes::StdPlane::ALL
                    .into_iter()
                    .find(|p| p.plane_ref() == d.plane)
                    .map_or("Custom plane".to_string(), |p| {
                        format!("{} plane", p.label())
                    })
            })
            .unwrap_or_default();
        let display = ed.display(model);
        let status = ed.status(model);
        let status_color = display
            .as_ref()
            .map_or(theme.text_muted, |d| sketch_paint::status_color(d, &theme));
        // Why the trim, extend or split under the pointer cannot be done.
        let hint = match &ed.curve_edit {
            Some(CurveEdit { plan: Err(why), .. }) => why.clone(),
            _ => ed.hint(),
        };
        let message = ed.message.clone();
        let tool = ed.tool;
        let mut sides = ed.polygon_sides;
        let mut radius_mm = ed.fillet_radius * 1e3;
        let mut count = ed.pattern_count;
        let mut offset_mm = ed.transform_offset * 1e3;
        let mut angle_deg = ed.transform_angle.to_degrees();
        let mut factor = ed.transform_factor;
        let mut copy = ed.transform_copy;
        let mut keep = ed.transform_keep;
        let mut apply = false;
        let mut text_content = ed.text_content.clone();
        let mut text_font = ed.text_font;
        let mut text_edit = ed.text_edit.clone();
        let mut apply_text = false;
        let (action, rect) = floating_dialog(
            ctx,
            &theme,
            "sketch",
            self.dialog_anchor(),
            Icon::Sketch,
            &name,
            |ui| {
                egui::Grid::new("sketch info")
                    .num_columns(2)
                    .spacing([10.0, 4.0])
                    .show(ui, |ui| {
                        ui.label(RichText::new("Plane").color(theme.text_muted));
                        ui.label(&plane);
                        ui.end_row();
                        ui.label(RichText::new("State").color(theme.text_muted));
                        ui.label(RichText::new(&status).color(status_color));
                        ui.end_row();
                        if tool == Tool::Polygon {
                            ui.label(RichText::new("Sides").color(theme.text_muted));
                            ui.add(egui::DragValue::new(&mut sides).range(MIN_SIDES..=MAX_SIDES));
                            ui.end_row();
                        }
                        if tool == Tool::Text {
                            ui.label(RichText::new("Text").color(theme.text_muted));
                            ui.add(
                                egui::TextEdit::singleline(&mut text_content)
                                    .desired_width(TEXT_FIELD_WIDTH),
                            );
                            ui.end_row();
                            ui.label(RichText::new("Font").color(theme.text_muted));
                            font_choice(ui, "text tool font", &mut text_font);
                            ui.end_row();
                        }
                        if let Some(e) = text_edit.as_mut() {
                            ui.label(RichText::new("Text").color(theme.text_muted));
                            let field = ui.add(
                                egui::TextEdit::singleline(&mut e.content)
                                    .desired_width(TEXT_FIELD_WIDTH),
                            );
                            if field.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                                apply_text = true;
                            }
                            ui.end_row();
                            ui.label(RichText::new("Font").color(theme.text_muted));
                            ui.horizontal(|ui| {
                                font_choice(ui, "selected text font", &mut e.font);
                                if ui.button("Apply").clicked() {
                                    apply_text = true;
                                }
                            });
                            ui.end_row();
                        }
                        if tool == Tool::SketchFillet {
                            ui.label(RichText::new("Radius").color(theme.text_muted));
                            ui.add(
                                egui::DragValue::new(&mut radius_mm)
                                    .range(0.001..=100_000.0)
                                    .speed(0.1)
                                    .suffix(" mm"),
                            );
                            ui.end_row();
                        }
                        if tool.is_pattern() {
                            ui.label(RichText::new("Instances").color(theme.text_muted));
                            ui.add(egui::DragValue::new(&mut count).range(
                                crate::sketch_ops::pattern::MIN_COUNT
                                    ..=crate::sketch_ops::pattern::MAX_COUNT,
                            ));
                            ui.end_row();
                        }
                        if tool.is_transform() {
                            // The typed amount, applied with the button on
                            // its (last) row.
                            let mut row = |name: &str, value: egui::DragValue<'_>, last: bool| {
                                ui.label(RichText::new(name).color(theme.text_muted));
                                ui.horizontal(|ui| {
                                    ui.add(value);
                                    if last && ui.button("Apply").clicked() {
                                        apply = true;
                                    }
                                });
                                ui.end_row();
                            };
                            match tool {
                                Tool::Move => {
                                    row(
                                        "X",
                                        egui::DragValue::new(&mut offset_mm.x)
                                            .speed(0.1)
                                            .suffix(" mm"),
                                        false,
                                    );
                                    row(
                                        "Y",
                                        egui::DragValue::new(&mut offset_mm.y)
                                            .speed(0.1)
                                            .suffix(" mm"),
                                        true,
                                    );
                                }
                                Tool::Rotate => row(
                                    "Angle",
                                    egui::DragValue::new(&mut angle_deg)
                                        .range(-360.0..=360.0)
                                        .speed(0.5)
                                        .suffix("°"),
                                    true,
                                ),
                                _ => row(
                                    "Factor",
                                    egui::DragValue::new(&mut factor)
                                        .range(0.001..=1000.0)
                                        .speed(0.01),
                                    true,
                                ),
                            }
                        }
                    });
                if tool.is_transform() {
                    ui.checkbox(&mut copy, "Copy (keep the original)");
                    ui.add_enabled_ui(!copy, |ui| {
                        ui.checkbox(&mut keep, "Keep constraints to the rest")
                            .on_hover_text(
                                "On: the rest of the sketch follows the selection. \
                                 Off: constraints tying them are removed.",
                            );
                    });
                }
                ui.add_space(4.0);
                ui.label(
                    RichText::new(hint)
                        .size(theme.small_size)
                        .color(theme.text_muted),
                );
                if let Some(m) = &message {
                    ui.add_space(2.0);
                    ui.label(RichText::new(m).size(theme.small_size).color(theme.text));
                }
            },
        );
        self.layout.dialog = Some(rect);
        if let Some(ed) = self.sketch.as_mut() {
            ed.text_content = text_content;
            ed.text_font = text_font;
            ed.text_edit = text_edit;
            if apply_text {
                ed.apply_text(&mut self.model);
            }
            ed.polygon_sides = sides;
            ed.pattern_count = count;
            ed.transform_copy = copy;
            ed.transform_keep = keep;
            if offset_mm != ed.transform_offset * 1e3 {
                ed.transform_offset = offset_mm * 1e-3;
            }
            if angle_deg != ed.transform_angle.to_degrees() {
                ed.transform_angle = angle_deg.to_radians();
            }
            ed.transform_factor = factor;
            if apply {
                ed.apply_typed_transform(&mut self.model);
            }
            if radius_mm != ed.fillet_radius * 1e3 {
                ed.fillet_radius = radius_mm * 1e-3;
            }
        }
        match action {
            DialogAction::Accept => self.actions.push(Action::Run(Command::FinishSketch)),
            DialogAction::Cancel => self.actions.push(Action::Run(Command::CancelSketch)),
            DialogAction::None => {}
        }
    }

    fn dimension_prompt(&mut self, ctx: &egui::Context) {
        let Some(ed) = self.sketch.as_mut() else {
            return;
        };
        let Some(mut prompt) = ed.prompt.clone() else {
            return;
        };
        let theme = self.theme.clone();
        let mut commit = false;
        let mut cancel = false;
        let mut kind = None;
        Area::new(Id::new("dimension prompt"))
            .order(Order::Foreground)
            .fixed_pos(prompt.at + Vec2::new(16.0, 16.0))
            .show(ctx, |ui| {
                egui::Frame::popup(ui.style()).show(ui, |ui| {
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
                        let r = ui
                            .add(egui::TextEdit::singleline(&mut prompt.text).desired_width(90.0));
                        if !r.has_focus() && !r.lost_focus() {
                            r.request_focus();
                        }
                        ui.label(RichText::new(prompt.unit()).color(theme.text_muted));
                        if r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                            commit = true;
                        }
                        if ui.input(|i| i.key_pressed(Key::Escape)) {
                            cancel = true;
                        }
                        if accept_button(ui, &theme, "Apply").clicked() {
                            commit = true;
                        }
                    });
                    if let Some(e) = &prompt.error {
                        ui.colored_label(theme.danger, e);
                    }
                });
            });
        if cancel {
            ed.prompt = None;
            return;
        }
        if let Some(p) = ed.prompt.as_mut() {
            p.text = prompt.text.clone();
        }
        if let Some(k) = kind {
            ed.set_prompt_kind(&self.model, k);
        } else if commit {
            ed.commit_prompt(&mut self.model);
        }
    }

    fn file_prompt(&mut self, ctx: &egui::Context) {
        let Some((action, mut path)) = self.file_prompt.take() else {
            return;
        };
        let theme = self.theme.clone();
        let title = match action {
            FileAction::Open => "Open document",
            FileAction::Save => "Save document",
            FileAction::ExportStl => "Export STL (millimetres)",
        };
        let mut keep = true;
        let mut ok = false;
        egui::Modal::new(Id::new("file prompt")).show(ctx, |ui| {
            ui.set_width(380.0);
            ui.label(RichText::new(title).strong().size(theme.heading_size));
            ui.add_space(6.0);
            ui.label(RichText::new("Path").color(theme.text_muted));
            let r = ui.add(egui::TextEdit::singleline(&mut path).desired_width(f32::INFINITY));
            if r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                ok = true;
            }
            ui.add_space(8.0);
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let label = match action {
                    FileAction::Open => "Open",
                    FileAction::Save => "Save",
                    FileAction::ExportStl => "Export",
                };
                if ui
                    .add(
                        egui::Button::new(RichText::new(label).color(theme.on_accent))
                            .fill(theme.accent),
                    )
                    .clicked()
                {
                    ok = true;
                }
                if ui.button("Cancel").clicked() || ui.input(|i| i.key_pressed(Key::Escape)) {
                    keep = false;
                }
            });
        });
        if ok {
            keep = false;
            let p = std::path::PathBuf::from(&path);
            let r = match action {
                FileAction::Open => Model::open(&p).map(|m| {
                    self.model = m;
                    self.selection = None;
                    self.hidden_parts.clear();
                    self.sketch = None;
                    self.pick = None;
                    self.dialog = None;
                    self.fit_pending = true;
                    self.doc_name = p
                        .file_stem()
                        .map_or(path.clone(), |s| s.to_string_lossy().into_owned());
                }),
                FileAction::Save => self.model.save(&p).map(|()| {
                    self.doc_name = p
                        .file_stem()
                        .map_or(path.clone(), |s| s.to_string_lossy().into_owned());
                }),
                FileAction::ExportStl => self.model.export_stl(&p),
            };
            self.notify(match r {
                Ok(()) => format!("{title}: {path}"),
                Err(e) => format!("{title} failed: {e}"),
            });
        }
        if keep {
            self.file_prompt = Some((action, path));
        }
    }

    pub(crate) fn command_palette(&mut self, ctx: &egui::Context) {
        let Some(mut palette) = self.palette.take() else {
            return;
        };
        let mode = self.mode();
        let commands: Vec<Command> = Command::all()
            .into_iter()
            .filter(|c| c.available_in(mode))
            .collect();
        let texts: Vec<(String, Option<String>)> = commands
            .iter()
            .map(|c| (c.label(), c.shortcut_text()))
            .collect();
        let entries: Vec<PaletteEntry<'_>> = commands
            .iter()
            .zip(&texts)
            .map(|(c, (label, shortcut))| PaletteEntry {
                label,
                detail: c.group(),
                icon: c.icon(),
                shortcut: shortcut.as_deref(),
                keywords: c.keywords(),
                enabled: self.enabled(*c),
            })
            .collect();
        let anchor = Pos2::new(
            self.layout
                .search
                .map_or(ctx.content_rect().center().x, |r| r.center().x),
            self.layout.search.map_or(60.0, |r| r.top()),
        );
        match palette.show(ctx, &self.theme, anchor, 440.0, &entries) {
            PaletteOutcome::Open => self.palette = Some(palette),
            PaletteOutcome::Run(i) => self.actions.push(Action::Run(commands[i])),
            PaletteOutcome::Closed => {}
        }
    }
}
