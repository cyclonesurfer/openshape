//! The feature list (ADR 0012): the default geometry, the features in
//! regeneration order with their status, the draggable rollback bar
//! (ADR 0013), the parts, and the selected part's properties.

use eframe::egui;
use egui::{
    Align2, Color32, CursorIcon, FontId, Pos2, Rect, Response, RichText, Sense, Stroke, Ui, Vec2,
};
use openshape_core_model::{Feature, FeatureId, InsertAt, ModelOp};
use openshape_features::FeatureStatus;
use openshape_kernel::Kernel;
use openshape_ui_kit::{paint_icon, Icon, Theme};

use crate::app::{Action, OpenShapeApp, PartProperties, Renaming, Selection};
use crate::commands::Command;
use crate::planes::StdPlane;
use crate::shell::menu_row;

/// Height of a list row, in points.
pub const ROW: f32 = 26.0;
/// Height of the rollback bar's hit area, in points.
pub const BAR: f32 = 10.0;

/// How a row is drawn.
struct RowLook<'a> {
    icon: Icon,
    label: &'a str,
    selected: bool,
    preselected: bool,
    muted: bool,
    struck: bool,
    status: Option<(Icon, Color32)>,
    /// `Some(visible)` draws an eye toggle.
    eye: Option<bool>,
}

/// What happened to a row.
struct RowResponse {
    response: Response,
    eye_clicked: bool,
}

fn row(ui: &mut Ui, theme: &Theme, look: RowLook<'_>) -> RowResponse {
    let width = ui.available_width();
    let (rect, response) = ui.allocate_exact_size(Vec2::new(width, ROW), Sense::click());
    let eye_rect = Rect::from_center_size(
        Pos2::new(rect.right() - 14.0, rect.center().y),
        Vec2::splat(20.0),
    );
    let eye = look
        .eye
        .map(|_| ui.interact(eye_rect, response.id.with("eye"), Sense::click()));
    let p = ui.painter();
    if look.selected {
        p.rect_filled(rect, theme.rounding, theme.accent_soft);
    } else if response.hovered() || look.preselected {
        p.rect_filled(rect, theme.rounding, theme.hover_soft);
    }
    let ink = if look.muted {
        theme.text_disabled
    } else {
        theme.text
    };
    let icon_ink = if look.selected && !look.muted {
        theme.accent
    } else {
        ink
    };
    paint_icon(
        p,
        Rect::from_center_size(
            Pos2::new(rect.left() + 16.0, rect.center().y),
            Vec2::splat(16.0),
        ),
        look.icon,
        icon_ink,
    );
    let font = FontId::proportional(theme.body_size);
    let galley = p.layout_no_wrap(look.label.to_string(), font, ink);
    let text_pos = Pos2::new(rect.left() + 32.0, rect.center().y - galley.size().y / 2.0);
    let text_w = galley.size().x;
    p.galley(text_pos, galley, ink);
    if look.struck {
        let y = rect.center().y;
        p.line_segment(
            [Pos2::new(text_pos.x, y), Pos2::new(text_pos.x + text_w, y)],
            Stroke::new(1.0, ink),
        );
    }
    let mut right = rect.right() - 6.0;
    if let (Some(visible), Some(eye)) = (look.eye, &eye) {
        right = eye_rect.left() - 2.0;
        if !visible || response.hovered() || eye.hovered() {
            let icon = if visible { Icon::Eye } else { Icon::EyeOff };
            let c = if eye.hovered() {
                theme.accent
            } else {
                theme.text_muted
            };
            paint_icon(p, eye_rect.shrink(2.0), icon, c);
        }
    }
    if let Some((icon, color)) = look.status {
        paint_icon(
            p,
            Rect::from_center_size(Pos2::new(right - 9.0, rect.center().y), Vec2::splat(16.0)),
            icon,
            color,
        );
    }
    RowResponse {
        eye_clicked: eye.is_some_and(|e| e.clicked()),
        response,
    }
}

fn section_header(ui: &mut Ui, theme: &Theme, title: &str, count: Option<usize>) {
    let (rect, _) = ui.allocate_exact_size(Vec2::new(ui.available_width(), 28.0), Sense::hover());
    let text = match count {
        Some(n) => format!("{title} ({n})"),
        None => title.to_string(),
    };
    ui.painter().text(
        Pos2::new(rect.left() + 6.0, rect.center().y),
        Align2::LEFT_CENTER,
        text,
        FontId::proportional(theme.heading_size),
        theme.text,
    );
}

/// The icon of a feature type.
pub fn feature_icon(f: &Feature) -> Icon {
    match f {
        Feature::Sketch(_) => Icon::Sketch,
        Feature::Box { .. } => Icon::Box,
        Feature::Extrude { .. } => Icon::Extrude,
        Feature::Revolve { .. } => Icon::Revolve,
    }
}

impl OpenShapeApp {
    pub(crate) fn feature_list(&mut self, ui: &mut Ui) {
        let theme = self.theme.clone();
        egui::Panel::left("feature list")
            .resizable(true)
            .default_size(250.0)
            .size_range(180.0..=480.0)
            .show_separator_line(true)
            .frame(theme.panel_frame())
            .show(ui, |ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        self.features_section(ui, &theme);
                        ui.add_space(10.0);
                        self.parts_section(ui, &theme);
                    });
            });
    }

    fn features_section(&mut self, ui: &mut Ui, theme: &Theme) {
        let studio = self.model.studio().clone();
        section_header(ui, theme, "Features", Some(studio.features.len()));
        // Default geometry.
        let origin = row(
            ui,
            theme,
            RowLook {
                icon: Icon::Origin,
                label: "Origin",
                selected: self.selection == Some(Selection::Origin),
                preselected: false,
                muted: false,
                struck: false,
                status: None,
                eye: None,
            },
        );
        if origin.response.clicked() {
            self.actions.push(Action::Select(Some(Selection::Origin)));
        }
        self.hover_plane = None;
        for p in StdPlane::ALL {
            let visible = self.planes_visible[p.index()];
            let r = row(
                ui,
                theme,
                RowLook {
                    icon: Icon::Plane,
                    label: p.label(),
                    selected: self.selection == Some(Selection::Plane(p)),
                    preselected: false,
                    muted: !visible,
                    struck: false,
                    status: None,
                    eye: Some(visible),
                },
            );
            self.layout.plane_rows.push((p, r.response.rect));
            if r.response.hovered() {
                self.hover_plane = Some(p);
            }
            if r.eye_clicked {
                self.actions.push(Action::TogglePlane(p));
            } else if r.response.clicked() {
                self.actions.push(Action::Select(Some(Selection::Plane(p))));
            }
            if r.response.double_clicked() {
                self.actions.push(Action::Run(Command::SketchOn(p)));
            }
            r.response.context_menu(|ui| {
                if menu_row(
                    ui,
                    theme,
                    Some(Icon::Sketch),
                    "New sketch on this plane",
                    None,
                    true,
                )
                .clicked()
                {
                    self.actions.push(Action::Run(Command::SketchOn(p)));
                    ui.close();
                }
                let (icon, label) = if visible {
                    (Icon::EyeOff, "Hide")
                } else {
                    (Icon::Eye, "Show")
                };
                if menu_row(ui, theme, Some(icon), label, None, true).clicked() {
                    self.actions.push(Action::TogglePlane(p));
                    ui.close();
                }
            });
        }
        // A hairline between the default geometry and the features.
        let (line, _) =
            ui.allocate_exact_size(Vec2::new(ui.available_width(), 7.0), Sense::hover());
        ui.painter().line_segment(
            [
                Pos2::new(line.left() + 6.0, line.center().y),
                Pos2::new(line.right() - 6.0, line.center().y),
            ],
            Stroke::new(1.0, theme.border),
        );

        let bar = studio.rollback_index();
        let mut rows: Vec<Rect> = Vec::with_capacity(studio.features.len());
        let mut bar_response: Option<Response> = None;
        for (i, f) in studio.features.iter().enumerate() {
            if i == bar {
                bar_response = Some(self.rollback_bar(ui, theme));
            }
            rows.push(self.feature_row(ui, theme, f, i, i >= bar));
        }
        if bar == studio.features.len() {
            bar_response = Some(self.rollback_bar(ui, theme));
        }
        if studio.features.is_empty() {
            ui.add_space(4.0);
            ui.label(
                RichText::new("Start with a sketch on a plane, or a template from the toolbar.")
                    .size(theme.small_size)
                    .color(theme.text_muted),
            );
        }
        // Dragging the bar: show where it would go, move it on release.
        if let Some(r) = bar_response {
            if r.hovered() || r.dragged() {
                ui.ctx().set_cursor_icon(CursorIcon::ResizeVertical);
            }
            if let Some(y) = r
                .interact_pointer_pos()
                .map(|p| p.y)
                .filter(|_| r.dragged() || r.drag_stopped())
            {
                let gap = nearest_gap(&rows, y);
                let gy = gap_y(&rows, gap, r.rect);
                if r.dragged() {
                    let x = (r.rect.left() + 6.0, r.rect.right() - 6.0);
                    ui.painter().line_segment(
                        [Pos2::new(x.0, gy), Pos2::new(x.1, gy)],
                        Stroke::new(3.0, theme.preselect),
                    );
                }
                if r.drag_stopped() && gap != bar {
                    let at = if gap >= studio.features.len() {
                        InsertAt::End
                    } else if gap == 0 {
                        InsertAt::Start
                    } else {
                        InsertAt::After(studio.features[gap - 1].id)
                    };
                    self.actions.push(Action::SetRollback(at));
                }
            }
        }
    }

    fn rollback_bar(&mut self, ui: &mut Ui, theme: &Theme) -> Response {
        let (rect, r) = ui.allocate_exact_size(Vec2::new(ui.available_width(), BAR), Sense::drag());
        self.layout.rollback_bar = Some(rect);
        let color = if r.hovered() || r.dragged() {
            theme.accent_hover
        } else {
            theme.accent
        };
        let y = rect.center().y;
        let p = ui.painter();
        p.line_segment(
            [
                Pos2::new(rect.left() + 6.0, y),
                Pos2::new(rect.right() - 6.0, y),
            ],
            Stroke::new(2.0, color),
        );
        // A grip at the left end.
        let grip = Rect::from_center_size(Pos2::new(rect.left() + 16.0, y), Vec2::new(18.0, 8.0));
        p.rect_filled(grip, 4.0, color);
        for k in [-4.0, 0.0, 4.0] {
            p.circle_filled(Pos2::new(grip.center().x + k, y), 1.0, theme.on_accent);
        }
        r.on_hover_text("Rollback bar: drag to roll the history back or forward")
    }

    fn feature_row(
        &mut self,
        ui: &mut Ui,
        theme: &Theme,
        f: &openshape_core_model::FeatureEntry,
        index: usize,
        rolled_back: bool,
    ) -> Rect {
        let status = self.model.status(f.id).cloned();
        if let Some(Renaming::Feature(id, text)) = self.renaming.as_mut() {
            if *id == f.id {
                let (rect, _) =
                    ui.allocate_exact_size(Vec2::new(ui.available_width(), ROW), Sense::hover());
                let edit_rect = Rect::from_min_max(
                    Pos2::new(rect.left() + 28.0, rect.top() + 2.0),
                    Pos2::new(rect.right() - 6.0, rect.bottom() - 2.0),
                );
                paint_icon(
                    ui.painter(),
                    Rect::from_center_size(
                        Pos2::new(rect.left() + 16.0, rect.center().y),
                        Vec2::splat(16.0),
                    ),
                    feature_icon(&f.feature),
                    theme.accent,
                );
                let r = ui.put(edit_rect, egui::TextEdit::singleline(text));
                if !r.has_focus() && !r.lost_focus() {
                    r.request_focus();
                }
                if r.lost_focus() {
                    let cancel = ui.input(|i| i.key_pressed(egui::Key::Escape));
                    self.actions.push(if cancel {
                        Action::CancelRename
                    } else {
                        Action::CommitRename
                    });
                }
                self.layout.feature_rows.push((f.id, rect));
                return rect;
            }
        }
        let status_icon = match &status {
            Some(FeatureStatus::Failed(_)) => Some((Icon::Error, theme.danger)),
            Some(FeatureStatus::Warning(_)) => Some((Icon::Warning, theme.warning)),
            _ => None,
        };
        let selected =
            matches!(self.selection, Some(Selection::Feature(s) | Selection::Part(s)) if s == f.id);
        let r = row(
            ui,
            theme,
            RowLook {
                icon: feature_icon(&f.feature),
                label: &f.name,
                selected,
                preselected: false,
                // Also greyed: features after the one being edited.
                muted: rolled_back
                    || f.suppressed
                    || matches!(status, Some(FeatureStatus::RolledBack)),
                struck: f.suppressed,
                status: status_icon,
                eye: None,
            },
        )
        .response;
        self.layout.feature_rows.push((f.id, r.rect));
        let r = match &status {
            Some(FeatureStatus::Failed(m)) | Some(FeatureStatus::Warning(m)) => {
                r.on_hover_text(format!("{}: {m}", f.feature.type_name()))
            }
            Some(FeatureStatus::RolledBack) => {
                r.on_hover_text(format!("{} (rolled back)", f.feature.type_name()))
            }
            Some(FeatureStatus::Suppressed) => {
                r.on_hover_text(format!("{} (suppressed)", f.feature.type_name()))
            }
            _ => r.on_hover_text(f.feature.type_name()),
        };
        if r.clicked() {
            self.actions
                .push(Action::Select(Some(Selection::Feature(f.id))));
        }
        if r.double_clicked() {
            self.actions.push(Action::EditFeature(f.id));
        }
        let element = self.model.studio_id();
        let up = self.move_feature_ops(f.id, true);
        let down = self.move_feature_ops(f.id, false);
        let bar = self.model.studio().rollback_index();
        r.context_menu(|ui| {
            if menu_row(ui, theme, Some(Icon::Sketch), "Edit…", None, !rolled_back).clicked() {
                self.actions.push(Action::EditFeature(f.id));
                ui.close();
            }
            if menu_row(ui, theme, None, "Rename…", None, true).clicked() {
                self.actions
                    .push(Action::StartRename(Renaming::Feature(f.id, f.name.clone())));
                ui.close();
            }
            let label = if f.suppressed {
                "Unsuppress"
            } else {
                "Suppress"
            };
            if menu_row(ui, theme, None, label, None, true).clicked() {
                self.actions.push(Action::Apply(ModelOp::SetSuppressed {
                    element,
                    feature: f.id,
                    suppressed: !f.suppressed,
                }));
                ui.close();
            }
            ui.separator();
            if menu_row(
                ui,
                theme,
                Some(Icon::Rollback),
                "Roll history to here",
                None,
                index + 1 != bar,
            )
            .clicked()
            {
                self.actions
                    .push(Action::SetRollback(InsertAt::After(f.id)));
                ui.close();
            }
            if menu_row(ui, theme, None, "Move up", None, up.is_some()).clicked() {
                if let Some(op) = up.clone() {
                    self.actions.push(Action::Apply(op));
                }
                ui.close();
            }
            if menu_row(ui, theme, None, "Move down", None, down.is_some()).clicked() {
                if let Some(op) = down.clone() {
                    self.actions.push(Action::Apply(op));
                }
                ui.close();
            }
            ui.separator();
            if menu_row(ui, theme, Some(Icon::Close), "Delete", Some("Delete"), true).clicked() {
                self.actions.push(Action::Apply(ModelOp::DeleteFeature {
                    element,
                    feature: f.id,
                }));
                ui.close();
            }
        });
        r.rect
    }

    fn parts_section(&mut self, ui: &mut Ui, theme: &Theme) {
        let owners: Vec<FeatureId> = self.model.mesh_owners().to_vec();
        section_header(ui, theme, "Parts", Some(owners.len()));
        for (i, f) in owners.iter().enumerate() {
            let visible = !self.hidden_parts.contains(f);
            let name = format!("Part {}", i + 1);
            let r = row(
                ui,
                theme,
                RowLook {
                    icon: Icon::Part,
                    label: &name,
                    selected: self.selection == Some(Selection::Part(*f)),
                    preselected: false,
                    muted: !visible,
                    struck: false,
                    status: None,
                    eye: Some(visible),
                },
            );
            self.layout.part_rows.push((*f, r.response.rect));
            if r.eye_clicked {
                self.actions.push(Action::TogglePart(*f));
            } else if r.response.clicked() {
                self.actions.push(Action::Select(Some(Selection::Part(*f))));
            }
            let feature_name = self
                .model
                .feature(*f)
                .map(|e| e.name.clone())
                .unwrap_or_default();
            let r = r.response.on_hover_text(format!("Made by {feature_name}"));
            r.context_menu(|ui| {
                let (icon, label) = if visible {
                    (Icon::EyeOff, "Hide")
                } else {
                    (Icon::Eye, "Show")
                };
                if menu_row(ui, theme, Some(icon), label, None, true).clicked() {
                    self.actions.push(Action::TogglePart(*f));
                    ui.close();
                }
                if menu_row(
                    ui,
                    theme,
                    Some(Icon::Eye),
                    "Show all parts",
                    None,
                    !self.hidden_parts.is_empty(),
                )
                .clicked()
                {
                    self.actions.push(Action::ShowAllParts);
                    ui.close();
                }
                if menu_row(
                    ui,
                    theme,
                    Some(Icon::ZoomFit),
                    "Zoom to part",
                    None,
                    visible,
                )
                .clicked()
                {
                    self.actions.push(Action::ZoomTo(*f));
                    ui.close();
                }
                if menu_row(
                    ui,
                    theme,
                    None,
                    &format!("Edit {feature_name}…"),
                    None,
                    true,
                )
                .clicked()
                {
                    self.actions.push(Action::EditFeature(*f));
                    ui.close();
                }
            });
        }
        if owners.is_empty() {
            ui.add_space(2.0);
            ui.label(
                RichText::new("No parts yet")
                    .size(theme.small_size)
                    .color(theme.text_muted),
            );
        }
        if let Some(Selection::Part(f)) = self.selection {
            let props = self.part_properties(f);
            ui.add_space(10.0);
            egui::Frame::new()
                .fill(theme.panel_header)
                .stroke(Stroke::new(1.0, theme.border))
                .corner_radius(theme.rounding)
                .inner_margin(8)
                .show(ui, |ui| {
                    ui.spacing_mut().item_spacing.y = 3.0;
                    ui.label(RichText::new("Properties").strong().color(theme.text));
                    egui::Grid::new("part properties")
                        .num_columns(2)
                        .spacing([12.0, 3.0])
                        .show(ui, |ui| {
                            for (k, v) in &props {
                                ui.label(
                                    RichText::new(k)
                                        .color(theme.text_muted)
                                        .size(theme.small_size),
                                );
                                ui.label(RichText::new(v).size(theme.small_size));
                                ui.end_row();
                            }
                        });
                });
        }
    }

    /// Volume, area, centroid and topology of a part (cached until the next
    /// regeneration).
    fn part_properties(&mut self, f: FeatureId) -> Vec<(String, String)> {
        let generation = self.model.generation;
        if let Some(p) = &self.properties {
            if p.generation == generation && p.feature == f {
                return p.lines.clone();
            }
        }
        let mut lines = Vec::new();
        if let Some(b) = self.model.bodies().bodies.iter().find(|b| b.feature == f) {
            let k = self.model.kernel();
            match k.mass_properties(&b.body, 1e-3) {
                Ok(m) => {
                    lines.push(("Volume".into(), format!("{:.3} cm³", m.volume * 1e6)));
                    lines.push(("Surface area".into(), format!("{:.3} cm²", m.area * 1e4)));
                    lines.push((
                        "Centroid".into(),
                        format!(
                            "{}, {}, {} mm",
                            fixed(m.centroid.x * 1e3, 2),
                            fixed(m.centroid.y * 1e3, 2),
                            fixed(m.centroid.z * 1e3, 2)
                        ),
                    ));
                }
                Err(e) => lines.push(("Mass properties".into(), e.to_string())),
            }
            let t = k.topology(&b.body);
            lines.push((
                "Topology".into(),
                format!(
                    "{} faces, {} edges, {} vertices",
                    t.faces, t.edges, t.vertices
                ),
            ));
        }
        self.properties = Some(PartProperties {
            generation,
            feature: f,
            lines: lines.clone(),
        });
        lines
    }
}

/// `v` with `decimals` decimals, without the sign of a value that rounds to
/// zero ("0.00", not "-0.00").
pub fn fixed(v: f64, decimals: usize) -> String {
    let s = format!("{v:.decimals$}");
    match s.strip_prefix('-') {
        Some(rest) if rest.chars().all(|c| c == '0' || c == '.') => rest.to_string(),
        _ => s,
    }
}

/// The gap between rows nearest to height `y`: gap `i` is above row `i`,
/// gap `rows.len()` below the last row.
pub fn nearest_gap(rows: &[Rect], y: f32) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let gaps = (0..=rows.len()).map(|g| (g, gap_y(rows, g, Rect::NOTHING)));
    gaps.min_by(|a, b| (a.1 - y).abs().total_cmp(&(b.1 - y).abs()))
        .map_or(0, |(g, _)| g)
}

fn gap_y(rows: &[Rect], gap: usize, fallback: Rect) -> f32 {
    if gap < rows.len() {
        rows[gap].top()
    } else if let Some(last) = rows.last() {
        last.bottom()
    } else {
        fallback.center().y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_that_round_to_zero_have_no_sign() {
        assert_eq!(fixed(-0.001, 2), "0.00");
        assert_eq!(fixed(-0.0, 2), "0.00");
        assert_eq!(fixed(-0.005_1, 2), "-0.01");
        assert_eq!(fixed(12.345_6, 2), "12.35");
        assert_eq!(fixed(-3.0, 0), "-3");
        assert_eq!(fixed(-0.4, 0), "0");
    }

    #[test]
    fn the_nearest_gap_follows_the_pointer() {
        let rows: Vec<Rect> = (0..3)
            .map(|i| {
                Rect::from_min_size(
                    Pos2::new(0.0, 100.0 + 30.0 * i as f32),
                    Vec2::new(200.0, 26.0),
                )
            })
            .collect();
        assert_eq!(nearest_gap(&rows, 0.0), 0);
        assert_eq!(nearest_gap(&rows, 128.0), 1);
        assert_eq!(nearest_gap(&rows, 162.0), 2);
        assert_eq!(nearest_gap(&rows, 500.0), 3);
        assert_eq!(nearest_gap(&[], 50.0), 0);
    }
}
