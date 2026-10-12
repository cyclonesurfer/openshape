//! The shell's bars (ADR 0012): the top bar (menu, name, undo and redo,
//! command search, display controls), the feature toolbar of the current
//! mode, and the element tabs along the bottom.

use eframe::egui;
use egui::{
    Align, Align2, Color32, CornerRadius, FontId, Layout, Pos2, Rect, Response, RichText, Sense,
    Stroke, StrokeKind, Ui, UiBuilder, Vec2,
};
use openshape_render::DisplayMode;
use openshape_ui_kit::{icon_button, paint_icon, palette::search_field, Icon, Theme};

use crate::app::{Action, OpenShapeApp, Renaming, Selection};
use crate::commands::Command;
use crate::dialogs::Template;
use crate::planes::StdPlane;
use crate::sketch_mode::{ConstraintTool, Tool};

/// Height of the top bar, in points.
pub const TOP_BAR: f32 = 40.0;
/// Height of the feature toolbar, in points.
pub const TOOLBAR: f32 = 38.0;
/// Height of the element tab bar, in points.
pub const TAB_BAR: f32 = 30.0;

/// A thin vertical divider between toolbar groups.
pub(crate) fn divider(ui: &mut Ui, theme: &Theme) {
    let (rect, _) = ui.allocate_exact_size(Vec2::new(9.0, theme.toolbar_button), Sense::hover());
    let x = rect.center().x;
    ui.painter().line_segment(
        [
            Pos2::new(x, rect.top() + 5.0),
            Pos2::new(x, rect.bottom() - 5.0),
        ],
        Stroke::new(1.0, theme.border),
    );
}

/// Where [`menu_rows`] are kept in egui's memory.
fn menu_rows_id() -> egui::Id {
    egui::Id::new("openshape menu rows")
}

/// The menu rows drawn in the last frame, by label: for tests and
/// automation, which cannot see inside egui's popups otherwise.
pub fn menu_rows(ctx: &egui::Context) -> Vec<(String, Rect)> {
    ctx.data(|d| d.get_temp::<Vec<(String, Rect)>>(menu_rows_id()))
        .unwrap_or_default()
}

/// Forgets the previous frame's menu rows; called at the start of a frame.
pub(crate) fn clear_menu_rows(ctx: &egui::Context) {
    ctx.data_mut(|d| d.insert_temp(menu_rows_id(), Vec::<(String, Rect)>::new()));
}

/// The width of a menu whose rows all fit in it.
pub(crate) const MENU_WIDTH: f32 = 210.0;

/// A menu row: icon, label and the shortcut right-aligned.
pub(crate) fn menu_row(
    ui: &mut Ui,
    theme: &Theme,
    icon: Option<Icon>,
    label: &str,
    shortcut: Option<&str>,
    enabled: bool,
) -> Response {
    // Wide enough for the icon, the label and the shortcut; every row of a
    // short menu gets the same width so the hover highlights line up.
    let text_width = |t: &str, size: f32| {
        ui.fonts_mut(|f| {
            f.layout_no_wrap(t.to_string(), FontId::proportional(size), theme.text)
                .size()
                .x
        })
    };
    let needed = 30.0
        + text_width(label, theme.body_size)
        + shortcut.map_or(0.0, |s| 24.0 + text_width(s, theme.small_size))
        + 12.0;
    let width = needed.max(MENU_WIDTH);
    let sense = if enabled {
        Sense::click()
    } else {
        Sense::hover()
    };
    let (rect, r) = ui.allocate_exact_size(Vec2::new(width, theme.interact_height), sense);
    ui.ctx().data_mut(|d| {
        d.get_temp_mut_or_default::<Vec<(String, Rect)>>(menu_rows_id())
            .push((label.to_string(), rect))
    });
    let p = ui.painter();
    if enabled && r.hovered() {
        p.rect_filled(rect, theme.rounding, theme.hover_soft);
    }
    let ink = if enabled {
        theme.text
    } else {
        theme.text_disabled
    };
    if let Some(icon) = icon {
        let ir = Rect::from_center_size(
            Pos2::new(rect.left() + 14.0, rect.center().y),
            Vec2::splat(16.0),
        );
        paint_icon(p, ir, icon, ink);
    }
    p.text(
        Pos2::new(rect.left() + 30.0, rect.center().y),
        Align2::LEFT_CENTER,
        label,
        FontId::proportional(theme.body_size),
        ink,
    );
    if let Some(s) = shortcut {
        p.text(
            Pos2::new(rect.right() - 8.0, rect.center().y),
            Align2::RIGHT_CENTER,
            s,
            FontId::proportional(theme.small_size),
            theme.text_muted,
        );
    }
    r
}

/// A filled accent button with an icon and a label (the primary action).
pub(crate) fn accent_button(ui: &mut Ui, theme: &Theme, icon: Icon, label: &str) -> Response {
    let font = FontId::proportional(theme.body_size);
    let text_w = ui
        .painter()
        .layout_no_wrap(label.to_string(), font.clone(), theme.on_accent)
        .size()
        .x;
    let size = Vec2::new(text_w + 40.0, theme.toolbar_button);
    let (rect, r) = ui.allocate_exact_size(size, Sense::click());
    let fill = if r.hovered() {
        theme.accent_hover
    } else {
        theme.accent
    };
    let p = ui.painter();
    p.rect_filled(rect, theme.rounding, fill);
    paint_icon(
        p,
        Rect::from_center_size(
            Pos2::new(rect.left() + 16.0, rect.center().y),
            Vec2::splat(18.0),
        ),
        icon,
        theme.on_accent,
    );
    p.text(
        Pos2::new(rect.left() + 30.0, rect.center().y),
        Align2::LEFT_CENTER,
        label,
        font,
        theme.on_accent,
    );
    r
}

/// The OpenShape mark: an accent tile with a white cube.
fn logo(ui: &mut Ui, theme: &Theme) {
    let (rect, _) = ui.allocate_exact_size(Vec2::splat(24.0), Sense::hover());
    ui.painter().rect_filled(rect, 5.0, theme.accent);
    paint_icon(ui.painter(), rect.shrink(3.0), Icon::Box, theme.on_accent);
}

fn bar_frame(fill: Color32) -> egui::Frame {
    egui::Frame::new()
        .fill(fill)
        .inner_margin(egui::Margin::symmetric(8, 0))
}

/// Paints a 1-point line along one edge of `rect`.
fn edge(ui: &Ui, rect: Rect, top: bool, color: Color32) {
    let y = if top {
        rect.top() + 0.5
    } else {
        rect.bottom() - 0.5
    };
    ui.painter().line_segment(
        [Pos2::new(rect.left(), y), Pos2::new(rect.right(), y)],
        Stroke::new(1.0, color),
    );
}

impl OpenShapeApp {
    /// An icon button for a command; records its rectangle and queues the
    /// command when clicked.
    pub(crate) fn command_button(&mut self, ui: &mut Ui, c: Command) -> Response {
        let icon = c.icon().unwrap_or(Icon::Info);
        let r = icon_button(
            ui,
            &self.theme,
            icon,
            &c.tooltip(),
            self.active(c),
            self.enabled(c),
        );
        self.layout.commands.push((c, r.rect));
        if r.clicked() {
            self.actions.push(Action::Run(c));
        }
        r
    }

    pub(crate) fn top_bar(&mut self, ui: &mut Ui) {
        let theme = self.theme.clone();
        egui::Panel::top("top bar")
            .exact_size(TOP_BAR)
            .resizable(false)
            .show_separator_line(false)
            .frame(bar_frame(theme.panel))
            .show(ui, |ui| {
                edge(
                    ui,
                    ui.max_rect().expand2(Vec2::new(8.0, 0.0)),
                    false,
                    theme.border,
                );
                let full = ui.max_rect();
                ui.horizontal_centered(|ui| {
                    ui.spacing_mut().item_spacing.x = 4.0;
                    let menu = icon_button(ui, &theme, Icon::Menu, "Menu", false, true);
                    egui::Popup::menu(&menu).show(|ui| {
                        for c in [
                            Command::NewDocument,
                            Command::OpenDocument,
                            Command::SaveDocument,
                            Command::ExportStl,
                        ] {
                            let s = c.shortcut_text();
                            if menu_row(ui, &theme, c.icon(), &c.label(), s.as_deref(), true)
                                .clicked()
                            {
                                self.actions.push(Action::Run(c));
                                ui.close();
                            }
                        }
                        ui.separator();
                        if menu_row(ui, &theme, Some(Icon::Close), "Quit", None, true).clicked() {
                            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                        }
                    });
                    ui.add_space(4.0);
                    logo(ui, &theme);
                    ui.label(
                        RichText::new("OpenShape")
                            .size(theme.heading_size)
                            .color(theme.text)
                            .strong(),
                    );
                    divider(ui, &theme);
                    ui.label(RichText::new(&self.doc_name).color(theme.text_muted));
                    divider(ui, &theme);
                    self.command_button(ui, Command::Undo);
                    self.command_button(ui, Command::Redo);
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        ui.spacing_mut().item_spacing.x = 2.0;
                        self.command_button(ui, Command::FitView);
                        divider(ui, &theme);
                        self.command_button(ui, Command::Orthographic);
                        self.command_button(ui, Command::Perspective);
                        divider(ui, &theme);
                        for m in [
                            DisplayMode::Wireframe,
                            DisplayMode::ShadedWithEdges,
                            DisplayMode::Shaded,
                        ] {
                            self.command_button(ui, Command::Display(m));
                        }
                    });
                });
                // The search field, centred on the bar.
                let w = 320.0f32.min(full.width() * 0.3);
                let rect = Rect::from_center_size(
                    full.center(),
                    Vec2::new(w, theme.interact_height + 4.0),
                );
                let mut child = ui.new_child(UiBuilder::new().max_rect(rect));
                let r = search_field(&mut child, &theme, "Search commands", "Ctrl+K", w);
                self.layout.search = Some(r.rect);
                if r.clicked() {
                    self.actions.push(Action::Run(Command::CommandSearch));
                }
            });
    }

    pub(crate) fn toolbar(&mut self, ui: &mut Ui) {
        let theme = self.theme.clone();
        egui::Panel::top("toolbar")
            .exact_size(TOOLBAR)
            .resizable(false)
            .show_separator_line(false)
            .frame(bar_frame(theme.panel))
            .show(ui, |ui| {
                edge(
                    ui,
                    ui.max_rect().expand2(Vec2::new(8.0, 0.0)),
                    false,
                    theme.border,
                );
                ui.horizontal_centered(|ui| {
                    ui.spacing_mut().item_spacing.x = 2.0;
                    if self.sketch.is_some() {
                        self.sketch_tools(ui);
                    } else {
                        self.feature_tools(ui);
                    }
                });
            });
    }

    fn feature_tools(&mut self, ui: &mut Ui) {
        let theme = self.theme.clone();
        let sketch = self.command_button(ui, Command::Sketch);
        let mut open = self.plane_menu;
        egui::Popup::from_response(&sketch)
            .open_bool(&mut open)
            .show(|ui| {
                ui.label(
                    RichText::new("Sketch on a plane")
                        .color(theme.text_muted)
                        .size(theme.small_size),
                );
                for p in StdPlane::ALL {
                    if menu_row(
                        ui,
                        &theme,
                        Some(Icon::Plane),
                        &format!("{} plane", p.label()),
                        None,
                        true,
                    )
                    .clicked()
                    {
                        self.actions.push(Action::Run(Command::SketchOn(p)));
                    }
                }
            });
        self.plane_menu = open && self.sketch.is_none();
        divider(ui, &theme);
        self.command_button(ui, Command::Extrude);
        self.command_button(ui, Command::Revolve);
        divider(ui, &theme);
        self.command_button(ui, Command::Template(Template::Box));
        self.command_button(ui, Command::Template(Template::Cylinder));
        let more = ui.add(
            egui::Button::new(RichText::new("Templates  ▾").color(theme.text))
                .frame(false)
                .min_size(Vec2::new(0.0, theme.toolbar_button)),
        );
        egui::Popup::menu(&more).show(|ui| {
            for t in Template::ALL {
                let c = Command::Template(t);
                if menu_row(ui, &theme, c.icon(), t.label(), None, true).clicked() {
                    self.actions.push(Action::Run(c));
                    ui.close();
                }
            }
        });
        if self.model.studio().rollback != openshape_core_model::InsertAt::End {
            divider(ui, &theme);
            self.command_button(ui, Command::RollToEnd);
        }
    }

    fn sketch_tools(&mut self, ui: &mut Ui) {
        let theme = self.theme.clone();
        self.command_button(ui, Command::Tool(Tool::Select));
        divider(ui, &theme);
        for t in Tool::ALL.into_iter().filter(|t| *t != Tool::Select) {
            self.command_button(ui, Command::Tool(t));
        }
        divider(ui, &theme);
        self.command_button(ui, Command::Construction);
        self.command_button(ui, Command::Dimension);
        divider(ui, &theme);
        for c in ConstraintTool::ALL {
            let r = self.command_button(ui, Command::Constrain(c));
            if !self.enabled(Command::Constrain(c)) {
                r.on_hover_text(format!("Select {}", c.needs()));
            }
        }
        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            let r = accent_button(ui, &theme, Icon::FinishSketch, "Finish");
            self.layout.commands.push((Command::FinishSketch, r.rect));
            if r.on_hover_text("Finish sketch").clicked() {
                self.actions.push(Action::Run(Command::FinishSketch));
            }
        });
    }

    pub(crate) fn tab_bar(&mut self, ui: &mut Ui) {
        let theme = self.theme.clone();
        egui::Panel::bottom("tabs")
            .exact_size(TAB_BAR)
            .resizable(false)
            .show_separator_line(false)
            .frame(bar_frame(theme.chrome))
            .show(ui, |ui| {
                edge(
                    ui,
                    ui.max_rect().expand2(Vec2::new(8.0, 0.0)),
                    true,
                    theme.border,
                );
                ui.horizontal_centered(|ui| {
                    ui.spacing_mut().item_spacing.x = 1.0;
                    let add = icon_button(ui, &theme, Icon::Plus, "Add a tab", false, true);
                    self.layout.add_tab = Some(add.rect);
                    egui::Popup::menu(&add).show(|ui| {
                        let c = Command::NewPartStudio;
                        if menu_row(ui, &theme, c.icon(), "Part Studio", None, true).clicked() {
                            self.actions.push(Action::Run(c));
                            ui.close();
                        }
                    });
                    ui.add_space(4.0);
                    let active = self.model.studio_id();
                    let elements: Vec<(openshape_core_model::ElementId, String)> = self
                        .model
                        .elements()
                        .iter()
                        .map(|e| (e.id, e.name.clone()))
                        .collect();
                    let can_delete = elements.len() > 1;
                    for (id, name) in elements {
                        self.tab(ui, &theme, id, &name, id == active, can_delete);
                    }
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        let parts = self.model.meshes().len();
                        ui.label(
                            RichText::new(format!(
                                "{parts} part{} · {}",
                                if parts == 1 { "" } else { "s" },
                                self.adapter
                            ))
                            .size(theme.small_size)
                            .color(theme.text_muted),
                        );
                    });
                });
            });
    }

    fn tab(
        &mut self,
        ui: &mut Ui,
        theme: &Theme,
        id: openshape_core_model::ElementId,
        name: &str,
        active: bool,
        can_delete: bool,
    ) {
        let font = FontId::proportional(theme.body_size);
        if let Some(Renaming::Element(rid, text)) = self.renaming.as_mut() {
            if *rid == id {
                let r = ui.add(egui::TextEdit::singleline(text).desired_width(140.0));
                if !r.has_focus() && !r.lost_focus() {
                    r.request_focus();
                }
                self.layout.tabs.push((id, r.rect));
                if r.lost_focus() {
                    let cancel = ui.input(|i| i.key_pressed(egui::Key::Escape));
                    self.actions.push(if cancel {
                        Action::CancelRename
                    } else {
                        Action::CommitRename
                    });
                }
                return;
            }
        }
        let text_w = ui
            .painter()
            .layout_no_wrap(name.to_string(), font.clone(), theme.text)
            .size()
            .x;
        let size = Vec2::new(text_w + 40.0, TAB_BAR - 4.0);
        let (rect, r) = ui.allocate_exact_size(size, Sense::click());
        self.layout.tabs.push((id, rect));
        let p = ui.painter();
        let top = CornerRadius {
            nw: theme.rounding as u8,
            ne: theme.rounding as u8,
            sw: 0,
            se: 0,
        };
        if active {
            p.rect_filled(rect, top, theme.panel);
            p.rect_stroke(
                rect,
                top,
                Stroke::new(1.0, theme.border),
                StrokeKind::Inside,
            );
            p.line_segment(
                [
                    rect.left_top() + Vec2::new(1.0, 1.0),
                    rect.right_top() + Vec2::new(-1.0, 1.0),
                ],
                Stroke::new(2.0, theme.accent),
            );
        } else if r.hovered() {
            p.rect_filled(rect, top, theme.hover_soft);
        }
        let ink = if active { theme.text } else { theme.text_muted };
        paint_icon(
            p,
            Rect::from_center_size(
                Pos2::new(rect.left() + 15.0, rect.center().y),
                Vec2::splat(16.0),
            ),
            openshape_ui_kit::Icon::PartStudio,
            ink,
        );
        p.text(
            Pos2::new(rect.left() + 28.0, rect.center().y),
            Align2::LEFT_CENTER,
            name,
            font,
            ink,
        );
        if r.clicked() {
            self.actions.push(Action::Activate(id));
        }
        if r.double_clicked() {
            self.actions
                .push(Action::StartRename(Renaming::Element(id, name.to_string())));
        }
        r.context_menu(|ui| {
            if menu_row(ui, theme, None, "Rename…", None, true).clicked() {
                self.actions
                    .push(Action::StartRename(Renaming::Element(id, name.to_string())));
                ui.close();
            }
            if menu_row(ui, theme, Some(Icon::Close), "Delete", None, can_delete).clicked() {
                self.actions.push(Action::DeleteElement(id));
                ui.close();
            }
        });
    }

    /// The selection, if it is a default plane.
    pub(crate) fn selected_plane(&self) -> Option<StdPlane> {
        match self.selection {
            Some(Selection::Plane(p)) => Some(p),
            _ => None,
        }
    }
}
