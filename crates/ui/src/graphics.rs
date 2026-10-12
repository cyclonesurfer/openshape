//! The graphics area (ADR 0012): the 3D view rendered offscreen and shown
//! through egui, the default planes, part and plane picking, sketches and
//! sketch-mode feedback, the view cube, the graphics-area context menu and
//! messages.

use eframe::egui;
use egui::{Align2, Color32, FontId, PointerButton, Pos2, Rect, Sense, Stroke, Ui, Vec2};
use openshape_core_model::{DVec2, FeatureId};
use openshape_kernel::Frame;
use openshape_render::{wgpu, DisplayMode, MeshData, SceneRenderer};
use openshape_ui_kit::{view_cube, Theme};

use crate::app::{Action, OpenShapeApp, Selection};
use crate::commands::Command;
use crate::picking::pick_mesh;
use crate::planes::{self, StdPlane};
use crate::shell::menu_row;
use crate::sketch_mode::{ScreenMap, SketchDisplay};
use crate::sketch_paint;
use crate::viewport::{scene_format, ViewportTexture};

/// Side of the view cube's area, in points.
pub const CUBE: f32 = 116.0;
/// How long a message stays, in seconds.
const MESSAGE_SECONDS: f64 = 6.0;

/// What the uploaded scene shows; it is re-uploaded when this changes.
#[derive(Clone, Debug, PartialEq)]
struct SceneKey {
    generation: u64,
    selected: Option<FeatureId>,
    hidden: Vec<FeatureId>,
    planes: [bool; 3],
    active_plane: Option<StdPlane>,
    /// The sketches drawn in 3D: shown, and the one being edited (drawn
    /// by the editor instead).
    sketches: Option<Option<FeatureId>>,
}

/// Offscreen target the 3D view renders into; shown by egui as a texture.
pub struct Viewport {
    scene: SceneRenderer,
    format: wgpu::TextureFormat,
    texture: Option<ViewportTexture>,
    texture_id: Option<egui::TextureId>,
    uploaded: Option<SceneKey>,
}

/// sRGB float triple of a colour.
fn rgb(c: Color32) -> [f32; 3] {
    [c.r(), c.g(), c.b()].map(|v| f32::from(v) / 255.0)
}

/// sRGB colour with linear alpha, as the renderer takes it, from an egui
/// (premultiplied) colour.
fn rgba(c: Color32) -> [f32; 4] {
    let [r, g, b, a] = c.to_srgba_unmultiplied();
    [r, g, b, a].map(|v| f32::from(v) / 255.0)
}

/// How far a selected part's colour moves towards the selection colour.
const SELECTION_TINT: f32 = 0.3;

/// Tints a part colour towards the selection colour.
fn tint(color: [f32; 4], toward: Color32, amount: f32) -> [f32; 4] {
    let t = rgb(toward);
    [
        color[0] + (t[0] - color[0]) * amount,
        color[1] + (t[1] - color[1]) * amount,
        color[2] + (t[2] - color[2]) * amount,
        color[3],
    ]
}

impl Viewport {
    /// GPU resources for the 3D view on this device.
    pub fn new(rs: &egui_wgpu::RenderState, theme: &Theme) -> Viewport {
        let format = scene_format(&rs.adapter);
        let mut scene = SceneRenderer::new(&rs.device, format);
        scene.set_background(rgb(theme.viewport_top), rgb(theme.viewport_bottom));
        Viewport {
            scene,
            format,
            texture: None,
            texture_id: None,
            uploaded: None,
        }
    }
}

impl OpenShapeApp {
    /// The default planes' half size for the current parts and sketches.
    pub(crate) fn plane_half(&self) -> f64 {
        planes::half_size(self.model.bounds(), self.sketch_bounds())
    }

    /// The plane to highlight in the view.
    fn active_plane(&self) -> Option<StdPlane> {
        self.hover_plane.or(self.selected_plane())
    }

    /// The meshes to draw: visible parts (the selected one tinted, with
    /// selection-coloured edges), the sketches not being edited (their
    /// non-construction curves as lines, hidden by the parts in front of
    /// them), then the visible default planes.
    fn scene_meshes(&self) -> Vec<MeshData> {
        let selected = match self.selection {
            Some(Selection::Part(f) | Selection::Feature(f)) => Some(f),
            _ => None,
        };
        // An open dialog's feature shows as selected.
        let previewed = self.model.preview().map(|p| p.feature());
        let mut out: Vec<MeshData> = self
            .model
            .meshes()
            .iter()
            .zip(self.model.mesh_owners())
            .filter(|(_, f)| !self.hidden_parts.contains(f))
            .map(|(m, f)| {
                if Some(*f) == selected || Some(*f) == previewed {
                    MeshData {
                        color: tint(m.color, self.theme.selection, SELECTION_TINT),
                        edge_color: Some(rgba(self.theme.selection)),
                        ..m.clone()
                    }
                } else {
                    m.clone()
                }
            })
            .collect();
        if let Some(editing) = self.sketches_in_3d() {
            for s in self.model.sketches() {
                if Some(s.feature) == editing {
                    continue;
                }
                let f = s.frame;
                let at = |p: DVec2| (f.origin + f.x * p.x + f.y * p.y).as_vec3().to_array();
                let mut lines = Vec::new();
                for c in SketchDisplay::from_result(s).curves {
                    if !c.construction {
                        for w in c.polyline.windows(2) {
                            lines.extend([at(w[0]), at(w[1])]);
                        }
                    }
                }
                out.push(MeshData {
                    edge_lines: lines,
                    edge_color: Some(rgba(self.theme.sketch.inactive)),
                    ..MeshData::default()
                });
            }
        }
        if self.sketch.is_none() {
            let h = self.plane_half();
            let k = &self.theme.sketch;
            for p in StdPlane::ALL {
                if self.planes_visible[p.index()] {
                    let (fill, edge) = if self.active_plane() == Some(p) {
                        (k.plane_active, k.plane_active_edge)
                    } else {
                        (k.plane, k.plane_edge)
                    };
                    out.push(planes::mesh(p, h, rgba(fill), rgba(edge)));
                }
            }
        }
        out
    }

    /// `Some(the sketch being edited)` when sketches are drawn in 3D.
    fn sketches_in_3d(&self) -> Option<Option<FeatureId>> {
        self.show_sketches
            .then(|| self.sketch.as_ref().map(|e| e.feature))
    }

    fn scene_key(&self) -> SceneKey {
        SceneKey {
            generation: self.model.generation,
            selected: match self.selection {
                Some(Selection::Part(f) | Selection::Feature(f)) => Some(f),
                _ => None,
            },
            hidden: self.hidden_parts.iter().copied().collect(),
            planes: if self.sketch.is_none() {
                self.planes_visible
            } else {
                [false; 3]
            },
            active_plane: self.active_plane(),
            sketches: self.sketches_in_3d(),
        }
    }

    /// The part or plane under a screen point of the graphics area.
    fn pick_at(&self, rect: Rect, pos: Pos2) -> Option<Selection> {
        let (o, d) = self.camera.pixel_ray(
            pos.x - rect.min.x,
            pos.y - rect.min.y,
            rect.width(),
            rect.height(),
        );
        let visible: Vec<(FeatureId, &MeshData)> = self
            .model
            .mesh_owners()
            .iter()
            .copied()
            .zip(self.model.meshes())
            .filter(|(f, _)| !self.hidden_parts.contains(f) && *f != crate::PREVIEW_FEATURE)
            .collect();
        let part = pick_mesh(visible.iter().map(|(_, m)| *m), o, d).map(|(i, t)| (visible[i].0, t));
        let h = self.plane_half();
        let plane = StdPlane::ALL
            .into_iter()
            .filter(|p| self.planes_visible[p.index()])
            .filter_map(|p| planes::hit(p, h, o, d).map(|t| (p, t)))
            .min_by(|a, b| a.1.total_cmp(&b.1));
        match (part, plane) {
            (Some((f, _)), _) => Some(Selection::Part(f)),
            (None, Some((p, _))) => Some(Selection::Plane(p)),
            (None, None) => None,
        }
    }

    pub(crate) fn graphics(&mut self, ui: &mut Ui, rs: Option<&egui_wgpu::RenderState>) {
        egui::CentralPanel::no_frame().show(ui, |ui| self.graphics_area(ui, rs));
    }

    fn graphics_area(&mut self, ui: &mut Ui, rs: Option<&egui_wgpu::RenderState>) {
        let avail = ui.available_size();
        let (rect, response) = ui.allocate_exact_size(avail, Sense::click_and_drag());
        self.layout.graphics = rect;
        if let Some(p) = response.hover_pos() {
            self.pointer = Some(p);
        }
        let ppp = ui.ctx().pixels_per_point();
        let sketching = self.sketch.is_some();
        // ADR 0012: the secondary button orbits and the middle button (or
        // Ctrl with the secondary) pans, in every mode. In sketch mode the
        // primary button draws and drags geometry; in a Part Studio it also
        // orbits until parts can be box-selected.
        let ctrl = ui.input(|i| i.modifiers.command);
        let secondary = response.dragged_by(PointerButton::Secondary);
        if (secondary && !ctrl) || (!sketching && response.dragged_by(PointerButton::Primary)) {
            let d = response.drag_delta() * ppp;
            self.camera.orbit(d.x, d.y);
        }
        if response.dragged_by(PointerButton::Middle) || (secondary && ctrl) {
            let d = response.drag_delta() * ppp;
            self.camera.pan(d.x, d.y, rect.height() * ppp);
        }
        self.sketch_input(ui, &response, rect);
        let picking_profile = self.pick.is_some();
        self.pick_input(&response, rect);
        if !sketching && !picking_profile {
            if response.clicked_by(PointerButton::Primary) {
                if let Some(pos) = response.interact_pointer_pos() {
                    self.actions.push(Action::Select(self.pick_at(rect, pos)));
                }
            }
            if response.double_clicked_by(PointerButton::Primary) {
                if let Some(Selection::Part(f)) = response
                    .interact_pointer_pos()
                    .and_then(|p| self.pick_at(rect, p))
                {
                    self.actions.push(Action::EditFeature(f));
                }
            }
            let theme = self.theme.clone();
            response.context_menu(|ui| self.view_menu(ui, &theme));
        }
        if self.sketch.as_ref().is_some_and(|ed| ed.menu_ready()) {
            let theme = self.theme.clone();
            response.context_menu(|ui| self.sketch_menu(ui, &theme));
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                self.camera.zoom((scroll as f64 * 0.002).exp());
            }
        }
        if self.fit_pending {
            self.fit();
            self.fit_pending = false;
        }
        self.render_scene(ui, rs, rect, ppp);
        let painter = ui.painter().with_clip_rect(rect);
        self.paint_overlays(&painter, rect);
        // The view cube, top right.
        let cube = Rect::from_min_size(
            Pos2::new(rect.right() - CUBE - 8.0, rect.top() + 8.0),
            Vec2::splat(CUBE),
        );
        self.layout.view_cube = Some(cube);
        if let Some(t) = view_cube::show(ui, cube, &self.camera, &self.theme) {
            self.actions.push(Action::SetView(t));
        }
        self.paint_message(&painter, rect);
    }

    fn view_menu(&mut self, ui: &mut Ui, theme: &Theme) {
        for c in [Command::FitView, Command::View(view_cube::CubeTarget::HOME)] {
            if menu_row(
                ui,
                theme,
                c.icon(),
                &c.label(),
                c.shortcut_text().as_deref(),
                true,
            )
            .clicked()
            {
                self.actions.push(Action::Run(c));
                ui.close();
            }
        }
        ui.separator();
        for m in [
            DisplayMode::Shaded,
            DisplayMode::ShadedWithEdges,
            DisplayMode::Wireframe,
        ] {
            let c = Command::Display(m);
            if menu_row(ui, theme, c.icon(), &c.label(), None, !self.active(c)).clicked() {
                self.actions.push(Action::Run(c));
                ui.close();
            }
        }
        ui.separator();
        let sketches = if self.show_sketches {
            "Hide sketches"
        } else {
            "Show sketches"
        };
        if menu_row(
            ui,
            theme,
            Some(openshape_ui_kit::Icon::Sketch),
            sketches,
            None,
            true,
        )
        .clicked()
        {
            self.actions.push(Action::Run(Command::ToggleSketches));
            ui.close();
        }
        let planes = if self.planes_visible.iter().any(|v| *v) {
            "Hide planes"
        } else {
            "Show planes"
        };
        if menu_row(
            ui,
            theme,
            Some(openshape_ui_kit::Icon::Plane),
            planes,
            None,
            true,
        )
        .clicked()
        {
            self.actions.push(Action::Run(Command::TogglePlanes));
            ui.close();
        }
        if !self.hidden_parts.is_empty()
            && menu_row(
                ui,
                theme,
                Some(openshape_ui_kit::Icon::Eye),
                "Show all parts",
                None,
                true,
            )
            .clicked()
        {
            self.actions.push(Action::ShowAllParts);
            ui.close();
        }
    }

    /// The context menu of sketch mode: commands on the selection, then
    /// fitting the view.
    fn sketch_menu(&mut self, ui: &mut Ui, theme: &Theme) {
        let comb = if self.active(Command::CurvatureComb) {
            "Hide curvature comb".to_string()
        } else {
            "Show curvature comb".to_string()
        };
        let rows = [
            (Command::CurvatureComb, comb),
            (Command::Construction, Command::Construction.label()),
            (Command::DeleteSelection, Command::DeleteSelection.label()),
        ];
        for (c, label) in rows {
            let shortcut = c.shortcut_text();
            if menu_row(
                ui,
                theme,
                c.icon(),
                &label,
                shortcut.as_deref(),
                self.enabled(c),
            )
            .clicked()
            {
                self.actions.push(Action::Run(c));
                ui.close();
            }
        }
        ui.separator();
        let c = Command::FitView;
        if menu_row(
            ui,
            theme,
            c.icon(),
            &c.label(),
            c.shortcut_text().as_deref(),
            true,
        )
        .clicked()
        {
            self.actions.push(Action::Run(c));
            ui.close();
        }
    }

    fn render_scene(&mut self, ui: &Ui, rs: Option<&egui_wgpu::RenderState>, rect: Rect, ppp: f32) {
        let key = self.scene_key();
        let meshes = (self.viewport.as_ref().and_then(|v| v.uploaded.as_ref()) != Some(&key))
            .then(|| self.scene_meshes());
        let (Some(vp), Some(rs)) = (self.viewport.as_mut(), rs) else {
            ui.painter().text(
                rect.center(),
                Align2::CENTER_CENTER,
                "No GPU device: the 3D view is unavailable",
                FontId::proportional(16.0),
                self.theme.danger,
            );
            return;
        };
        let size = (
            ((rect.width() * ppp).round() as u32).max(1),
            ((rect.height() * ppp).round() as u32).max(1),
        );
        if let Some(meshes) = meshes {
            vp.scene.set_meshes(&rs.device, &meshes);
            vp.uploaded = Some(key);
        }
        vp.scene.set_display_mode(self.display_mode);
        // Same apparent line weight on high-DPI screens.
        vp.scene
            .set_edge_width(openshape_render::DEFAULT_EDGE_WIDTH_PX * ppp);
        if vp.texture.as_ref().map(|t| t.size) != Some(size) {
            let t = ViewportTexture::new(&rs.device, vp.format, size);
            let mut renderer = rs.renderer.write();
            match vp.texture_id {
                Some(id) => renderer.update_egui_texture_from_wgpu_texture(
                    &rs.device,
                    &t.egui_view,
                    wgpu::FilterMode::Linear,
                    id,
                ),
                None => {
                    vp.texture_id = Some(renderer.register_native_texture(
                        &rs.device,
                        &t.egui_view,
                        wgpu::FilterMode::Linear,
                    ))
                }
            }
            vp.texture = Some(t);
        }
        if let Some(t) = &vp.texture {
            vp.scene
                .render(&rs.device, &rs.queue, &self.camera, &t.render_view, size);
        }
        if let Some(id) = vp.texture_id {
            ui.painter().image(
                id,
                rect,
                Rect::from_min_max(Pos2::new(0.0, 0.0), Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
        }
    }

    /// Routes viewport input to sketch mode.
    fn sketch_input(&mut self, ui: &Ui, response: &egui::Response, rect: Rect) {
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
            // With nothing in progress the click opens the sketch's context
            // menu, for what is under the pointer; otherwise it cancels.
            if ed.menu_ready() {
                ed.select_hovered();
            } else {
                ed.escape();
            }
        }
        if !response.dragged() {
            ed.hover(model, &map, response.hover_pos());
        }
    }

    /// Routes viewport clicks to the profile picker.
    fn pick_input(&mut self, response: &egui::Response, rect: Rect) {
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

    fn paint_overlays(&mut self, painter: &egui::Painter, rect: Rect) {
        let theme = &self.theme;
        let view = ScreenMap::new(&self.camera, Frame::WORLD, rect);
        // Plane labels at a corner of each visible plane not seen edge-on.
        if self.sketch.is_none() {
            let h = self.plane_half();
            let forward = self.camera.forward();
            for p in StdPlane::ALL {
                if !self.planes_visible[p.index()] || planes::edge_on(p, forward) {
                    continue;
                }
                let corner = planes::corners(p, h)[3];
                if let Some(s) = view.world_to_screen(corner) {
                    let active = self.active_plane() == Some(p);
                    painter.text(
                        s + Vec2::new(4.0, 2.0),
                        Align2::LEFT_TOP,
                        p.label(),
                        FontId::proportional(theme.small_size),
                        if active {
                            theme.accent
                        } else {
                            theme.sketch.plane_label
                        },
                    );
                }
            }
        }
        // While picking profiles every sketch shows on top, even behind
        // parts; otherwise the 3D view draws them (`scene_meshes`).
        let editing = self.sketch.as_ref().map(|e| e.feature);
        if self.pick.is_some() {
            for s in self.model.sketches() {
                if Some(s.feature) != editing {
                    let d = SketchDisplay::from_result(s);
                    sketch_paint::paint_inactive(painter, &view.with_frame(s.frame), &d, theme);
                }
            }
        }
        if let Some(p) = self.pick.as_mut() {
            let mut regions: Vec<(FeatureId, usize, Color32)> = Vec::new();
            if let Some(f) = p.sketch {
                regions.extend(
                    p.picked(&self.model)
                        .into_iter()
                        .map(|i| (f, i, theme.sketch.region)),
                );
            }
            if let Some((f, i)) = p.hover {
                if !regions.iter().any(|r| r.0 == f && r.1 == i) {
                    regions.push((f, i, theme.sketch.region_hover));
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
                            painter.add(egui::Shape::line(run, Stroke::new(2.5, theme.selection)));
                        }
                    }
                }
            }
        }
        if let Some(ed) = &self.sketch {
            if let (Some(d), Some(frame)) = (ed.display(&self.model), ed.frame(&self.model)) {
                sketch_paint::paint_editor(painter, &view.with_frame(frame), &d, ed, theme);
            }
        }
    }

    fn paint_message(&mut self, painter: &egui::Painter, rect: Rect) {
        let Some((text, at)) = &self.message else {
            return;
        };
        if self.time - at > MESSAGE_SECONDS {
            self.message = None;
            return;
        }
        let theme = &self.theme;
        let font = FontId::proportional(theme.body_size);
        let galley = painter.layout(text.clone(), font, theme.text, 420.0);
        let size = galley.size() + Vec2::new(20.0, 12.0);
        let r = Rect::from_min_size(
            Pos2::new(rect.left() + 12.0, rect.bottom() - 12.0 - size.y),
            size,
        );
        painter.rect_filled(r, theme.rounding, theme.panel);
        painter.rect_stroke(
            r,
            theme.rounding,
            Stroke::new(1.0, theme.border),
            egui::StrokeKind::Inside,
        );
        painter.galley(r.min + Vec2::new(10.0, 6.0), galley, theme.text);
        // Keep repainting until it fades.
        painter
            .ctx()
            .request_repaint_after(std::time::Duration::from_millis(500));
    }
}
