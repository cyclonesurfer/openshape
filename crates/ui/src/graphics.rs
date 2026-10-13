//! The graphics area (ADR 0012): the 3D view rendered offscreen and shown
//! through egui, the default planes, picking faces, edges, vertices and
//! planes (ADR 0019) with their highlights, sketches and sketch-mode
//! feedback, the view cube, the graphics-area context menu and messages.

use eframe::egui;
use egui::{Align2, Color32, FontId, PointerButton, Pos2, Rect, Sense, Stroke, Ui, Vec2};
use opendatum_core_model::{DVec2, FeatureId};
use opendatum_kernel::Frame;
use opendatum_naming::Pick;
use opendatum_render::{wgpu, DisplayMode, MeshData, SceneRenderer};
use opendatum_ui_kit::{view_cube, Theme};

use crate::app::{Action, OpenDatumApp, Selection};
use crate::commands::Command;
use crate::entity_pick::{logical_pieces, pick_entity, EntityHit, PickKinds};
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
    /// Highlighted faces, edges and vertices, with their colours' roles.
    highlights: Vec<(FeatureId, Pick, Highlight)>,
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

/// How far a selected or hovered face's colour moves towards the
/// selection or preselection colour.
const FACE_TINT: f32 = 0.55;

/// Why an entity is highlighted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Highlight {
    /// Selected (or the face an open dialog uses).
    Selected,
    /// Under the pointer.
    Hovered,
}

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

impl OpenDatumApp {
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
        let highlights = self.highlights();
        let color = |h: Highlight| match h {
            Highlight::Selected => self.theme.selection,
            Highlight::Hovered => self.theme.preselect,
        };
        let mut out: Vec<MeshData> = Vec::new();
        let parts = self
            .model
            .meshes()
            .iter()
            .zip(self.model.body_meshes())
            .zip(self.model.mesh_owners());
        for ((m, body), f) in parts {
            if self.hidden_parts.contains(f) {
                continue;
            }
            let m = if Some(*f) == selected || Some(*f) == previewed {
                MeshData {
                    color: tint(m.color, self.theme.selection, SELECTION_TINT),
                    edge_color: Some(rgba(self.theme.selection)),
                    ..m.clone()
                }
            } else {
                m.clone()
            };
            // Highlighted faces are drawn in their own meshes, so no face
            // is drawn twice; the part's edges stay with the rest.
            let faces: Vec<(usize, Highlight)> = highlights
                .iter()
                .filter(|(g, _, _)| g == f)
                .filter_map(|(_, p, h)| match p {
                    Pick::Face(ix) => Some((ix.0 as usize, *h)),
                    _ => None,
                })
                .collect();
            if faces.is_empty() {
                out.push(m);
            } else {
                let lit = |i: usize| faces.iter().find(|(j, _)| *j == i).map(|(_, h)| *h);
                out.push(MeshData {
                    edge_color: m.edge_color,
                    ..crate::to_mesh_data_where(body, m.color, |i| lit(i).is_none(), true)
                });
                for h in [Highlight::Hovered, Highlight::Selected] {
                    if faces.iter().any(|(_, k)| *k == h) {
                        out.push(crate::to_mesh_data_where(
                            body,
                            tint(m.color, color(h), FACE_TINT),
                            |i| lit(i) == Some(h),
                            false,
                        ));
                    }
                }
            }
            // Highlighted edges, drawn again over the part's in their
            // colour.
            for h in [Highlight::Hovered, Highlight::Selected] {
                let mut lines = Vec::new();
                for (_, p, _) in highlights.iter().filter(|(g, _, k)| g == f && *k == h) {
                    if let Pick::Edge(e) = p {
                        if let Some(line) = body.edges.get(e.0 as usize) {
                            for w in line.windows(2) {
                                lines.push(w[0].as_vec3().to_array());
                                lines.push(w[1].as_vec3().to_array());
                            }
                        }
                    }
                }
                if !lines.is_empty() {
                    out.push(MeshData {
                        edge_lines: lines,
                        edge_color: Some(rgba(color(h))),
                        ..MeshData::default()
                    });
                }
            }
        }
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

    /// The kernel faces, edges and vertices to highlight: every piece of
    /// the selected entity, of the one under the pointer, and of the face
    /// an open extrude ends at. A hovered entity that is also selected
    /// shows as selected.
    fn highlights(&self) -> Vec<(FeatureId, Pick, Highlight)> {
        let mut out: Vec<(FeatureId, Pick, Highlight)> = Vec::new();
        let mut add = |f: FeatureId, pieces: Vec<Pick>, h: Highlight| {
            for p in pieces {
                if !out.iter().any(|(g, q, _)| *g == f && *q == p) {
                    out.push((f, p, h));
                }
            }
        };
        let pieces = |hit: &EntityHit| {
            let names = self
                .model
                .part(hit.feature)
                .and_then(|p| p.names.as_ref().ok());
            logical_pieces(names, hit.pick)
        };
        if let Some(Selection::Entity(hit)) = &self.selection {
            add(hit.feature, pieces(hit), Highlight::Selected);
        }
        if let Some((f, picks)) = self
            .pick
            .as_ref()
            .and_then(|p| p.up_to.as_ref())
            .and_then(|r| self.model.resolve(r))
        {
            add(f, picks, Highlight::Selected);
        }
        if let Some(hit) = &self.hover_entity {
            add(hit.feature, pieces(hit), Highlight::Hovered);
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
            highlights: self.highlights(),
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

    /// The face, edge or vertex of a shown part under a screen point of
    /// the graphics area (ADR 0019).
    pub(crate) fn entity_at(&self, rect: Rect, pos: Pos2, kinds: PickKinds) -> Option<EntityHit> {
        let parts = self.model.pick_parts(|f| self.hidden_parts.contains(&f));
        pick_entity(&parts, &self.camera, rect, pos, kinds)
    }

    /// The face, edge, vertex or plane under a screen point of the
    /// graphics area. Parts hide the planes behind them.
    fn pick_at(&self, rect: Rect, pos: Pos2) -> Option<Selection> {
        if let Some(hit) = self.entity_at(rect, pos, PickKinds::All) {
            return Some(Selection::Entity(hit));
        }
        let (o, d) = self.camera.pixel_ray(
            pos.x - rect.min.x,
            pos.y - rect.min.y,
            rect.width(),
            rect.height(),
        );
        let h = self.plane_half();
        StdPlane::ALL
            .into_iter()
            .filter(|p| self.planes_visible[p.index()])
            .filter_map(|p| planes::hit(p, h, o, d).map(|t| (p, t)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(p, _)| Selection::Plane(p))
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
        // Faces, edges and vertices under the pointer light up outside
        // sketch mode; while an extrude picks the face it ends at, faces
        // only.
        let kinds = match &self.pick {
            None => Some(PickKinds::All),
            Some(p) if p.picking_face => Some(PickKinds::Faces),
            Some(_) => None,
        };
        let hover = match (sketching, kinds, response.hover_pos()) {
            (false, Some(k), Some(pos)) if !response.dragged() => self.entity_at(rect, pos, k),
            _ => None,
        };
        self.hover_entity = hover;
        self.pick_input(&response, rect);
        if !sketching && !picking_profile {
            if response.clicked_by(PointerButton::Primary) {
                if let Some(pos) = response.interact_pointer_pos() {
                    self.actions.push(Action::Select(self.pick_at(rect, pos)));
                }
            }
            if response.double_clicked_by(PointerButton::Primary) {
                if let Some(Selection::Entity(hit)) = response
                    .interact_pointer_pos()
                    .and_then(|p| self.pick_at(rect, p))
                {
                    self.actions.push(Action::EditFeature(hit.feature));
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
            Some(opendatum_ui_kit::Icon::Sketch),
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
            Some(opendatum_ui_kit::Icon::Plane),
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
                Some(opendatum_ui_kit::Icon::Eye),
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
            .set_edge_width(opendatum_render::DEFAULT_EDGE_WIDTH_PX * ppp);
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

    /// Routes viewport clicks to the profile picker: regions and axes, or
    /// the face an extrude ends at.
    fn pick_input(&mut self, response: &egui::Response, rect: Rect) {
        if self.pick.as_ref().is_some_and(|p| p.picking_face) {
            if response.clicked_by(PointerButton::Primary) {
                let hit = response
                    .interact_pointer_pos()
                    .and_then(|pos| self.entity_at(rect, pos, PickKinds::Faces));
                if let Some(hit) = hit {
                    let picked = self.model.planar_reference(&hit);
                    if let Some(p) = self.pick.as_mut() {
                        p.set_up_to(picked);
                    }
                }
            }
            return;
        }
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
        // Highlighted vertices: dots over the view.
        if self.sketch.is_none() {
            for (f, p, h) in self.highlights() {
                let Pick::Vertex(v) = p else { continue };
                let Some(q) = self
                    .model
                    .body_mesh(f)
                    .and_then(|m| m.vertices.get(v.0 as usize))
                else {
                    continue;
                };
                if let Some(s) = view.world_to_screen(*q) {
                    let c = match h {
                        Highlight::Selected => theme.selection,
                        Highlight::Hovered => theme.preselect,
                    };
                    painter.circle(s, 4.0, c, Stroke::new(1.5, theme.panel));
                }
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
