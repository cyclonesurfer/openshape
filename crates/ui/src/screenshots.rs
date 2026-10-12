//! The application rendered offscreen: [`AppDriver`] runs the real
//! [`OpenShapeApp`] in a [`UiSnapshot`] and feeds it synthetic pointer and
//! keyboard input, for end-to-end tests of the shell and for the review
//! screenshots that [`run`] writes (`openshape --screenshots <dir>`).

use std::path::{Path, PathBuf};

use egui::{Event, Key, Modifiers, PointerButton, Pos2, Vec2};
use openshape_core_model::sketch::Geometry;
use openshape_core_model::{DVec2, DVec3, FeatureId, PlaneRef};
use openshape_kernel::Frame;
use openshape_render::{write_png, Camera, Projection, RgbaImage};
use openshape_ui_kit::UiSnapshot;

use crate::app::{OpenShapeApp, Selection};
use crate::dialogs::{FeatureDialog, Template};
use crate::profile_pick::{ProfilePick, SweepKind};
use crate::sketch_mode::{normal_view, ScreenMap, SketchEditor, Tool};
use crate::Model;

/// Seconds between synthetic frames.
const FRAME: f64 = 1.0 / 60.0;

/// The application running offscreen with synthetic input.
pub struct AppDriver {
    snap: UiSnapshot,
    rs: egui_wgpu::RenderState,
    /// The application.
    pub app: OpenShapeApp,
    time: f64,
    pointer: Pos2,
}

impl AppDriver {
    /// The app on `model` in a `size` pixel screen at `pixels_per_point`.
    pub fn new(size: [u32; 2], pixels_per_point: f32, model: Model) -> Result<Self, String> {
        let snap = UiSnapshot::new(size, pixels_per_point).map_err(|e| e.to_string())?;
        let rs = snap.render_state();
        let app = OpenShapeApp::with_render_state(snap.context(), Some(&rs), model);
        let mut d = AppDriver {
            snap,
            rs,
            app,
            time: 0.0,
            pointer: Pos2::new(-10.0, -10.0),
        };
        // egui sizes some widgets from the previous frame.
        d.settle()?;
        Ok(d)
    }

    /// The adapter the app renders with.
    pub fn adapter(&self) -> String {
        let info = self.snap.adapter_info();
        format!("{} ({:?})", info.name, info.backend)
    }

    /// Runs one frame with `events` and returns the image.
    pub fn frame_with(&mut self, events: Vec<Event>) -> Result<RgbaImage, String> {
        self.time += FRAME;
        let input = egui::RawInput {
            time: Some(self.time),
            predicted_dt: FRAME as f32,
            events,
            focused: true,
            ..Default::default()
        };
        let rs = self.rs.clone();
        let app = &mut self.app;
        self.snap
            .frame_ui(input, |ui| app.draw(ui, Some(&rs)))
            .map_err(|e| e.to_string())
    }

    /// Runs one frame without input.
    pub fn frame(&mut self) -> Result<RgbaImage, String> {
        self.frame_with(Vec::new())
    }

    /// Runs frames until layout has settled and popups have faded in
    /// (egui's default animation time is 1/12 s).
    pub fn settle(&mut self) -> Result<RgbaImage, String> {
        for _ in 0..7 {
            self.frame()?;
        }
        self.frame()
    }

    /// Moves the pointer.
    pub fn hover(&mut self, p: Pos2) -> Result<RgbaImage, String> {
        self.pointer = p;
        self.frame_with(vec![Event::PointerMoved(p)])
    }

    fn button(
        &mut self,
        p: Pos2,
        button: PointerButton,
        pressed: bool,
    ) -> Result<RgbaImage, String> {
        self.pointer = p;
        self.frame_with(vec![
            Event::PointerMoved(p),
            Event::PointerButton {
                pos: p,
                button,
                pressed,
                modifiers: Modifiers::NONE,
            },
        ])
    }

    /// Lets a second pass, so egui does not count the next click as part of
    /// a double or triple click with earlier ones (the frames of a script
    /// are only 1/60 s apart).
    fn pause(&mut self) {
        self.time += 1.0;
    }

    /// A click of `button` at `p`.
    pub fn click_with(&mut self, p: Pos2, button: PointerButton) -> Result<RgbaImage, String> {
        self.pause();
        self.hover(p)?;
        self.button(p, button, true)?;
        self.button(p, button, false)?;
        self.frame()
    }

    /// A primary click at `p`.
    pub fn click(&mut self, p: Pos2) -> Result<RgbaImage, String> {
        self.click_with(p, PointerButton::Primary)
    }

    /// A secondary (context-menu) click at `p`.
    pub fn right_click(&mut self, p: Pos2) -> Result<RgbaImage, String> {
        self.click_with(p, PointerButton::Secondary)
    }

    /// A primary double-click at `p`.
    pub fn double_click(&mut self, p: Pos2) -> Result<RgbaImage, String> {
        self.pause();
        self.hover(p)?;
        for _ in 0..2 {
            self.button(p, PointerButton::Primary, true)?;
            self.button(p, PointerButton::Primary, false)?;
        }
        self.frame()
    }

    /// A primary drag from `from` to `to` in `steps` moves.
    pub fn drag(&mut self, from: Pos2, to: Pos2, steps: usize) -> Result<RgbaImage, String> {
        self.pause();
        self.hover(from)?;
        self.button(from, PointerButton::Primary, true)?;
        let steps = steps.max(1);
        for i in 1..=steps {
            let p = from + (to - from) * (i as f32 / steps as f32);
            self.hover(p)?;
        }
        self.button(to, PointerButton::Primary, false)?;
        self.frame()
    }

    /// Presses and releases a key.
    pub fn key(&mut self, key: Key, modifiers: Modifiers) -> Result<RgbaImage, String> {
        let ev = |pressed| Event::Key {
            key,
            physical_key: None,
            pressed,
            repeat: false,
            modifiers,
        };
        self.frame_with(vec![ev(true)])?;
        self.frame_with(vec![ev(false)])
    }

    /// Types text into the focused field.
    pub fn type_text(&mut self, text: &str) -> Result<RgbaImage, String> {
        self.frame_with(vec![Event::Text(text.to_string())])
    }

    /// The pointer's last position.
    pub fn pointer(&self) -> Pos2 {
        self.pointer
    }

    /// The row of an open menu labelled `label`, as drawn in the last frame.
    pub fn menu_row(&self, label: &str) -> Option<egui::Rect> {
        crate::shell::menu_rows(self.snap.context())
            .into_iter()
            .find(|(l, _)| l == label)
            .map(|(_, r)| r)
    }
}

const MM: f64 = 1e-3;

/// A top view of ±150 mm for scripting sketches with screen clicks.
fn sketch_view(frame: Frame) -> ScreenMap {
    let (yaw, pitch) = normal_view(&frame);
    let camera = Camera {
        target: DVec3::ZERO,
        distance: 0.12 / (openshape_render::DEFAULT_FOV_Y_RAD as f64 * 0.5).tan(),
        yaw,
        pitch,
        projection: Projection::Orthographic,
        scene_bounds: None,
    };
    ScreenMap::new(
        &camera,
        frame,
        egui::Rect::from_min_size(Pos2::ZERO, Vec2::new(1000.0, 800.0)),
    )
}

fn at(map: &ScreenMap, x: f64, y: f64) -> Pos2 {
    map.to_screen(DVec2::new(x, y) * MM).expect("on screen")
}

fn click(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, x: f64, y: f64) {
    let p = at(map, x, y);
    ed.hover(m, map, Some(p));
    ed.click(m, map, p, false);
}

fn dimension(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, x: f64, y: f64, value: &str) {
    ed.dimension(m, at(map, x, y));
    if let Some(p) = ed.prompt.as_mut() {
        p.text = value.into();
    }
    ed.commit_prompt(m);
}

/// The document the screenshots show: a plate with a hole (sketch and
/// extrude), a cylinder and a vase, and a second, empty Part Studio.
pub fn demo_model() -> Result<(Model, FeatureId, FeatureId), String> {
    let mut m = Model::new();
    let sketch = m.add_sketch(PlaneRef::Top);
    let map = sketch_view(Frame::top());
    let mut ed = SketchEditor::new(sketch);
    ed.set_tool(Tool::Rectangle);
    click(&mut m, &mut ed, &map, -40.0, -25.0);
    click(&mut m, &mut ed, &map, 40.0, 25.0);
    ed.set_tool(Tool::Circle);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    click(&mut m, &mut ed, &map, 12.0, 0.0);
    ed.set_tool(Tool::Select);
    // Width, height and the hole's diameter.
    click(&mut m, &mut ed, &map, 0.0, -25.0);
    dimension(&mut m, &mut ed, &map, 0.0, -38.0, "80");
    click(&mut m, &mut ed, &map, 40.0, 0.0);
    dimension(&mut m, &mut ed, &map, 52.0, 0.0, "50");
    click(&mut m, &mut ed, &map, 0.0, 12.0);
    dimension(&mut m, &mut ed, &map, 16.0, 16.0, "24");
    ed.escape();
    let mut pick = ProfilePick::new(SweepKind::Extrude);
    let view = sketch_view(Frame::WORLD);
    pick.click(&m, &view, at(&view, 30.0, 15.0));
    pick.distance = "12".into();
    let extrude = pick.build()?;
    let extrude = m.add_feature("Extrude 1", extrude);
    let cylinder = FeatureDialog::new(Template::Cylinder).build()?;
    m.add_feature("Cylinder 1", cylinder);
    m.add_feature("Revolve 1", crate::demo::vase_at(-0.075));
    m.add_part_studio();
    let first = m.elements()[0].id;
    m.set_active(first);
    Ok((m, sketch, extrude))
}

/// The ids of a sketch's line entities.
fn lines(m: &Model, sketch: FeatureId) -> Vec<openshape_core_model::SketchId> {
    m.sketch_def(sketch)
        .map(|d| {
            d.sketch
                .entities
                .iter()
                .filter(|e| matches!(e.geometry, Geometry::Line { .. }))
                .map(|e| e.id)
                .collect()
        })
        .unwrap_or_default()
}

/// Writes the review screenshots of the shell into `out` and returns their
/// paths.
pub fn run(out: &Path) -> Result<Vec<PathBuf>, String> {
    std::fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let size = [1600, 1000];
    let mut written = Vec::new();
    let mut save = |name: &str, img: &RgbaImage| -> Result<(), String> {
        let p = out.join(name);
        write_png(&p, img).map_err(|e| e.to_string())?;
        written.push(p);
        Ok(())
    };

    // A new document: the planes and the empty lists.
    let mut d = AppDriver::new(size, 1.0, Model::new())?;
    save("01-new-document.png", &d.settle()?)?;

    // The demo document with a part selected.
    let (model, sketch, extrude) = demo_model()?;
    let mut d = AppDriver::new(size, 1.0, model)?;
    d.app.set_document_name("Bracket study");
    d.app.select(Some(Selection::Part(extrude)));
    save("02-part-studio.png", &d.settle()?)?;

    // A feature's context menu.
    if let Some(r) = d.app.layout().feature_row(extrude) {
        d.right_click(r.center())?;
        save("03-context-menu.png", &d.settle()?)?;
        d.key(Key::Escape, Modifiers::NONE)?;
        d.click(d.app.layout().graphics.left_top() + Vec2::new(300.0, 500.0))?;
    }

    // Editing the extrude: its dialog and the picked region.
    d.app.edit_feature(extrude);
    save("04-extrude-dialog.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;

    // Sketch mode with a selected line.
    d.app.enter_sketch(sketch);
    let first_line = lines(d.app.model(), sketch).first().copied();
    if let (Some(ed), Some(line)) = (d.app.sketch_editor_mut(), first_line) {
        ed.selection = vec![line];
    }
    save("05-sketch-mode.png", &d.settle()?)?;
    d.app.run(crate::commands::Command::FinishSketch);

    // The history rolled back to the sketch.
    d.app.select(None);
    d.app.roll_to_gap(1);
    d.app
        .set_view(openshape_ui_kit::view_cube::CubeTarget::HOME);
    save("06-rolled-back.png", &d.settle()?)?;
    d.app.roll_to_gap(usize::MAX);

    // The command palette with a query.
    d.app
        .set_view(openshape_ui_kit::view_cube::CubeTarget::HOME);
    d.key(Key::K, Modifiers::COMMAND)?;
    d.type_text("ex")?;
    save("07-command-palette.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;

    // Trimming: a line across the plate, and the piece under the pointer
    // that a click would remove.
    d.app.enter_sketch(sketch);
    d.settle()?;
    let on_plane = |d: &AppDriver, x: f64, y: f64| {
        ScreenMap::new(&d.app.camera, Frame::top(), d.app.layout().graphics)
            .to_screen(DVec2::new(x, y) * 1e-3)
            .ok_or("the sketch plane is off screen")
    };
    d.app.run(crate::commands::Command::Tool(Tool::Line));
    d.click(on_plane(&d, -60.0, 0.0)?)?;
    d.click(on_plane(&d, 60.0, 0.0)?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    d.app.run(crate::commands::Command::Tool(Tool::Trim));
    d.hover(on_plane(&d, 26.0, 0.0)?)?;
    save("08-trim-preview.png", &d.settle()?)?;
    // A fillet of the plate's corner, previewed.
    d.app
        .run(crate::commands::Command::Tool(Tool::SketchFillet));
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.fillet_radius = 0.008;
    }
    d.hover(on_plane(&d, 40.0, 25.0)?)?;
    save("09-fillet-preview.png", &d.settle()?)?;
    // An offset of the plate's outline, chosen by its bottom side and
    // previewed outwards.
    d.app.run(crate::commands::Command::Tool(Tool::Offset));
    d.click(on_plane(&d, 0.0, -25.0)?)?;
    d.hover(on_plane(&d, 10.0, -31.0)?)?;
    save("10-offset-preview.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    // A triangle mirrored about the line across the plate, previewed.
    d.app.run(crate::commands::Command::Tool(Tool::Line));
    for (x, y) in [(-34.0, 6.0), (-14.0, 6.0), (-24.0, 18.0), (-34.0, 6.0)] {
        d.click(on_plane(&d, x, y)?)?;
    }
    d.key(Key::Escape, Modifiers::NONE)?;
    let triangle: Vec<_> = lines(d.app.model(), sketch)
        .into_iter()
        .rev()
        .take(3)
        .collect();
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.selection = triangle.clone();
    }
    d.app.run(crate::commands::Command::Tool(Tool::Mirror));
    d.hover(on_plane(&d, -55.0, 0.0)?)?;
    save("11-mirror-preview.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    // The triangle patterned downwards, previewed with the step's end
    // under the pointer.
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.selection = triangle;
    }
    d.app
        .run(crate::commands::Command::Tool(Tool::LinearPattern));
    d.click(on_plane(&d, -34.0, 6.0)?)?;
    d.hover(on_plane(&d, -34.0, -12.0)?)?;
    save("12-linear-pattern-preview.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    // A small circle turned about the plate's centre: previewed, then
    // made, with the angle open for typing.
    d.app.run(crate::commands::Command::Tool(Tool::Circle));
    d.click(on_plane(&d, 0.0, 18.0)?)?;
    d.click(on_plane(&d, 3.0, 18.0)?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    let small = d.app.model().sketch_def(sketch).and_then(|def| {
        def.sketch
            .entities
            .iter()
            .rev()
            .find(|e| matches!(e.geometry, Geometry::Circle { .. }))
            .map(|e| e.id)
    });
    if let (Some(ed), Some(c)) = (d.app.sketch_editor_mut(), small) {
        ed.selection = vec![c];
        ed.pattern_count = 4;
    }
    d.app
        .run(crate::commands::Command::Tool(Tool::CircularPattern));
    d.hover(on_plane(&d, 0.0, 0.0)?)?;
    save("13-circular-pattern-preview.png", &d.settle()?)?;
    d.click(on_plane(&d, 0.0, 0.0)?)?;
    save("14-circular-pattern.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    // Dimension labels dragged to where the user wants them: the width
    // further below, the height out to the right, the pattern's angle out
    // past the copies.
    for (text, x, y) in [("80", -22.0, -42.0), ("50", 60.0, 8.0), ("90°", 24.0, 24.0)] {
        let from = label_centre(&mut d, sketch, text)?;
        d.drag(from, on_plane(&d, x, y)?, 6)?;
    }
    d.app.select(None);
    d.hover(on_plane(&d, 120.0, -80.0)?)?;
    save("15-placed-labels.png", &d.settle()?)?;
    d.app.run(crate::commands::Command::CancelSketch);
    // Editing the extrude: a new depth typed and previewed before OK.
    d.app.edit_feature(extrude);
    if let Some(p) = d.app.profile_pick_mut() {
        p.distance = "30".into();
    }
    d.app
        .set_view(openshape_ui_kit::view_cube::CubeTarget::HOME);
    save("16-extrude-preview.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    // A control-point spline left of the plate: being drawn (its control
    // polygon and the curve to the pointer), then finished.
    d.app.enter_sketch(sketch);
    d.settle()?;
    d.app
        .run(crate::commands::Command::Tool(Tool::ControlSpline));
    for (x, y) in [(-80.0, -38.0), (-76.0, 18.0), (-64.0, -30.0), (-56.0, 12.0)] {
        d.click(on_plane(&d, x, y)?)?;
    }
    d.hover(on_plane(&d, -48.0, -26.0)?)?;
    save("17-control-spline-drawing.png", &d.settle()?)?;
    d.click(on_plane(&d, -48.0, -26.0)?)?;
    d.key(Key::Enter, Modifiers::NONE)?;
    d.hover(on_plane(&d, 120.0, -80.0)?)?;
    save("18-control-spline.png", &d.settle()?)?;
    // Right of the plate: an ellipse at an angle, then an elliptical arc
    // being drawn the long way round (the pointer turned from its start
    // through the top), then finished.
    d.app.run(crate::commands::Command::Tool(Tool::Ellipse));
    for (x, y) in [(62.0, 24.0), (75.0, 30.0), (60.0, 32.0)] {
        d.click(on_plane(&d, x, y)?)?;
    }
    d.app
        .run(crate::commands::Command::Tool(Tool::EllipticalArc));
    for (x, y) in [(62.0, -28.0), (76.0, -28.0), (62.0, -19.0)] {
        d.click(on_plane(&d, x, y)?)?;
    }
    let on = |deg: f64| {
        let t = deg.to_radians();
        (62.0 + 14.0 * t.cos(), -28.0 + 9.0 * t.sin())
    };
    let (x, y) = on(-30.0);
    d.click(on_plane(&d, x, y)?)?;
    for deg in [0.0, 40.0, 80.0, 120.0, 150.0] {
        let (x, y) = on(deg);
        d.hover(on_plane(&d, x, y)?)?;
    }
    save("19-elliptical-arc-drawing.png", &d.settle()?)?;
    for deg in [180.0, 210.0] {
        let (x, y) = on(deg);
        d.hover(on_plane(&d, x, y)?)?;
    }
    let (x, y) = on(210.0);
    d.click(on_plane(&d, x, y)?)?;
    d.hover(on_plane(&d, 120.0, -80.0)?)?;
    save("20-ellipses.png", &d.settle()?)?;
    // Below the plate: a conic being drawn (its ends, then the corner where
    // its end tangents meet; the curve reaches up to the pointer), then a
    // second one, flatter, with its ρ dimensioned.
    d.app.run(crate::commands::Command::Tool(Tool::Conic));
    for (x, y) in [(-34.0, -50.0), (-6.0, -50.0), (-24.0, -33.0)] {
        d.click(on_plane(&d, x, y)?)?;
    }
    // ρ = 0.7: 0.7 of the corner's height above the chord.
    d.hover(on_plane(&d, -12.0, -38.1)?)?;
    save("21-conic-drawing.png", &d.settle()?)?;
    d.click(on_plane(&d, -12.0, -38.1)?)?;
    for (x, y) in [(6.0, -50.0), (34.0, -50.0), (20.0, -33.0), (26.0, -44.05)] {
        d.click(on_plane(&d, x, y)?)?;
    }
    d.app.run(crate::commands::Command::Tool(Tool::Select));
    // Selected at its shoulder (ρ = 0.35), dimensioned, and set to 0.42.
    d.click(on_plane(&d, 20.0, -44.05)?)?;
    d.hover(on_plane(&d, 30.0, -36.0)?)?;
    d.app.run(crate::commands::Command::Dimension);
    if let Some(p) = d.app.sketch_editor_mut().and_then(|ed| ed.prompt.as_mut()) {
        p.text = "0.42".into();
    }
    d.settle()?;
    d.key(Key::Enter, Modifiers::NONE)?;
    d.app.select(None);
    d.hover(on_plane(&d, 120.0, -80.0)?)?;
    save("22-conics.png", &d.settle()?)?;
    // The conics turned about the first one's start: previewed with the
    // reference at its end and the pointer 40° up, the dialog showing the
    // angle and the options.
    let conics: Vec<_> = d
        .app
        .model()
        .sketch_def(sketch)
        .map(|def| {
            def.sketch
                .entities
                .iter()
                .filter(|e| matches!(e.geometry, Geometry::Conic { .. }))
                .map(|e| e.id)
                .collect()
        })
        .unwrap_or_default();
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.selection = conics.clone();
    }
    d.app.run(crate::commands::Command::Tool(Tool::Rotate));
    d.click(on_plane(&d, 6.0, -50.0)?)?;
    d.click(on_plane(&d, 34.0, -50.0)?)?;
    let turn = 40f64.to_radians();
    d.hover(on_plane(
        &d,
        6.0 + 28.0 * turn.cos(),
        -50.0 + 28.0 * turn.sin(),
    )?)?;
    save("23-rotate-preview.png", &d.settle()?)?;
    d.key(Key::Escape, Modifiers::NONE)?;
    // Copied above the plate, with their ρ dimension.
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.selection = conics;
    }
    d.app.run(crate::commands::Command::Tool(Tool::Move));
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.transform_copy = true;
    }
    d.click(on_plane(&d, 6.0, -50.0)?)?;
    d.click(on_plane(&d, 6.0, 32.0)?)?;
    d.hover(on_plane(&d, 120.0, -80.0)?)?;
    save("24-moved-copy.png", &d.settle()?)?;
    // Bottom left: an arc turning down to its end, a spline drawn on from
    // it bending the other way, then the two made curvature continuous (the
    // spline turned round, the arc unchanged, the joint marked C).
    d.app.run(crate::commands::Command::Tool(Tool::CenterArc));
    d.click(on_plane(&d, -70.0, -44.0)?)?;
    d.click(on_plane(&d, -76.0, -44.0)?)?;
    d.hover(on_plane(&d, -74.2, -48.2)?)?;
    d.click(on_plane(&d, -70.0, -50.0)?)?;
    d.app
        .run(crate::commands::Command::Tool(Tool::ControlSpline));
    for (x, y) in [
        (-70.0, -50.0),
        (-63.0, -48.0),
        (-57.0, -52.0),
        (-50.0, -47.0),
    ] {
        d.click(on_plane(&d, x, y)?)?;
    }
    d.key(Key::Enter, Modifiers::NONE)?;
    d.app.run(crate::commands::Command::Tool(Tool::Select));
    let joined: Vec<_> = d
        .app
        .model()
        .sketch_def(sketch)
        .map(|def| {
            let last = |pick: fn(&Geometry<_>) -> bool| {
                def.sketch
                    .entities
                    .iter()
                    .rev()
                    .find(|e| pick(&e.geometry))
                    .map(|e| e.id)
            };
            [
                last(|g| matches!(g, Geometry::ControlSpline { .. })),
                last(|g| matches!(g, Geometry::Arc { .. })),
            ]
            .into_iter()
            .flatten()
            .collect()
        })
        .unwrap_or_default();
    if let Some(ed) = d.app.sketch_editor_mut() {
        ed.selection = joined;
    }
    d.app.run(crate::commands::Command::Constrain(
        crate::sketch_mode::ConstraintTool::Curvature,
    ));
    d.app.select(None);
    d.hover(on_plane(&d, 120.0, -80.0)?)?;
    save("25-curvature.png", &d.settle()?)?;
    d.app.run(crate::commands::Command::CancelSketch);

    // Text in a new part: being placed (the box growing towards the
    // pointer, the letters in it), then selected with its characters and
    // font in the dialog, then its letters extruded.
    let mut model = Model::new();
    let sketch = model.add_sketch(PlaneRef::Top);
    let mut t = AppDriver::new(size, 1.0, model)?;
    t.app.enter_sketch(sketch);
    t.settle()?;
    t.app.run(crate::commands::Command::Tool(Tool::Text));
    if let Some(ed) = t.app.sketch_editor_mut() {
        ed.text_content = "OpenShape".into();
        ed.text_font = openshape_text::Font::SansBold;
    }
    t.click(on_plane(&t, -60.0, -8.0)?)?;
    t.hover(on_plane(&t, -30.0, 12.0)?)?;
    save("26-text-placing.png", &t.settle()?)?;
    t.click(on_plane(&t, -30.0, 12.0)?)?;
    t.app.run(crate::commands::Command::Tool(Tool::Select));
    let text = t
        .app
        .model()
        .sketch_def(sketch)
        .and_then(|def| {
            def.sketch
                .entities
                .iter()
                .find(|e| matches!(e.geometry, Geometry::Text { .. }))
                .map(|e| e.id)
        })
        .ok_or("the text was not placed")?;
    if let Some(ed) = t.app.sketch_editor_mut() {
        ed.selection = vec![text];
    }
    t.hover(on_plane(&t, 90.0, -60.0)?)?;
    save("27-text.png", &t.settle()?)?;
    t.app.run(crate::commands::Command::FinishSketch);
    // Every letter's ink (not the counters inside O, p, a and e).
    let interiors: Vec<DVec2> = t
        .app
        .model()
        .sketch_result(sketch)
        .and_then(|r| r.regions.as_ref().ok())
        .map(|regions| {
            let a = &regions.arrangement;
            let edges = |l: &openshape_kernel::RegionLoop| -> Vec<usize> {
                let mut e: Vec<usize> = l.edges.iter().map(|(e, _)| *e).collect();
                e.sort_unstable();
                e
            };
            let holes: Vec<Vec<usize>> = a
                .regions
                .iter()
                .flat_map(|r| r.holes.iter().map(edges))
                .collect();
            a.regions
                .iter()
                .filter(|r| !holes.contains(&edges(&r.outer)))
                .map(|r| r.interior_point)
                .collect()
        })
        .unwrap_or_default();
    let view = ScreenMap::new(&t.app.camera, Frame::top(), t.app.layout().graphics);
    let mut pick = ProfilePick::new(SweepKind::Extrude);
    for p in interiors {
        let at = view.to_screen(p).ok_or("a letter is off screen")?;
        pick.click(t.app.model(), &view, at);
    }
    pick.distance = "4".into();
    let letters = pick.build()?;
    t.app.model_mut().add_feature("Extrude 1", letters);
    // Closer than "fit", which frames the default planes too.
    t.app
        .set_view(openshape_ui_kit::view_cube::CubeTarget::HOME);
    t.app.camera.target = DVec3::new(-0.011, 0.002, 0.002);
    t.app.camera.distance *= 0.4;
    save("28-text-extruded.png", &t.settle()?)?;
    Ok(written)
}

/// The screen centre of the dimension label of `sketch` (being edited)
/// that reads `text`.
fn label_centre(d: &mut AppDriver, sketch: FeatureId, text: &str) -> Result<Pos2, String> {
    let rect = d.app.layout().graphics;
    let camera = d.app.camera.clone();
    let model = d.app.model();
    let frame = model
        .sketch_result(sketch)
        .ok_or("the sketch does not solve")?
        .frame;
    let shown = d
        .app
        .sketch_editor()
        .and_then(|ed| ed.display(model))
        .ok_or("not editing the sketch")?;
    crate::sketch_mode::annotations(&shown, &ScreenMap::new(&camera, frame, rect))
        .into_iter()
        .find(|a| a.dimension && a.text.starts_with(text))
        .map(|a| a.rect.center())
        .ok_or_else(|| format!("no dimension reads {text}"))
}
