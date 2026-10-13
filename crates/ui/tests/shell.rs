//! End-to-end tests of the application shell (ADR 0012): the real app runs
//! offscreen with synthetic pointer and keyboard input ([`AppDriver`]) and
//! the tests click what a user would click, found through the shell's
//! layout record. They need a Vulkan or GL adapter (lavapipe in CI).

use egui::{Key, Modifiers, Pos2, Rect, Vec2};
use opendatum_app::app::Selection;
use opendatum_app::commands::{Command, Mode};
use opendatum_app::screenshots::{demo_model, AppDriver};
use opendatum_app::sketch_mode::{annotations, normal_view, ScreenMap, SketchDisplay};
use opendatum_app::{Model, PREVIEW_FEATURE};
use opendatum_core_model::{DVec3, Feature, FeatureId, InsertAt, ModelOp};
use opendatum_kernel::{Frame, Kernel};
use opendatum_render::{Camera, Projection, RgbaImage};

const SIZE: [u32; 2] = [1280, 800];

fn driver(model: Model) -> AppDriver {
    AppDriver::new(SIZE, 1.0, model).expect("an adapter (Vulkan or GL)")
}

fn demo() -> (AppDriver, FeatureId, FeatureId) {
    let (model, sketch, extrude) = demo_model().expect("demo model builds");
    (driver(model), sketch, extrude)
}

fn feature_ids(d: &AppDriver) -> Vec<FeatureId> {
    d.app
        .model()
        .studio()
        .features
        .iter()
        .map(|f| f.id)
        .collect()
}

/// Where a model-space point appears in the graphics area.
fn on_screen(d: &mut AppDriver, p: DVec3) -> Pos2 {
    let rect = d.app.layout().graphics;
    let camera = d.app.camera_mut().clone();
    ScreenMap::new(&camera, Frame::WORLD, rect)
        .world_to_screen(p)
        .expect("in front of the camera")
}

fn pixel(img: &RgbaImage, p: Pos2) -> [u8; 4] {
    img.pixel(p.x.round() as u32, p.y.round() as u32)
}

#[test]
fn the_shell_lays_out_its_regions() {
    let (mut d, sketch, extrude) = demo();
    let img = d.settle().unwrap();
    let l = d.app.layout().clone();
    let screen = Rect::from_min_size(Pos2::ZERO, Vec2::new(SIZE[0] as f32, SIZE[1] as f32));
    for c in [
        Command::Sketch,
        Command::Extrude,
        Command::Revolve,
        Command::Undo,
    ] {
        let r = l.command(c).unwrap_or_else(|| panic!("{c:?} has a button"));
        assert!(
            screen.contains_rect(r) && r.max.y < l.graphics.min.y,
            "{c:?} at {r:?}"
        );
    }
    assert!(l.graphics.width() > 800.0 && l.graphics.height() > 500.0);
    assert_eq!(l.feature_rows.len(), 4);
    assert_eq!(l.plane_rows.len(), 3);
    assert_eq!(l.part_rows.len(), 3);
    assert_eq!(l.tabs.len(), 2);
    assert!(l.feature_row(sketch).unwrap().max.y <= l.feature_row(extrude).unwrap().min.y);
    let bar = l.rollback_bar.expect("a rollback bar");
    assert!(bar.min.y >= l.feature_rows.last().unwrap().1.max.y - 1.0);
    let cube = l.view_cube.expect("a view cube");
    assert!(l.graphics.contains_rect(cube));
    // The panels are light, the graphics area is not blank.
    let panel = pixel(&img, Pos2::new(120.0, l.graphics.max.y - 20.0));
    assert!(panel[..3].iter().all(|v| *v > 235), "panel {panel:?}");
    let centre = pixel(&img, l.graphics.center());
    assert!(centre[..3] != [255, 255, 255], "graphics centre {centre:?}");
}

#[test]
fn toolbar_buttons_open_their_dialogs() {
    let (mut d, _, _) = demo();
    let extrude = d.app.layout().command(Command::Extrude).unwrap();
    d.click(extrude.center()).unwrap();
    d.settle().unwrap();
    assert!(
        d.app.profile_pick_mut().is_some(),
        "the extrude dialog is open"
    );
    let dialog = d.app.layout().dialog.expect("a floating dialog");
    assert!(d.app.layout().graphics.contains_rect(dialog));
    assert!(d.app.active(Command::Extrude));
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    assert!(d.app.profile_pick_mut().is_none());

    let cylinder = Command::Template(opendatum_app::dialogs::Template::Cylinder);
    let before = feature_ids(&d).len();
    d.app.run(cylinder);
    d.settle().unwrap();
    assert!(d.app.feature_dialog().is_some());
    d.key(Key::Enter, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert!(d.app.feature_dialog().is_none(), "Enter accepts");
    assert_eq!(feature_ids(&d).len(), before + 1);
    assert_eq!(d.app.model().meshes().len(), 4);
}

#[test]
fn a_new_sketch_from_a_plane_row_and_cancelling_it() {
    let mut d = driver(Model::new());
    let top = d
        .app
        .layout()
        .plane_rows
        .iter()
        .find(|(p, _)| *p == opendatum_app::planes::StdPlane::Top)
        .map(|(_, r)| *r)
        .expect("a Top row");
    d.double_click(top.center()).unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.mode(), Mode::Sketch);
    assert_eq!(
        feature_ids(&d).len(),
        1,
        "the sketch feature exists while editing"
    );
    assert!(d.app.layout().command(Command::FinishSketch).is_some());
    // Draw a line, then cancel: the sketch and the line are gone.
    d.key(Key::L, Modifiers::NONE).unwrap();
    let g = d.app.layout().graphics;
    d.click(g.center() - Vec2::new(100.0, 0.0)).unwrap();
    d.click(g.center() + Vec2::new(100.0, 50.0)).unwrap();
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    let id = feature_ids(&d)[0];
    let entities = d.app.model().sketch_def(id).unwrap().sketch.entities.len();
    assert!(entities >= 1, "a line was drawn");
    d.app.run(Command::CancelSketch);
    d.settle().unwrap();
    assert_eq!(d.app.mode(), Mode::PartStudio);
    assert!(feature_ids(&d).is_empty(), "cancel removes the new sketch");
}

#[test]
fn cancelling_an_edit_restores_the_sketch() {
    let (mut d, sketch, _) = demo();
    let before = d.app.model().sketch_def(sketch).unwrap().clone();
    d.app.enter_sketch(sketch);
    d.settle().unwrap();
    // Select a line by clicking it and delete it.
    let line = before
        .sketch
        .entities
        .iter()
        .find(|e| {
            matches!(
                e.geometry,
                opendatum_core_model::sketch::Geometry::Line { .. }
            )
        })
        .unwrap()
        .id;
    d.app.sketch_editor_mut().unwrap().selection = vec![line];
    d.key(Key::Delete, Modifiers::NONE).unwrap();
    assert_ne!(d.app.model().sketch_def(sketch).unwrap(), &before);
    d.app.run(Command::CancelSketch);
    d.settle().unwrap();
    assert_eq!(d.app.model().sketch_def(sketch).unwrap(), &before);
    assert_eq!(d.app.model().meshes().len(), 3);
}

#[test]
fn editing_shows_the_history_up_to_the_feature() {
    let (mut d, sketch, extrude) = demo();
    assert_eq!(d.app.model().meshes().len(), 3);
    d.app.edit_feature(extrude);
    d.frame().unwrap();
    assert_eq!(d.app.model().edit_point(), Some(extrude));
    assert_eq!(d.app.model().meshes().len(), 1, "only the edited extrude");
    assert_eq!(
        d.app.model().studio().rollback,
        InsertAt::End,
        "the bar stays"
    );
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    assert_eq!(d.app.model().edit_point(), None);
    assert_eq!(d.app.model().meshes().len(), 3);
    d.app.enter_sketch(sketch);
    d.frame().unwrap();
    assert!(
        d.app.model().meshes().is_empty(),
        "nothing before the sketch"
    );
    d.app.run(Command::FinishSketch);
    d.frame().unwrap();
    assert_eq!(d.app.model().meshes().len(), 3);
}

#[test]
fn dragging_the_rollback_bar() {
    let (mut d, sketch, extrude) = demo();
    let bar = d.app.layout().rollback_bar.unwrap();
    let above_extrude = d.app.layout().feature_row(extrude).unwrap().min.y + 2.0;
    d.drag(bar.center(), Pos2::new(bar.center().x, above_extrude), 6)
        .unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.model().studio().rollback, InsertAt::After(sketch));
    assert!(d.app.model().meshes().is_empty());
    assert!(d.app.enabled(Command::RollToEnd));
    let bar = d.app.layout().rollback_bar.unwrap();
    let bottom = d.app.layout().feature_rows.last().unwrap().1.max.y + 4.0;
    d.drag(bar.center(), Pos2::new(bar.center().x, bottom), 6)
        .unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.model().studio().rollback, InsertAt::End);
    assert_eq!(d.app.model().meshes().len(), 3);
}

#[test]
fn moving_features_treats_the_bar_as_a_row() {
    let (mut d, sketch, extrude) = demo();
    let ids = feature_ids(&d);
    let (cylinder, revolve) = (ids[2], ids[3]);
    // The last feature moves up past its neighbour.
    let op = d.app.move_feature_ops(revolve, true).unwrap();
    d.app.model_mut().apply(op).unwrap();
    assert_eq!(feature_ids(&d), vec![sketch, extrude, revolve, cylinder]);
    assert!(
        d.app.move_feature_ops(sketch, true).is_none(),
        "already first"
    );
    // With the bar after the extrude, moving the extrude down rolls it
    // back, and moving the feature below the bar up rolls it forward.
    d.app.roll_to_gap(2);
    let op = d.app.move_feature_ops(extrude, false).unwrap();
    assert!(matches!(op, ModelOp::SetRollback { at: InsertAt::After(f), .. } if f == sketch));
    d.app.model_mut().apply(op).unwrap();
    assert_eq!(feature_ids(&d), vec![sketch, extrude, revolve, cylinder]);
    assert!(d.app.model().studio().is_rolled_back(extrude));
    let op = d.app.move_feature_ops(extrude, true).unwrap();
    d.app.model_mut().apply(op).unwrap();
    assert_eq!(d.app.model().studio().rollback, InsertAt::After(extrude));
}

#[test]
fn context_menus_suppress_and_delete() {
    let (mut d, _, _) = demo();
    let cylinder = feature_ids(&d)[2];
    let row = d.app.layout().feature_row(cylinder).unwrap();
    d.right_click(row.center()).unwrap();
    d.settle().unwrap();
    let suppress = d.menu_row("Suppress").expect("the context menu is open");
    assert!(d.menu_row("Roll history to here").is_some());
    d.click(suppress.center()).unwrap();
    d.settle().unwrap();
    assert!(d.app.model().feature(cylinder).unwrap().suppressed);
    assert_eq!(d.app.model().meshes().len(), 2);
    assert!(d.menu_row("Suppress").is_none(), "the menu closed");

    let row = d.app.layout().feature_row(cylinder).unwrap();
    d.right_click(row.center()).unwrap();
    d.settle().unwrap();
    let delete = d.menu_row("Delete").unwrap();
    d.click(delete.center()).unwrap();
    d.settle().unwrap();
    assert!(d.app.model().feature(cylinder).is_none());
    assert_eq!(feature_ids(&d).len(), 3);
}

#[test]
fn escape_closes_an_open_menu_and_nothing_else() {
    let (mut d, _, _) = demo();
    let cylinder = feature_ids(&d)[2];
    let selected = Some(Selection::Feature(cylinder));
    d.app.select(selected);
    d.settle().unwrap();
    let row = d.app.layout().feature_row(cylinder).unwrap();
    d.right_click(row.center()).unwrap();
    d.settle().unwrap();
    assert!(d.menu_row("Suppress").is_some());
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert!(d.menu_row("Suppress").is_none(), "the menu closed");
    assert_eq!(d.app.selection(), selected, "the selection stays");
    // With no menu open, Escape clears the selection as before.
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.selection(), None);
}

#[test]
fn tabs_add_switch_rename_and_delete_part_studios() {
    let (mut d, _, _) = demo();
    let elements = |d: &AppDriver| d.app.model().elements().len();
    assert_eq!(elements(&d), 2);
    let add = d.app.layout().add_tab.unwrap();
    d.click(add.center()).unwrap();
    d.settle().unwrap();
    let item = d.menu_row("Part Studio").expect("the add menu");
    d.click(item.center()).unwrap();
    d.settle().unwrap();
    assert_eq!(elements(&d), 3);
    let third = d.app.model().elements()[2].id;
    assert_eq!(d.app.model().studio_id(), third, "the new tab is active");
    assert!(d.app.model().meshes().is_empty());

    let first = d.app.model().elements()[0].id;
    d.click(d.app.layout().tab(first).unwrap().center())
        .unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.model().studio_id(), first);
    assert_eq!(d.app.model().meshes().len(), 3);

    // Rename by double-click and typing.
    d.double_click(d.app.layout().tab(third).unwrap().center())
        .unwrap();
    d.settle().unwrap();
    d.key(Key::A, Modifiers::COMMAND).unwrap();
    d.type_text("Fixtures").unwrap();
    d.key(Key::Enter, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    let name = |d: &AppDriver, e| {
        d.app
            .model()
            .elements()
            .iter()
            .find(|x| x.id == e)
            .map(|x| x.name.clone())
    };
    assert_eq!(name(&d, third).as_deref(), Some("Fixtures"));

    d.right_click(d.app.layout().tab(third).unwrap().center())
        .unwrap();
    d.settle().unwrap();
    let delete = d.menu_row("Delete").unwrap();
    d.click(delete.center()).unwrap();
    d.settle().unwrap();
    assert_eq!(elements(&d), 2);
    assert!(name(&d, third).is_none());
}

#[test]
fn the_palette_runs_a_command() {
    let (mut d, _, _) = demo();
    let before = feature_ids(&d).len();
    d.key(Key::K, Modifiers::COMMAND).unwrap();
    d.settle().unwrap();
    assert!(d.app.palette_open());
    d.type_text("box").unwrap();
    d.settle().unwrap();
    d.key(Key::Enter, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert!(!d.app.palette_open());
    assert!(d.app.feature_dialog().is_some(), "the Box dialog opened");
    d.key(Key::Enter, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert_eq!(feature_ids(&d).len(), before + 1);
    let last = d.app.model().studio().features.last().unwrap();
    assert_eq!(last.name, "Box 1");
}

#[test]
fn clicking_faces_and_planes_selects_them() {
    let (mut d, _, extrude) = demo();
    d.app.select(None);
    d.settle().unwrap();
    // A point on the plate's top face, away from the hole.
    let p = on_screen(&mut d, DVec3::new(0.03, -0.015, 0.012));
    let unselected = pixel(&d.frame().unwrap(), p);
    let img = d.click(p).unwrap();
    match d.app.selection() {
        Some(Selection::Entity(hit)) => {
            assert_eq!(hit.feature, extrude);
            assert!(
                matches!(hit.pick, opendatum_naming::Pick::Face(_)),
                "{hit:?}"
            );
            assert!(
                (hit.point.z - 0.012).abs() < 1e-6,
                "on the top face: {hit:?}"
            );
        }
        other => panic!("expected the plate's top face, got {other:?}"),
    }
    // Move the pointer away, so the face shows selected, not hovered.
    let g = d.app.layout().graphics;
    d.hover(g.left_top() + Vec2::new(30.0, g.height() - 30.0))
        .unwrap();
    let selected = pixel(&d.settle().unwrap(), p);
    assert!(
        i32::from(selected[2]) - i32::from(selected[0])
            > i32::from(unselected[2]) - i32::from(unselected[0]),
        "the selected face turns towards the accent: {unselected:?} → {selected:?}"
    );
    drop(img);
    // A click on empty space clears the selection.
    let g = d.app.layout().graphics;
    d.click(g.left_top() + Vec2::new(30.0, g.height() - 30.0))
        .unwrap();
    assert_eq!(d.app.selection(), None);

    // In an empty document, a click on the planes picks one.
    let mut d = driver(Model::new());
    let p = on_screen(&mut d, DVec3::new(0.015, -0.015, 0.0));
    d.click(p).unwrap();
    assert!(
        matches!(d.app.selection(), Some(Selection::Plane(_))),
        "{:?}",
        d.app.selection()
    );
}

#[test]
fn dimensions_sit_outside_the_sketch() {
    let (model, sketch, _) = demo_model().unwrap();
    let r = model.sketch_result(sketch).unwrap();
    let frame = r.frame;
    let (yaw, pitch) = normal_view(&frame);
    let camera = Camera {
        target: DVec3::ZERO,
        distance: 0.1 / (opendatum_render::DEFAULT_FOV_Y_RAD as f64 * 0.5).tan(),
        yaw,
        pitch,
        projection: Projection::Orthographic,
        scene_bounds: None,
    };
    let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(1000.0, 800.0));
    let map = ScreenMap::new(&camera, frame, rect);
    let d = SketchDisplay::from_result(r);
    let corner = |x: f64, y: f64| {
        map.to_screen(opendatum_core_model::DVec2::new(x, y) * 1e-3)
            .unwrap()
    };
    let outline = Rect::from_two_pos(corner(-40.0, -25.0), corner(40.0, 25.0));
    let dims: Vec<_> = annotations(&d, &map)
        .into_iter()
        .filter(|a| a.dimension && a.text.contains("mm") && !a.text.contains('Ø'))
        .collect();
    assert_eq!(
        dims.len(),
        2,
        "{:?}",
        dims.iter().map(|a| &a.text).collect::<Vec<_>>()
    );
    for a in dims {
        assert!(
            !outline.intersects(a.rect),
            "{} at {:?} overlaps the rectangle {outline:?}",
            a.text,
            a.rect
        );
    }
}

#[test]
fn dragging_a_dimension_label_in_the_viewport() {
    let (mut d, sketch, _) = demo();
    d.app.enter_sketch(sketch);
    d.settle().unwrap();
    let labels = |d: &mut AppDriver| {
        let rect = d.app.layout().graphics;
        let camera = d.app.camera_mut().clone();
        let model = d.app.model();
        let frame = model.sketch_result(sketch).unwrap().frame;
        let shown = d.app.sketch_editor().unwrap().display(model).unwrap();
        annotations(&shown, &ScreenMap::new(&camera, frame, rect))
    };
    let a = labels(&mut d)
        .into_iter()
        .find(|a| a.dimension && a.text.contains("mm"))
        .expect("a length dimension");
    let to = a.rect.center() + Vec2::new(-40.0, -70.0);
    let ops = d.app.model().session().log().len();
    d.drag(a.rect.center(), to, 8).unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.model().session().log().len(), ops + 1);
    assert!(d
        .app
        .model()
        .sketch_def(sketch)
        .unwrap()
        .label(a.constraint)
        .is_some());
    let moved = labels(&mut d)
        .into_iter()
        .find(|b| b.constraint == a.constraint)
        .unwrap();
    assert!(
        moved.rect.center().distance(to) < 0.5,
        "{:?} not at {to:?}",
        moved.rect.center()
    );
    assert!(d.app.enabled(Command::ResetLabels));
    d.key(Key::Z, Modifiers::COMMAND).unwrap();
    assert_eq!(
        d.app
            .model()
            .sketch_def(sketch)
            .unwrap()
            .label(a.constraint),
        None,
        "undo returns it to automatic placement"
    );
    assert!(!d.app.enabled(Command::ResetLabels));
    // Redone, then reset from the command instead.
    d.key(Key::Z, Modifiers::COMMAND | Modifiers::SHIFT)
        .unwrap();
    assert!(d.app.enabled(Command::ResetLabels));
    d.app.run(Command::ResetLabels);
    assert!(d.app.model().sketch_def(sketch).unwrap().labels.is_empty());
}

/// The volume of the body `feature` made, as regenerated now.
fn volume(d: &AppDriver, feature: FeatureId) -> f64 {
    let m = d.app.model();
    let b = m
        .bodies()
        .bodies
        .iter()
        .find(|b| b.feature == feature)
        .expect("the feature made a body");
    m.kernel().mass_properties(&b.body, 1e-9).unwrap().volume
}

#[test]
fn a_new_feature_dialog_previews_its_values_live() {
    use opendatum_app::dialogs::Template;
    let mut d = driver(Model::new());
    let ops = d.app.model().session().log().len();
    d.app.run(Command::Template(Template::Box));
    d.frame().unwrap();
    // Shown as soon as the dialog opens, with its default values.
    assert_eq!(d.app.model().mesh_owners(), &[PREVIEW_FEATURE]);
    assert!((volume(&d, PREVIEW_FEATURE) - 40e-3 * 30e-3 * 20e-3).abs() < 1e-12);
    // One frame after a value changes, the preview has it.
    d.app.feature_dialog_mut().unwrap().fields[0].text = "80 mm".into();
    d.frame().unwrap();
    assert!((volume(&d, PREVIEW_FEATURE) - 80e-3 * 30e-3 * 20e-3).abs() < 1e-12);
    // A value that does not parse shows nothing new.
    d.app.feature_dialog_mut().unwrap().fields[1].text = "30 deg".into();
    d.frame().unwrap();
    assert!(d.app.model().meshes().is_empty());
    assert_eq!(d.app.model().session().log().len(), ops, "nothing recorded");
    // Cancelling leaves no trace.
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert!(d.app.feature_dialog().is_none());
    assert!(d.app.model().meshes().is_empty());
    assert_eq!(d.app.model().session().log().len(), ops);
    // Accepting makes the feature the preview showed.
    d.app.run(Command::Template(Template::Box));
    d.frame().unwrap();
    d.app.feature_dialog_mut().unwrap().fields[2].text = "5 mm".into();
    d.frame().unwrap();
    d.key(Key::Enter, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    let ids = feature_ids(&d);
    assert_eq!(ids.len(), 1);
    assert_eq!(d.app.model().mesh_owners(), &[ids[0]]);
    assert!(d.app.model().preview().is_none());
    assert!((volume(&d, ids[0]) - 40e-3 * 30e-3 * 5e-3).abs() < 1e-12);
}

#[test]
fn editing_an_extrude_previews_the_new_depth_until_cancelled() {
    let (mut d, _, extrude) = demo();
    let before = d.app.model().studio().clone();
    let ops = d.app.model().session().log().len();
    let v = volume(&d, extrude);
    d.app.edit_feature(extrude);
    d.frame().unwrap();
    let depth = match &d.app.model().feature(extrude).unwrap().feature {
        Feature::Extrude { distance, .. } => *distance,
        f => panic!("{f:?}"),
    };
    d.app.profile_pick_mut().unwrap().distance = format!("{}", depth * 2e3);
    d.frame().unwrap();
    assert!(
        (volume(&d, extrude) - 2.0 * v).abs() < 1e-9 * v,
        "previewed"
    );
    assert_eq!(d.app.model().studio(), &before, "the document is unchanged");
    d.key(Key::Escape, Modifiers::NONE).unwrap();
    d.settle().unwrap();
    assert!((volume(&d, extrude) - v).abs() < 1e-12);
    assert_eq!(d.app.model().session().log().len(), ops);
}

#[test]
fn the_view_cube_turns_the_view() {
    use opendatum_ui_kit::view_cube::{self, CubeFace, CubeTarget};
    let (mut d, _, _) = demo();
    let cube = d.app.layout().view_cube.unwrap();
    let camera = d.app.camera_mut().clone();
    let top = CubeTarget::Face(CubeFace::Top);
    let n = 48;
    let p = (0..n * n)
        .map(|k| {
            cube.min
                + Vec2::new(
                    cube.width() * (k % n) as f32,
                    cube.height() * (k / n) as f32,
                ) / n as f32
        })
        .find(|p| view_cube::pick(&camera, cube, *p) == Some(top))
        .expect("the top face is visible from the home view");
    d.click(p).unwrap();
    d.settle().unwrap();
    let (yaw, pitch) = top.yaw_pitch();
    let c = d.app.camera_mut();
    assert!(
        (c.pitch - pitch).abs() < 1e-9,
        "pitch {} for {pitch}",
        c.pitch
    );
    assert!((c.yaw - yaw).abs() < 1e-9, "yaw {} for {yaw}", c.yaw);
}
