//! Picking faces, edges and vertices in the view, sketching on a picked
//! face and extruding up to one (ADR 0019), end to end in the real app
//! offscreen. They need a Vulkan or GL adapter (lavapipe in CI).

use egui::{Pos2, Vec2};
use opendatum_app::app::Selection;
use opendatum_app::commands::Command;
use opendatum_app::entity_pick::EntityHit;
use opendatum_app::profile_pick::EndCondition;
use opendatum_app::screenshots::{demo_model, AppDriver};
use opendatum_app::sketch_mode::ScreenMap;
use opendatum_app::Model;
use opendatum_core_model::naming::{FaceName, Name, Role};
use opendatum_core_model::{DVec3, Feature, FeatureId, PlaneRef};
use opendatum_kernel::{Frame, Kernel};
use opendatum_naming::Pick;
use opendatum_render::RgbaImage;

const SIZE: [u32; 2] = [1280, 800];

fn demo() -> (AppDriver, FeatureId, FeatureId) {
    let (model, sketch, extrude) = demo_model().expect("demo model builds");
    let d = AppDriver::new(SIZE, 1.0, model).expect("an adapter (Vulkan or GL)");
    (d, sketch, extrude)
}

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

/// How far a pixel leans towards blue (the accent and preselection hues).
fn blueness(p: [u8; 4]) -> i32 {
    i32::from(p[2]) - i32::from(p[0])
}

/// The name of the entity a hit picked.
fn name(d: &AppDriver, hit: &EntityHit) -> Name {
    d.app
        .model()
        .part(hit.feature)
        .unwrap()
        .names
        .as_ref()
        .unwrap()
        .name(hit.pick)
        .expect("a named entity")
}

fn top_cap(n: &FaceName) -> bool {
    n.role == Role::Cap { end: true }
}

/// The cylinder template's part (the demo's third feature).
fn cylinder(d: &AppDriver) -> FeatureId {
    d.app.model().studio().features[2].id
}

#[test]
fn the_pointer_picks_vertices_then_edges_then_faces() {
    let (mut d, _, extrude) = demo();
    d.app.select(None);
    let away = {
        let g = d.app.layout().graphics;
        g.left_top() + Vec2::new(30.0, g.height() - 30.0)
    };
    d.hover(away).unwrap();
    let face_point = on_screen(&mut d, DVec3::new(0.03, -0.015, 0.012));
    let plain = pixel(&d.settle().unwrap(), face_point);

    // Inside the plate's top face.
    d.hover(face_point).unwrap();
    let lit = pixel(&d.settle().unwrap(), face_point);
    let hit = d
        .app
        .hovered_entity()
        .expect("the top face is under the pointer");
    assert_eq!(hit.feature, extrude);
    match name(&d, &hit) {
        Name::Face(f) => assert!(top_cap(&f), "{f:?}"),
        other => panic!("expected a face, got {other:?}"),
    }
    assert!(
        blueness(lit) > blueness(plain),
        "the hovered face lights up: {plain:?} → {lit:?}"
    );

    // Two points from the middle of its front edge.
    let edge = on_screen(&mut d, DVec3::new(0.0, -0.025, 0.012)) + Vec2::new(0.0, 2.0);
    d.hover(edge).unwrap();
    let hit = d.app.hovered_entity().expect("an entity under the pointer");
    assert!(matches!(hit.pick, Pick::Edge(_)), "{hit:?}");
    match name(&d, &hit) {
        Name::Edge(e) => {
            assert_eq!(e.faces.len(), 2);
            assert!(e.faces.iter().any(top_cap), "{e:?}");
        }
        other => panic!("expected an edge, got {other:?}"),
    }

    // Near a corner the vertex wins over its edges and faces.
    let corner = on_screen(&mut d, DVec3::new(0.04, -0.025, 0.012)) + Vec2::new(-3.0, 1.0);
    d.hover(corner).unwrap();
    let hit = d.app.hovered_entity().expect("an entity under the pointer");
    assert!(matches!(hit.pick, Pick::Vertex(_)), "{hit:?}");
    // The demo's corners were placed by clicking: within a micron.
    let off = (hit.point - DVec3::new(0.04, -0.025, 0.012)).length();
    assert!(off < 1e-6, "the corner is {off} away");
    match name(&d, &hit) {
        Name::Vertex(v) => assert_eq!(v.faces.len(), 3, "{v:?}"),
        other => panic!("expected a vertex, got {other:?}"),
    }

    // Clicking selects what is under the pointer; nothing is hovered off
    // the parts.
    d.click(corner).unwrap();
    assert!(matches!(
        d.app.selection(),
        Some(Selection::Entity(EntityHit {
            pick: Pick::Vertex(_),
            ..
        }))
    ));
    d.hover(away).unwrap();
    assert_eq!(d.app.hovered_entity(), None);
}

#[test]
fn a_sketch_starts_on_a_picked_planar_face() {
    let (mut d, _, extrude) = demo();
    let p = on_screen(&mut d, DVec3::new(0.03, -0.015, 0.012));
    d.click(p).unwrap();
    let before = d.app.model().studio().features.len();
    d.app.run(Command::Sketch);
    d.settle().unwrap();
    let ed = d.app.sketch_editor().expect("sketch mode");
    let sketch = ed.feature;
    assert_eq!(d.app.model().studio().features.len(), before + 1);
    match &d.app.model().sketch_def(sketch).unwrap().plane {
        PlaneRef::Face(r) => match &r.name {
            Name::Face(f) => {
                assert_eq!(f.feature, extrude);
                assert!(top_cap(f));
                assert!(
                    (r.hint.point.z - 0.012).abs() < 1e-12,
                    "the hint is on the face"
                );
            }
            other => panic!("expected a face, got {other:?}"),
        },
        other => panic!("expected a face plane, got {other:?}"),
    }
    let frame = d.app.model().sketch_result(sketch).unwrap().frame;
    assert!((frame.origin - DVec3::new(0.0, 0.0, 0.012)).length() < 1e-12);
    assert!((frame.z - DVec3::Z).length() < 1e-12);
    // Cancelling removes the sketch again.
    d.app.run(Command::CancelSketch);
    d.settle().unwrap();
    assert_eq!(d.app.model().studio().features.len(), before);

    // A curved face (the cylinder's side) is refused.
    let c = cylinder(&d);
    let side = on_screen(&mut d, DVec3::new(0.06, 0.005, 0.015));
    d.click(side).unwrap();
    match d.app.selection() {
        Some(Selection::Entity(hit)) => assert_eq!(hit.feature, c, "{hit:?}"),
        other => panic!("expected the cylinder's side, got {other:?}"),
    }
    d.app.run(Command::Sketch);
    d.settle().unwrap();
    assert!(d.app.sketch_editor().is_none());
    assert_eq!(d.app.model().studio().features.len(), before);
}

#[test]
fn an_extrude_ends_at_a_picked_face() {
    let (mut d, sketch, _) = demo();
    let c = cylinder(&d);
    d.app.run(Command::Extrude);
    d.settle().unwrap();
    // The plate's region, on the sketch plane under the plate.
    let region = on_screen(&mut d, DVec3::new(0.03, -0.015, 0.0));
    d.click(region).unwrap();
    {
        let p = d.app.profile_pick_mut().expect("the extrude dialog");
        assert_eq!(p.sketch, Some(sketch));
        p.set_end(EndCondition::UpToFace);
        assert!(p.picking_face);
    }
    // A curved face is refused with a message; the top of the cylinder
    // is taken.
    let side = on_screen(&mut d, DVec3::new(0.06, 0.005, 0.015));
    d.click(side).unwrap();
    {
        let p = d.app.profile_pick_mut().unwrap();
        assert_eq!(p.error.as_deref(), Some("Pick a planar face"));
        assert!(p.up_to.is_none() && p.picking_face);
    }
    let top = on_screen(&mut d, DVec3::new(0.06, 0.015, 0.03));
    d.click(top).unwrap();
    {
        let p = d.app.profile_pick_mut().unwrap();
        assert!(!p.picking_face && p.error.is_none());
        match &p.up_to.as_ref().expect("a face").name {
            Name::Face(f) => {
                assert_eq!(f.feature, c);
                assert!(top_cap(f));
            }
            other => panic!("expected a face, got {other:?}"),
        }
    }
    d.key(egui::Key::Enter, egui::Modifiers::NONE).unwrap();
    d.settle().unwrap();
    let last = d.app.model().studio().features.last().unwrap().clone();
    match &last.feature {
        Feature::Extrude { up_to, .. } => assert!(up_to.is_some()),
        other => panic!("expected an extrude, got {other:?}"),
    }
    // The plate's profile (80 × 50 mm less the 24 mm hole) up to the
    // cylinder's top, 30 mm up.
    let m = d.app.model();
    let body = &m.part(last.id).expect("the extrude regenerates").body;
    let volume = m.kernel().mass_properties(body, 1e-10).unwrap().volume;
    let area = 0.08 * 0.05 - std::f64::consts::PI * 0.012f64.powi(2);
    assert!((volume - area * 0.03).abs() < 1e-9 * volume, "{volume}");
    let bounds = d.app.model().body_mesh(last.id).unwrap().bounds().unwrap();
    assert!((bounds.1.z - 0.03).abs() < 1e-9);
}

#[test]
fn parts_without_entities_under_the_pointer_pick_nothing() {
    let mut d = AppDriver::new(SIZE, 1.0, Model::new()).expect("an adapter (Vulkan or GL)");
    d.settle().unwrap();
    let g = d.app.layout().graphics;
    d.hover(g.center()).unwrap();
    assert_eq!(d.app.hovered_entity(), None);
}
