//! The 3D view shown through egui has the same colours as the scene
//! renderer's own output (crates/ui/src/viewport.rs): egui samples native
//! textures as gamma-encoded bytes, so it must not read the sRGB view.

use openshape_app::viewport::{scene_format, ViewportTexture};
use openshape_render::{read_texture_rgba8, wgpu, Camera, MeshData, SceneRenderer};
use openshape_ui_kit::snapshot::UiSnapshot;

const W: u32 = 160;
const H: u32 = 120;

fn scene() -> (UiSnapshot, ViewportTexture, Vec<u8>) {
    let snap = UiSnapshot::new([W, H], 1.0).expect("a wgpu adapter (install Mesa lavapipe)");
    let rs = snap.render_state();
    let format = scene_format(&rs.adapter);
    let mut scene = SceneRenderer::new(&rs.device, format);
    let parts = [
        MeshData::cuboid([0.0, 0.0, 0.0], [1.0, 0.6, 0.4], [0.55, 0.62, 0.72, 1.0]),
        MeshData::cuboid([0.2, 0.1, 0.4], [0.6, 0.4, 0.9], [0.85, 0.45, 0.2, 1.0]),
    ];
    scene.set_meshes(&rs.device, &parts);
    let mut camera = Camera::isometric();
    camera.fit([0.0, 0.0, 0.0], [1.0, 0.6, 0.9]);
    let t = ViewportTexture::new(&rs.device, format, (W, H));
    scene.render(&rs.device, &rs.queue, &camera, &t.render_view, (W, H));
    let direct = read_texture_rgba8(&rs.device, &rs.queue, &t.texture)
        .unwrap()
        .pixels;
    (snap, t, direct)
}

fn shown(snap: &mut UiSnapshot, view: &wgpu::TextureView) -> Vec<u8> {
    let rs = snap.render_state();
    let id =
        rs.renderer
            .write()
            .register_native_texture(&rs.device, view, wgpu::FilterMode::Nearest);
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            ctx.layer_painter(egui::LayerId::background()).image(
                id,
                egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(W as f32, H as f32)),
                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                egui::Color32::WHITE,
            );
        })
        .unwrap();
    image.pixels
}

fn worst(a: &[u8], b: &[u8]) -> u8 {
    assert_eq!(a.len(), b.len());
    a.chunks(4)
        .zip(b.chunks(4))
        .flat_map(|(p, q)| (0..3).map(move |c| p[c].abs_diff(q[c])))
        .max()
        .unwrap()
}

#[test]
fn the_view_shown_by_egui_matches_the_rendered_scene() {
    let (mut snap, t, direct) = scene();
    // The scene is not blank: several distinct colours.
    let mut colours: Vec<[u8; 3]> = direct.chunks(4).map(|p| [p[0], p[1], p[2]]).collect();
    colours.sort();
    colours.dedup();
    assert!(colours.len() > 20, "{} colours", colours.len());
    let egui = shown(&mut snap, &t.egui_view);
    let w = worst(&direct, &egui);
    assert!(w <= 1, "egui shows the view up to {w} codes off");
}

#[test]
fn reading_the_srgb_view_would_darken_it() {
    // The defect the separate egui view prevents: through the sRGB view the
    // shader receives linear values and treats them as gamma-encoded.
    let (mut snap, t, direct) = scene();
    if t.texture.format() != wgpu::TextureFormat::Rgba8UnormSrgb {
        return; // GL: the texture is not sRGB, so there is no such view.
    }
    let wrong = shown(&mut snap, &t.render_view);
    assert!(worst(&direct, &wrong) > 40);
}
