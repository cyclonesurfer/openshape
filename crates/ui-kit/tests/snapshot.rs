//! `UiSnapshot` renders egui offscreen with exact colours, egui's intended
//! (gamma-space) blending, correct pixel scaling, text, and native textures
//! registered through its render state. Runs on the GPU (lavapipe in CI).

mod common;

use common::{max_channel_diff, maybe_dump, rgba, snapshot};
use egui::{Color32, LayerId, Pos2, Rect, Vec2};
use opendatum_render::RenderError;
use opendatum_ui_kit::{Theme, UiSnapshot};

/// Paints `cells` side by side as 16x16-point squares on `background`.
fn paint_cells(ctx: &egui::Context, background: Color32, cells: &[Color32]) {
    let painter = ctx.layer_painter(LayerId::background());
    painter.rect_filled(ctx.content_rect(), 0.0, background);
    for (i, c) in cells.iter().enumerate() {
        let min = Pos2::new((i % 16) as f32 * 16.0, (i / 16) as f32 * 16.0);
        painter.rect_filled(Rect::from_min_size(min, Vec2::splat(16.0)), 0.0, *c);
    }
}

fn palette() -> Vec<Color32> {
    let t = Theme::light();
    let mut colors = vec![
        t.chrome,
        t.panel,
        t.panel_header,
        t.border,
        t.border_strong,
        t.text,
        t.text_muted,
        t.text_disabled,
        t.accent,
        t.accent_hover,
        t.accent_soft,
        t.on_accent,
        t.hover_soft,
        t.selection,
        t.preselect,
        t.success,
        t.warning,
        t.danger,
        t.viewport_top,
        t.viewport_bottom,
        Color32::BLACK,
        Color32::WHITE,
    ];
    // A spread of other values, including the low end where sRGB is linear.
    for i in 0..10u8 {
        colors.push(Color32::from_rgb(
            i.wrapping_mul(29),
            1 + i.wrapping_mul(97),
            255 - i.wrapping_mul(13),
        ));
    }
    colors
}

fn check_round_trip(mut snap: UiSnapshot, name: &str) {
    let colors = palette();
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            paint_cells(ctx, Color32::from_gray(77), &colors)
        })
        .unwrap();
    maybe_dump(name, &image);
    assert_eq!((image.width, image.height), (256, 48));
    for (i, c) in colors.iter().enumerate() {
        let (x, y) = ((i % 16) as u32 * 16 + 8, (i / 16) as u32 * 16 + 8);
        let px = image.pixel(x, y);
        assert!(
            max_channel_diff(px, rgba(*c)) <= 1 && px[3] == 255,
            "{name}: colour {i} {c:?} read back as {px:?}"
        );
    }
}

#[test]
fn solid_colours_round_trip_as_srgb_bytes() {
    check_round_trip(snapshot([256, 48], 1.0), "snapshot_colours");
}

#[test]
fn solid_colours_round_trip_on_an_srgb_target_too() {
    let snap = UiSnapshot::with_format([256, 48], 1.0, wgpu::TextureFormat::Rgba8UnormSrgb)
        .expect("adapter");
    assert_eq!(snap.format(), wgpu::TextureFormat::Rgba8UnormSrgb);
    check_round_trip(snap, "snapshot_colours_srgb");
}

/// Half-transparent shapes over white, as egui composites them in the app's
/// (non-sRGB) window surface.
fn translucent(snap: &mut UiSnapshot) -> (opendatum_render::RgbaImage, Vec<[u8; 4]>) {
    let cells = [
        Color32::WHITE.gamma_multiply(0.5),
        Color32::BLACK.gamma_multiply(0.5),
        Theme::light().preselect.gamma_multiply(0.6),
    ];
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            paint_cells(ctx, Color32::WHITE, &cells)
        })
        .unwrap();
    let px = (0..cells.len())
        .map(|i| image.pixel(i as u32 * 16 + 8, 8))
        .collect();
    (image, px)
}

#[test]
fn translucent_colours_blend_like_the_app_window() {
    let mut snap = snapshot([256, 16], 1.0);
    let (image, px) = translucent(&mut snap);
    maybe_dump("snapshot_translucent", &image);
    // Premultiplied blending in gamma space: src + (1 - a) * dst.
    assert!(
        max_channel_diff(px[0], [255, 255, 255, 255]) <= 1,
        "{:?}",
        px[0]
    );
    assert!(
        max_channel_diff(px[1], [128, 128, 128, 255]) <= 1,
        "{:?}",
        px[1]
    );
    let p = Theme::light().preselect;
    let over_white = |c: u8| ((c as f32 * 0.6) + 0.4 * 255.0).round() as u8;
    let expected = [over_white(p.r()), over_white(p.g()), over_white(p.b()), 255];
    assert!(
        max_channel_diff(px[2], expected) <= 2,
        "{:?} vs {expected:?}",
        px[2]
    );
}

#[test]
fn an_srgb_target_would_darken_translucent_shapes() {
    // Why the default target is not sRGB: the GPU then blends egui's
    // premultiplied gamma colours in linear space.
    let mut snap =
        UiSnapshot::with_format([256, 16], 1.0, wgpu::TextureFormat::Rgba8UnormSrgb).unwrap();
    let (_, px) = translucent(&mut snap);
    assert!(px[0][0] < 235, "white at 50% over white: {:?}", px[0]);
    assert!(px[1][0] > 170, "black at 50% over white: {:?}", px[1]);
}

#[test]
fn pixels_per_point_scales_the_ui() {
    let mut snap = snapshot([80, 60], 2.0);
    assert_eq!(snap.screen_rect().size(), Vec2::new(40.0, 30.0));
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            assert_eq!(ctx.pixels_per_point(), 2.0);
            let painter = ctx.layer_painter(LayerId::background());
            painter.rect_filled(ctx.content_rect(), 0.0, Color32::WHITE);
            painter.rect_filled(
                Rect::from_min_max(Pos2::new(10.0, 10.0), Pos2::new(20.0, 20.0)),
                0.0,
                Color32::RED,
            );
        })
        .unwrap();
    let red = [255, 0, 0, 255];
    let white = [255, 255, 255, 255];
    assert_eq!(image.pixel(20, 20), red);
    assert_eq!(image.pixel(39, 39), red);
    assert_eq!(image.pixel(19, 19), white);
    assert_eq!(image.pixel(40, 40), white);
}

#[test]
fn text_renders_across_frames() {
    let mut snap = snapshot([160, 40], 1.0);
    let theme = Theme::light();
    theme.apply(snap.context());
    let mut dark_pixels = Vec::new();
    for _ in 0..3 {
        let image = snap
            .frame_ui(egui::RawInput::default(), |ui| {
                ui.painter().rect_filled(ui.max_rect(), 0.0, theme.panel);
                ui.label("Extrude 1 \u{2016} \u{27C2}");
            })
            .unwrap();
        let dark = image
            .pixels
            .chunks_exact(4)
            .filter(|p| p[0] < 128 && p[1] < 128)
            .count();
        dark_pixels.push(dark);
    }
    maybe_dump("snapshot_text", &{
        snap.frame_ui(egui::RawInput::default(), |ui| {
            ui.label("Extrude 1 \u{2016} \u{27C2}");
        })
        .unwrap()
    });
    assert!(dark_pixels[0] > 30, "no text drawn: {dark_pixels:?}");
    assert!(
        dark_pixels.windows(2).all(|w| w[0] == w[1]),
        "frames differ: {dark_pixels:?}"
    );
}

#[test]
fn native_textures_register_through_the_render_state() {
    let mut snap = snapshot([64, 64], 1.0);
    let rs = snap.render_state();
    assert_eq!(rs.target_format, UiSnapshot::FORMAT);
    let color = [12u8, 150, 222, 255];
    let texture = rs.device.create_texture(&wgpu::TextureDescriptor {
        label: Some("test native texture"),
        size: wgpu::Extent3d {
            width: 4,
            height: 4,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Rgba8Unorm,
        usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
        view_formats: &[],
    });
    let data: Vec<u8> = color.repeat(16);
    rs.queue.write_texture(
        texture.as_image_copy(),
        &data,
        wgpu::TexelCopyBufferLayout {
            offset: 0,
            bytes_per_row: Some(16),
            rows_per_image: Some(4),
        },
        texture.size(),
    );
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let id =
        rs.renderer
            .write()
            .register_native_texture(&rs.device, &view, wgpu::FilterMode::Nearest);
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            ctx.layer_painter(LayerId::background()).image(
                id,
                Rect::from_min_size(Pos2::new(16.0, 16.0), Vec2::splat(32.0)),
                Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
        })
        .unwrap();
    assert!(max_channel_diff(image.pixel(32, 32), color) <= 1);
    // Outside the image: the panel fill the frame is cleared to.
    let panel = rgba(snap.context().global_style().visuals.panel_fill);
    assert_eq!(image.pixel(4, 4), panel);
}

#[test]
fn invalid_sizes_and_formats_are_errors() {
    assert!(matches!(
        UiSnapshot::new([0, 10], 1.0),
        Err(RenderError::InvalidSize { .. })
    ));
    assert!(matches!(
        UiSnapshot::new([10, 1 << 20], 1.0),
        Err(RenderError::InvalidSize { .. })
    ));
    assert!(matches!(
        UiSnapshot::with_format([10, 10], 1.0, wgpu::TextureFormat::Bgra8Unorm),
        Err(RenderError::UnsupportedTexture(_))
    ));
}
