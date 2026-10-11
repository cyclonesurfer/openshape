//! Shared helpers for the GPU tests (they run on Mesa lavapipe in CI).
#![allow(dead_code)]

use openshape_render::{write_png, RgbaImage};
use openshape_ui_kit::UiSnapshot;

/// Creates a snapshot renderer or panics loudly: these tests must never pass
/// silently on a machine without a usable adapter.
pub fn snapshot(size: [u32; 2], pixels_per_point: f32) -> UiSnapshot {
    match UiSnapshot::new(size, pixels_per_point) {
        Ok(s) => {
            let info = s.adapter_info();
            eprintln!("adapter: {} ({:?})", info.name, info.backend);
            s
        }
        Err(e) => panic!(
            "no wgpu adapter available ({e}); install Mesa lavapipe (mesa-vulkan-drivers) \
             or set OPENSHAPE_WGPU_BACKEND"
        ),
    }
}

/// Writes `image` to `$OPENSHAPE_RENDER_TEST_OUT/<name>.png` when that
/// variable is set (for eyeballing test renders).
pub fn maybe_dump(name: &str, image: &RgbaImage) {
    if let Ok(dir) = std::env::var("OPENSHAPE_RENDER_TEST_OUT") {
        std::fs::create_dir_all(&dir).expect("create dump dir");
        let path = std::path::Path::new(&dir).join(format!("{name}.png"));
        write_png(&path, image).expect("write png");
    }
}

/// Largest per-channel RGB difference.
pub fn max_channel_diff(a: [u8; 4], b: [u8; 4]) -> u8 {
    (0..3).map(|i| a[i].abs_diff(b[i])).max().unwrap()
}

/// The RGBA bytes of an opaque colour.
pub fn rgba(c: egui::Color32) -> [u8; 4] {
    [c.r(), c.g(), c.b(), 255]
}
