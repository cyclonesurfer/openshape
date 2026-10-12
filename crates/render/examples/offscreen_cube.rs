//! Renders a shaded, edge-outlined unit cube headlessly and writes a PNG.
//!
//! ```sh
//! cargo run -p opendatum-render --example offscreen_cube -- out/cube.png [width height]
//! ```
//!
//! The backend can be forced with `OPENDATUM_WGPU_BACKEND=vulkan|gl|any`.

use opendatum_render::{write_png, Camera, MeshData, OffscreenRenderer};
use std::process::ExitCode;

fn main() -> ExitCode {
    env_logger::init();
    let mut args = std::env::args().skip(1);
    let path = args.next().unwrap_or_else(|| "cube.png".to_string());
    let width = args.next().and_then(|s| s.parse().ok()).unwrap_or(640);
    let height = args.next().and_then(|s| s.parse().ok()).unwrap_or(480);

    match run(&path, width, height) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("offscreen_cube: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(path: &str, width: u32, height: u32) -> Result<(), opendatum_render::RenderError> {
    let cube = MeshData::cuboid([-0.5; 3], [0.5; 3], MeshData::DEFAULT_COLOR);
    let (min, max) = cube.bounds().expect("cube has bounds");
    let mut camera = Camera::isometric();
    camera.fit_with_aspect(min, max, width as f32 / height as f32);

    let mut renderer = OffscreenRenderer::new()?;
    let info = renderer.adapter_info();
    println!(
        "adapter: {} ({:?}, {:?}, {} {})",
        info.name, info.backend, info.device_type, info.driver, info.driver_info
    );
    let image = renderer.render_meshes_rgba(&[cube], &camera, width, height)?;
    write_png(path, &image)?;
    println!("wrote {path} ({width}x{height})");
    Ok(())
}
