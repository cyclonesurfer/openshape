//! Offscreen egui rendering for screenshots and pixel tests.
//!
//! [`UiSnapshot`] owns an egui [`egui::Context`], a headless wgpu device
//! chosen like the 3D renderer's ([`opendatum_render::BACKEND_ENV_VAR`],
//! Vulkan preferred, Mesa lavapipe in CI) and an [`egui_wgpu::Renderer`].
//! Each [`UiSnapshot::frame`] runs one egui pass, renders it into an RGBA8
//! texture and reads it back as an [`opendatum_render::RgbaImage`], so a
//! [`Color32`](egui::Color32) painted as a solid fill reads back as the same
//! sRGB bytes.
//!
//! # Target format
//!
//! The target is `Rgba8Unorm` holding sRGB-encoded bytes, which is what the
//! desktop app's window surface uses: egui-wgpu prefers a non-sRGB surface
//! format and then blends in gamma space, which is how egui's premultiplied
//! colours are meant to be composited. An `Rgba8UnormSrgb` target would
//! make the GPU blend in linear space instead, so every translucent shape
//! (disabled widgets, shadows, antialiased edges and text) would come out
//! darker than in the app: for example a disabled button's white fill at
//! 50 % opacity over white reads back as grey 220 instead of 255.
//! [`UiSnapshot::with_format`] can still select `Rgba8UnormSrgb` (for
//! instance to compare the two); solid colours read back the same in both.
//!
//! The renderer uses egui's predictable options (no MSAA, no dithering,
//! texture filtering in the shader) so the same input gives the same pixels
//! on every GPU as far as the hardware allows.

use egui::mutex::RwLock;
use opendatum_render::offscreen::HeadlessGpu;
use opendatum_render::{read_texture_rgba8, RenderError, RgbaImage};
use std::sync::Arc;

/// Renders egui frames offscreen and reads them back.
pub struct UiSnapshot {
    ctx: egui::Context,
    gpu: HeadlessGpu,
    renderer: Arc<RwLock<egui_wgpu::Renderer>>,
    size: [u32; 2],
    pixels_per_point: f32,
    target: wgpu::Texture,
    view: wgpu::TextureView,
}

impl UiSnapshot {
    /// Default colour format of the render target: sRGB-encoded bytes in a
    /// non-sRGB format, like the app's window surface (see the
    /// [module documentation](self#target-format)).
    pub const FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;

    /// Creates a snapshot renderer for `size` physical pixels at
    /// `pixels_per_point` (the UI's logical size is `size /
    /// pixels_per_point` points), rendering into [`UiSnapshot::FORMAT`].
    ///
    /// Fails if no GPU adapter is available or the size is zero or larger
    /// than the device allows. Panics if `pixels_per_point` is not positive
    /// and finite.
    pub fn new(size: [u32; 2], pixels_per_point: f32) -> Result<Self, RenderError> {
        Self::with_format(size, pixels_per_point, Self::FORMAT)
    }

    /// Like [`UiSnapshot::new`] with an explicit target format:
    /// `Rgba8Unorm` (the default) or `Rgba8UnormSrgb`. Any other format is
    /// a [`RenderError::UnsupportedTexture`] error.
    pub fn with_format(
        size: [u32; 2],
        pixels_per_point: f32,
        format: wgpu::TextureFormat,
    ) -> Result<Self, RenderError> {
        if !matches!(
            format,
            wgpu::TextureFormat::Rgba8Unorm | wgpu::TextureFormat::Rgba8UnormSrgb
        ) {
            return Err(RenderError::UnsupportedTexture(format!(
                "snapshot target format {format:?} (expected Rgba8Unorm or Rgba8UnormSrgb)"
            )));
        }
        assert!(
            pixels_per_point.is_finite() && pixels_per_point > 0.0,
            "pixels_per_point must be positive and finite, got {pixels_per_point}"
        );
        let gpu = HeadlessGpu::new("opendatum-ui-kit snapshot")?;
        let max = gpu.device.limits().max_texture_dimension_2d;
        let [width, height] = size;
        if width == 0 || height == 0 || width > max || height > max {
            return Err(RenderError::InvalidSize { width, height, max });
        }
        let renderer =
            egui_wgpu::Renderer::new(&gpu.device, format, egui_wgpu::RendererOptions::PREDICTABLE);
        let target = gpu.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("opendatum-ui-kit snapshot target"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = target.create_view(&wgpu::TextureViewDescriptor::default());
        let info = gpu.adapter.get_info();
        log::info!(
            "opendatum-ui-kit: snapshot on {} ({:?})",
            info.name,
            info.backend
        );
        Ok(UiSnapshot {
            ctx: egui::Context::default(),
            gpu,
            renderer: Arc::new(RwLock::new(renderer)),
            size,
            pixels_per_point,
            target,
            view,
        })
    }

    /// The egui context the frames run in (apply a theme to it, keep state
    /// across frames).
    pub fn context(&self) -> &egui::Context {
        &self.ctx
    }

    /// Image size in physical pixels.
    pub fn size(&self) -> [u32; 2] {
        self.size
    }

    /// Colour format of the render target.
    pub fn format(&self) -> wgpu::TextureFormat {
        self.target.format()
    }

    /// Physical pixels per egui point.
    pub fn pixels_per_point(&self) -> f32 {
        self.pixels_per_point
    }

    /// The UI's area in points (what egui's screen rect will be).
    pub fn screen_rect(&self) -> egui::Rect {
        egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(self.size[0] as f32, self.size[1] as f32) / self.pixels_per_point,
        )
    }

    /// Information about the adapter in use.
    pub fn adapter_info(&self) -> wgpu::AdapterInfo {
        self.gpu.adapter.get_info()
    }

    /// The device, queue and egui renderer as an [`egui_wgpu::RenderState`],
    /// so an app can register native wgpu textures (for example its 3D
    /// viewport image) through the same renderer that draws the snapshots.
    pub fn render_state(&self) -> egui_wgpu::RenderState {
        egui_wgpu::RenderState {
            adapter: self.gpu.adapter.clone(),
            available_adapters: vec![self.gpu.adapter.clone()],
            instance: self.gpu.instance.clone(),
            device: self.gpu.device.clone(),
            queue: self.gpu.queue.clone(),
            target_format: self.target.format(),
            renderer: Arc::clone(&self.renderer),
            surface_config: egui_wgpu::SurfaceConfig::LOW_LATENCY,
        }
    }

    /// Runs one egui pass with `input` (screen size and pixels per point are
    /// filled in when missing), calling `run` with the context (use
    /// [`egui::Area`], [`egui::Window`] or a root [`egui::Ui`] built from
    /// it), then renders and reads back the frame.
    pub fn frame(
        &mut self,
        input: egui::RawInput,
        mut run: impl FnMut(&egui::Context),
    ) -> Result<RgbaImage, RenderError> {
        let input = self.prepare_input(input);
        let output = self.ctx.run_ui(input, |ui| run(ui.ctx()));
        self.render(output)
    }

    /// Like [`UiSnapshot::frame`], but `run` gets the root [`egui::Ui`]
    /// covering the whole screen (for panel-based layouts, as in an
    /// application's `ui` method).
    pub fn frame_ui(
        &mut self,
        input: egui::RawInput,
        run: impl FnMut(&mut egui::Ui),
    ) -> Result<RgbaImage, RenderError> {
        let input = self.prepare_input(input);
        let output = self.ctx.run_ui(input, run);
        self.render(output)
    }

    fn prepare_input(&self, mut input: egui::RawInput) -> egui::RawInput {
        if input.screen_rect.is_none() {
            input.screen_rect = Some(self.screen_rect());
        }
        if input.max_texture_side.is_none() {
            input.max_texture_side =
                Some(self.gpu.device.limits().max_texture_dimension_2d as usize);
        }
        let viewport = input.viewports.entry(input.viewport_id).or_default();
        if viewport.native_pixels_per_point.is_none() {
            viewport.native_pixels_per_point = Some(self.pixels_per_point);
        }
        input
    }

    fn render(&mut self, mut output: egui::FullOutput) -> Result<RgbaImage, RenderError> {
        // Applied (set before painting, freed after) below; `TexturesDelta`
        // must be emptied before it is dropped.
        let mut textures = std::mem::take(&mut output.textures_delta);
        let device = &self.gpu.device;
        let queue = &self.gpu.queue;
        let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
        let screen = egui_wgpu::ScreenDescriptor {
            size_in_pixels: self.size,
            pixels_per_point: output.pixels_per_point,
        };
        let primitives = self
            .ctx
            .tessellate(std::mem::take(&mut output.shapes), output.pixels_per_point);
        let clear = clear_color(
            self.ctx.global_style().visuals.panel_fill,
            self.target.format().is_srgb(),
        );

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("opendatum-ui-kit snapshot frame"),
        });
        let user_buffers = {
            let mut renderer = self.renderer.write();
            for (id, deltas) in textures.set.drain() {
                for delta in &deltas {
                    renderer.update_texture(device, queue, id, delta);
                }
            }
            renderer.update_buffers(device, queue, &mut encoder, &primitives, &screen)
        };
        {
            let renderer = self.renderer.read();
            let mut pass = encoder
                .begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("opendatum-ui-kit snapshot pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: &self.view,
                        depth_slice: None,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(clear),
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: None,
                    timestamp_writes: None,
                    occlusion_query_set: None,
                    multiview_mask: None,
                })
                .forget_lifetime();
            renderer.render(&mut pass, &primitives, &screen);
        }
        queue.submit(user_buffers.into_iter().chain(Some(encoder.finish())));
        {
            let mut renderer = self.renderer.write();
            for id in textures.free.drain() {
                renderer.free_texture(&id);
            }
        }
        if let Some(e) = pollster::block_on(scope.pop()) {
            return Err(RenderError::Gpu(e.to_string()));
        }
        read_texture_rgba8(device, queue, &self.target)
    }
}

/// The clear value that stores `color` (opaque) in a target: linear values
/// for an sRGB format (the GPU encodes them), the sRGB bytes themselves
/// otherwise.
fn clear_color(color: egui::Color32, srgb_target: bool) -> wgpu::Color {
    if srgb_target {
        let linear = egui::Rgba::from(color.to_opaque());
        wgpu::Color {
            r: linear.r() as f64,
            g: linear.g() as f64,
            b: linear.b() as f64,
            a: 1.0,
        }
    } else {
        let byte = |b: u8| b as f64 / 255.0;
        wgpu::Color {
            r: byte(color.r()),
            g: byte(color.g()),
            b: byte(color.b()),
            a: 1.0,
        }
    }
}
