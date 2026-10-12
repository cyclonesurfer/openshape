//! Headless rendering: [`OffscreenRenderer`] owns its own wgpu instance,
//! adapter and device and reads rendered frames back as [`RgbaImage`]s.

use crate::camera::Camera;
use crate::mesh::MeshData;
use crate::scene::SceneRenderer;
use crate::RenderError;
use std::path::Path;

/// Environment variable selecting the backend of [`OffscreenRenderer::new`]
/// and of the application window: `vulkan`, `gl` or `any`
/// (case-insensitive). Unset means [`BackendPreference::PreferVulkan`].
pub const BACKEND_ENV_VAR: &str = "OPENSHAPE_WGPU_BACKEND";

/// Which wgpu backend an [`OffscreenRenderer`] may use.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum BackendPreference {
    /// Vulkan if it has any adapter (hardware or lavapipe), otherwise the
    /// other backends (the default).
    #[default]
    PreferVulkan,
    /// Vulkan only (`OPENSHAPE_WGPU_BACKEND=vulkan`).
    Vulkan,
    /// OpenGL / GLES only (`OPENSHAPE_WGPU_BACKEND=gl`).
    Gl,
    /// Whatever wgpu picks among all backends (`OPENSHAPE_WGPU_BACKEND=any`).
    Any,
}

impl BackendPreference {
    /// Parses `vulkan` / `vk`, `gl` / `gles` / `opengl`, or `any` (case-insensitive).
    pub fn parse(value: &str) -> Option<BackendPreference> {
        match value.trim().to_ascii_lowercase().as_str() {
            "vulkan" | "vk" => Some(BackendPreference::Vulkan),
            "gl" | "gles" | "opengl" => Some(BackendPreference::Gl),
            "any" | "all" => Some(BackendPreference::Any),
            _ => None,
        }
    }

    /// Reads [`BACKEND_ENV_VAR`]; unset or empty gives the default, an
    /// unrecognised value is an error (so typos fail loudly).
    pub fn from_env() -> Result<BackendPreference, RenderError> {
        match std::env::var(BACKEND_ENV_VAR) {
            Ok(v) if v.trim().is_empty() => Ok(BackendPreference::default()),
            Ok(v) => BackendPreference::parse(&v).ok_or(RenderError::InvalidBackendEnv(v)),
            Err(_) => Ok(BackendPreference::default()),
        }
    }

    /// The backends a window must restrict itself to under this
    /// preference: `None` for [`BackendPreference::PreferVulkan`], which
    /// leaves the windowing toolkit's own choice (Vulkan where present) in
    /// place.
    pub fn forced_backends(self) -> Option<wgpu::Backends> {
        match self {
            BackendPreference::PreferVulkan => None,
            BackendPreference::Vulkan => Some(wgpu::Backends::VULKAN),
            BackendPreference::Gl => Some(wgpu::Backends::GL),
            BackendPreference::Any => Some(wgpu::Backends::all()),
        }
    }

    /// Backend sets to try, in order.
    fn attempts(self) -> &'static [wgpu::Backends] {
        const VULKAN: wgpu::Backends = wgpu::Backends::VULKAN;
        const GL: wgpu::Backends = wgpu::Backends::GL;
        const ALL: wgpu::Backends = wgpu::Backends::all();
        const NOT_VULKAN: wgpu::Backends = ALL.difference(VULKAN);
        match self {
            BackendPreference::PreferVulkan => &[VULKAN, NOT_VULKAN],
            BackendPreference::Vulkan => &[VULKAN],
            BackendPreference::Gl => &[GL],
            BackendPreference::Any => &[ALL],
        }
    }
}

/// An 8-bit RGBA image, rows top to bottom, tightly packed (no padding),
/// sRGB-encoded colour.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RgbaImage {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// `width * height * 4` bytes, row-major from the top-left corner.
    pub pixels: Vec<u8>,
}

impl RgbaImage {
    /// The pixel at `(x, y)` (origin top-left). Panics if out of bounds.
    pub fn pixel(&self, x: u32, y: u32) -> [u8; 4] {
        assert!(
            x < self.width && y < self.height,
            "pixel ({x}, {y}) out of bounds"
        );
        let i = (y as usize * self.width as usize + x as usize) * 4;
        [
            self.pixels[i],
            self.pixels[i + 1],
            self.pixels[i + 2],
            self.pixels[i + 3],
        ]
    }
}

/// Writes `image` as an 8-bit RGBA PNG file.
pub fn write_png(path: impl AsRef<Path>, image: &RgbaImage) -> Result<(), RenderError> {
    let expected = image.width as usize * image.height as usize * 4;
    if image.pixels.len() != expected || image.width == 0 || image.height == 0 {
        return Err(RenderError::BadImage {
            width: image.width,
            height: image.height,
            len: image.pixels.len(),
        });
    }
    let file = std::fs::File::create(path)?;
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), image.width, image.height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_source_srgb(png::SrgbRenderingIntent::Perceptual);
    let mut writer = encoder.write_header()?;
    writer.write_image_data(&image.pixels)?;
    writer.finish()?;
    Ok(())
}

struct Target {
    size: (u32, u32),
    texture: wgpu::Texture,
    view: wgpu::TextureView,
    readback: wgpu::Buffer,
    padded_bytes_per_row: u32,
}

/// A headless wgpu instance, adapter, device and queue, chosen by a
/// [`BackendPreference`]: the device behind [`OffscreenRenderer`], also
/// usable on its own (for example to render a UI offscreen on the same
/// backend selection).
pub struct HeadlessGpu {
    /// The instance the adapter came from; keep it alive with the device.
    pub instance: wgpu::Instance,
    /// The adapter in use.
    pub adapter: wgpu::Adapter,
    /// The device, created with the adapter's limits and no extra features.
    pub device: wgpu::Device,
    /// The device's queue.
    pub queue: wgpu::Queue,
}

impl HeadlessGpu {
    /// Creates a headless device using the backend selected by
    /// [`BACKEND_ENV_VAR`] (Vulkan preferred by default); `label` names the
    /// device in debug tools.
    ///
    /// Fails if no adapter is available; it never silently falls back to a
    /// no-op device.
    pub fn new(label: &str) -> Result<Self, RenderError> {
        Self::with_backend(BackendPreference::from_env()?, label)
    }

    /// Creates a headless device restricted to `preference`.
    ///
    /// The backend sets of the preference are tried in order and the first
    /// adapter found wins, so with [`BackendPreference::PreferVulkan`] any
    /// Vulkan adapter wins, including Mesa's software lavapipe; set
    /// [`BACKEND_ENV_VAR`]`=gl` on machines whose only hardware driver is
    /// OpenGL.
    pub fn with_backend(preference: BackendPreference, label: &str) -> Result<Self, RenderError> {
        let mut last_error = None;
        for &backends in preference.attempts() {
            let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
                backends,
                flags: wgpu::InstanceFlags::from_build_config().with_env(),
                backend_options: wgpu::BackendOptions::from_env_or_default(),
                ..wgpu::InstanceDescriptor::new_without_display_handle()
            });
            let options = wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::from_env()
                    .unwrap_or(wgpu::PowerPreference::HighPerformance),
                force_fallback_adapter: false,
                compatible_surface: None,
                apply_limit_buckets: false,
            };
            match pollster::block_on(instance.request_adapter(&options)) {
                Ok(adapter) => return Self::from_adapter(instance, adapter, label),
                Err(e) => {
                    log::info!("openshape-render: no adapter for {backends:?}: {e}");
                    last_error = Some((backends, e));
                }
            }
        }
        let (backends, source) = last_error.expect("at least one backend attempt");
        Err(RenderError::NoAdapter {
            preference,
            backends,
            source,
        })
    }

    fn from_adapter(
        instance: wgpu::Instance,
        adapter: wgpu::Adapter,
        label: &str,
    ) -> Result<Self, RenderError> {
        let info = adapter.get_info();
        log::info!(
            "openshape-render: using {} ({:?}, {:?}, driver {} {})",
            info.name,
            info.backend,
            info.device_type,
            info.driver,
            info.driver_info
        );
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
                label: Some(label),
                required_features: wgpu::Features::empty(),
                required_limits: adapter.limits(),
                ..Default::default()
            }))?;
        Ok(HeadlessGpu {
            instance,
            adapter,
            device,
            queue,
        })
    }
}

/// Headless renderer with its own wgpu device, for tests, thumbnails and
/// command-line export.
///
/// Renders into an `Rgba8UnormSrgb` texture through a [`SceneRenderer`]
/// (4x MSAA) and copies the result back to the CPU.
pub struct OffscreenRenderer {
    // Kept alive for the lifetime of the device.
    _instance: wgpu::Instance,
    adapter: wgpu::Adapter,
    device: wgpu::Device,
    queue: wgpu::Queue,
    scene: SceneRenderer,
    target: Option<Target>,
}

impl OffscreenRenderer {
    /// Colour format of the offscreen target (and of [`RgbaImage`] data).
    pub const COLOR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;

    /// Creates a headless renderer using the backend selected by
    /// [`BACKEND_ENV_VAR`] (Vulkan preferred by default).
    ///
    /// Fails if no adapter is available; it never silently falls back to a
    /// no-op device.
    pub fn new() -> Result<Self, RenderError> {
        Self::with_backend(BackendPreference::from_env()?)
    }

    /// Creates a headless renderer restricted to `preference` (see
    /// [`HeadlessGpu::with_backend`]).
    ///
    /// With [`BackendPreference::PreferVulkan`] any Vulkan adapter wins,
    /// including Mesa's software lavapipe; set [`BACKEND_ENV_VAR`]`=gl` on
    /// machines whose only hardware driver is OpenGL.
    pub fn with_backend(preference: BackendPreference) -> Result<Self, RenderError> {
        let HeadlessGpu {
            instance,
            adapter,
            device,
            queue,
        } = HeadlessGpu::with_backend(preference, "openshape-render offscreen")?;
        let scene = SceneRenderer::new(&device, Self::COLOR_FORMAT);
        Ok(OffscreenRenderer {
            _instance: instance,
            adapter,
            device,
            queue,
            scene,
            target: None,
        })
    }

    /// Information about the adapter in use (name, backend, driver, type).
    pub fn adapter_info(&self) -> wgpu::AdapterInfo {
        self.adapter.get_info()
    }

    /// The device this renderer created.
    pub fn device(&self) -> &wgpu::Device {
        &self.device
    }

    /// The queue this renderer created.
    pub fn queue(&self) -> &wgpu::Queue {
        &self.queue
    }

    /// The scene renderer (display mode, colours, bounds, stats).
    pub fn scene(&self) -> &SceneRenderer {
        &self.scene
    }

    /// Mutable access to the scene renderer.
    pub fn scene_mut(&mut self) -> &mut SceneRenderer {
        &mut self.scene
    }

    /// Uploads the meshes to draw; see [`SceneRenderer::set_meshes`].
    pub fn set_meshes(&mut self, meshes: &[MeshData]) -> Result<(), RenderError> {
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        self.scene.set_meshes(&self.device, meshes);
        check_scope(scope)
    }

    /// Uploads `meshes` and renders them; see [`OffscreenRenderer::render_rgba`].
    pub fn render_meshes_rgba(
        &mut self,
        meshes: &[MeshData],
        camera: &Camera,
        width: u32,
        height: u32,
    ) -> Result<RgbaImage, RenderError> {
        self.set_meshes(meshes)?;
        self.render_rgba(camera, width, height)
    }

    /// Renders the meshes from the last [`OffscreenRenderer::set_meshes`]
    /// at `width` x `height` and reads the image back (blocking).
    pub fn render_rgba(
        &mut self,
        camera: &Camera,
        width: u32,
        height: u32,
    ) -> Result<RgbaImage, RenderError> {
        let max = self.device.limits().max_texture_dimension_2d;
        if width == 0 || height == 0 || width > max || height > max {
            return Err(RenderError::InvalidSize { width, height, max });
        }
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        self.ensure_target(width, height);
        let target = self.target.as_ref().expect("target created above");

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("openshape-render offscreen frame"),
            });
        self.scene.render_to_encoder(
            &self.device,
            &self.queue,
            &mut encoder,
            camera,
            &target.view,
            (width, height),
        );
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &target.texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &target.readback,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(target.padded_bytes_per_row),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );
        self.queue.submit(Some(encoder.finish()));
        check_scope(scope)?;

        let pixels = map_and_unpad(
            &self.device,
            &target.readback,
            target.padded_bytes_per_row,
            width,
            height,
        )?;
        Ok(RgbaImage {
            width,
            height,
            pixels,
        })
    }

    /// Renders the meshes from the last [`OffscreenRenderer::set_meshes`]
    /// at `width` x `height` without reading the image back, and blocks
    /// until the GPU has finished. Used to measure frame times.
    pub fn render_frame(
        &mut self,
        camera: &Camera,
        width: u32,
        height: u32,
    ) -> Result<(), RenderError> {
        let max = self.device.limits().max_texture_dimension_2d;
        if width == 0 || height == 0 || width > max || height > max {
            return Err(RenderError::InvalidSize { width, height, max });
        }
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        self.ensure_target(width, height);
        let target = self.target.as_ref().expect("target created above");
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("openshape-render timed frame"),
            });
        self.scene.render_to_encoder(
            &self.device,
            &self.queue,
            &mut encoder,
            camera,
            &target.view,
            (width, height),
        );
        self.queue.submit(Some(encoder.finish()));
        check_scope(scope)?;
        self.device.poll(wgpu::PollType::wait_indefinitely())?;
        Ok(())
    }

    fn ensure_target(&mut self, width: u32, height: u32) {
        if self
            .target
            .as_ref()
            .is_some_and(|t| t.size == (width, height))
        {
            return;
        }
        let texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("openshape-render offscreen target"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: Self::COLOR_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
        let padded_bytes_per_row = padded_bytes_per_row(width);
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("openshape-render readback"),
            size: padded_bytes_per_row as u64 * height as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        self.target = Some(Target {
            size: (width, height),
            texture,
            view,
            readback,
            padded_bytes_per_row,
        });
    }
}

/// Bytes per row of a `width`-pixel RGBA8 row rounded up to wgpu's
/// 256-byte copy alignment.
pub fn padded_bytes_per_row(width: u32) -> u32 {
    let unpadded = width * 4;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    unpadded.div_ceil(align) * align
}

/// Reads back a single-sampled 2D texture of format `Rgba8Unorm[Srgb]` or
/// `Bgra8Unorm[Srgb]` with `COPY_SRC` usage (e.g. a viewport target rendered
/// by a [`SceneRenderer`] on a shared device) as an [`RgbaImage`], blocking
/// until the GPU is done. BGRA data is swizzled to RGBA; bytes are returned
/// as stored (sRGB-encoded for the formats [`SceneRenderer`] produces).
pub fn read_texture_rgba8(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    texture: &wgpu::Texture,
) -> Result<RgbaImage, RenderError> {
    use wgpu::TextureFormat as F;
    let format = texture.format();
    let bgra = match format {
        F::Rgba8Unorm | F::Rgba8UnormSrgb => false,
        F::Bgra8Unorm | F::Bgra8UnormSrgb => true,
        other => {
            return Err(RenderError::UnsupportedTexture(format!(
                "format {other:?} (expected RGBA8 or BGRA8)"
            )))
        }
    };
    if !texture.usage().contains(wgpu::TextureUsages::COPY_SRC)
        || texture.sample_count() != 1
        || texture.dimension() != wgpu::TextureDimension::D2
    {
        return Err(RenderError::UnsupportedTexture(
            "texture must be 2D, single-sampled and have COPY_SRC usage".into(),
        ));
    }
    let (width, height) = (texture.width(), texture.height());
    let padded = padded_bytes_per_row(width);
    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("openshape-render texture readback"),
        size: padded as u64 * height as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });
    let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("openshape-render texture readback"),
    });
    encoder.copy_texture_to_buffer(
        texture.as_image_copy(),
        wgpu::TexelCopyBufferInfo {
            buffer: &buffer,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(padded),
                rows_per_image: Some(height),
            },
        },
        texture.size(),
    );
    queue.submit(Some(encoder.finish()));
    check_scope(scope)?;
    let mut pixels = map_and_unpad(device, &buffer, padded, width, height)?;
    if bgra {
        for px in pixels.chunks_exact_mut(4) {
            px.swap(0, 2);
        }
    }
    Ok(RgbaImage {
        width,
        height,
        pixels,
    })
}

/// Maps `buffer` (filled by a texture copy with `padded` bytes per row),
/// waits for the GPU and returns the tightly packed rows.
fn map_and_unpad(
    device: &wgpu::Device,
    buffer: &wgpu::Buffer,
    padded: u32,
    width: u32,
    height: u32,
) -> Result<Vec<u8>, RenderError> {
    let (tx, rx) = std::sync::mpsc::channel();
    buffer
        .slice(..)
        .map_async(wgpu::MapMode::Read, move |result| {
            let _ = tx.send(result);
        });
    device.poll(wgpu::PollType::wait_indefinitely())?;
    rx.recv()
        .map_err(|_| RenderError::Gpu("buffer map callback was dropped".into()))??;

    let row_bytes = width as usize * 4;
    let mut pixels = Vec::with_capacity(row_bytes * height as usize);
    {
        let mapped = buffer.get_mapped_range(..)?;
        for row in mapped.chunks(padded as usize).take(height as usize) {
            pixels.extend_from_slice(&row[..row_bytes]);
        }
    }
    buffer.unmap();
    Ok(pixels)
}

fn check_scope(scope: wgpu::ErrorScopeGuard) -> Result<(), RenderError> {
    match pollster::block_on(scope.pop()) {
        Some(e) => Err(RenderError::Gpu(e.to_string())),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn padding_is_256_aligned_and_minimal() {
        assert_eq!(padded_bytes_per_row(64), 256);
        assert_eq!(padded_bytes_per_row(65), 512);
        assert_eq!(padded_bytes_per_row(333), 1536);
        assert_eq!(padded_bytes_per_row(1), 256);
    }

    #[test]
    fn backend_parsing() {
        assert_eq!(
            BackendPreference::parse("Vulkan"),
            Some(BackendPreference::Vulkan)
        );
        assert_eq!(
            BackendPreference::parse(" gl "),
            Some(BackendPreference::Gl)
        );
        assert_eq!(
            BackendPreference::parse("ANY"),
            Some(BackendPreference::Any)
        );
        assert_eq!(BackendPreference::parse("metal"), None);
    }

    #[test]
    fn a_window_is_restricted_only_when_a_backend_is_named() {
        assert_eq!(BackendPreference::PreferVulkan.forced_backends(), None);
        assert_eq!(
            BackendPreference::Vulkan.forced_backends(),
            Some(wgpu::Backends::VULKAN)
        );
        assert_eq!(
            BackendPreference::Gl.forced_backends(),
            Some(wgpu::Backends::GL)
        );
        assert_eq!(
            BackendPreference::Any.forced_backends(),
            Some(wgpu::Backends::all())
        );
    }
}
