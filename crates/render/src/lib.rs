//! OpenShape 3D viewport renderer (wgpu).
//!
//! This crate turns neutral mesh data into shaded, edge-outlined images. It
//! depends neither on the CAD kernel nor on the UI toolkit:
//!
//! - [`MeshData`]: triangles + normals + B-rep edge lines + colour, produced
//!   by whoever tessellates the model.
//! - [`Camera`]: an orbit (turntable) camera, Z up, reverse-Z projections
//!   with clip planes derived from the scene bounds.
//! - [`SceneRenderer`]: draws meshes with a caller-provided `wgpu::Device` and
//!   `Queue` (e.g. egui-wgpu's) into a caller-provided colour target, with its
//!   own depth buffer and 4x MSAA.
//! - [`OffscreenRenderer`]: creates its own headless device (Vulkan preferred,
//!   see [`BACKEND_ENV_VAR`]) and reads frames back as [`RgbaImage`]s;
//!   [`write_png`] saves them.
//!
//! ```no_run
//! use openshape_render::{Camera, MeshData, OffscreenRenderer, write_png};
//!
//! let cube = MeshData::cuboid([-0.5; 3], [0.5; 3], MeshData::DEFAULT_COLOR);
//! let camera = Camera::isometric().fitted([-0.5; 3], [0.5; 3]);
//! let mut renderer = OffscreenRenderer::new()?;
//! let image = renderer.render_meshes_rgba(&[cube], &camera, 640, 480)?;
//! write_png("cube.png", &image)?;
//! # Ok::<(), openshape_render::RenderError>(())
//! ```

pub mod camera;
pub mod mesh;
pub mod offscreen;
pub mod scene;

pub use camera::{Camera, Projection, ViewPreset, DEFAULT_FOV_Y_RAD, ORBIT_RADIANS_PER_PIXEL};
pub use mesh::MeshData;
pub use offscreen::{
    padded_bytes_per_row, read_texture_rgba8, write_png, BackendPreference, OffscreenRenderer,
    RgbaImage, BACKEND_ENV_VAR,
};
pub use scene::{DisplayMode, SceneRenderer, SceneStats, DEFAULT_SAMPLE_COUNT, DEPTH_FORMAT};

/// Re-export of the wgpu version this crate is built against, so callers can
/// name matching types.
pub use wgpu;

/// Errors from headless rendering and image output.
#[derive(Debug, thiserror::Error)]
pub enum RenderError {
    /// [`BACKEND_ENV_VAR`] holds an unrecognised value.
    #[error("{BACKEND_ENV_VAR}={0:?} is not one of: vulkan, gl, any")]
    InvalidBackendEnv(String),
    /// No GPU adapter (hardware or software) was found.
    #[error("no GPU adapter found ({preference:?}; last tried {backends:?}): {source}")]
    NoAdapter {
        /// The requested backend preference.
        preference: BackendPreference,
        /// The backend set of the last attempt.
        backends: wgpu::Backends,
        /// wgpu's explanation.
        #[source]
        source: wgpu::RequestAdapterError,
    },
    /// The adapter refused to create a device.
    #[error("failed to create GPU device: {0}")]
    RequestDevice(#[from] wgpu::RequestDeviceError),
    /// Requested image size is zero or exceeds the device limit.
    #[error("invalid image size {width}x{height} (each side must be 1..={max})")]
    InvalidSize {
        /// Requested width.
        width: u32,
        /// Requested height.
        height: u32,
        /// Device's maximum 2D texture dimension.
        max: u32,
    },
    /// A texture that [`read_texture_rgba8`] cannot read back.
    #[error("unsupported texture for readback: {0}")]
    UnsupportedTexture(String),
    /// A wgpu validation error or other GPU-side failure.
    #[error("GPU error: {0}")]
    Gpu(String),
    /// Mapping the readback buffer failed.
    #[error("readback buffer mapping failed: {0}")]
    BufferMap(#[from] wgpu::BufferAsyncError),
    /// Accessing the mapped readback buffer failed.
    #[error("readback buffer access failed: {0}")]
    MapRange(#[from] wgpu::MapRangeError),
    /// Waiting for the device failed.
    #[error("device poll failed: {0}")]
    Poll(#[from] wgpu::PollError),
    /// An [`RgbaImage`] whose pixel buffer does not match its dimensions.
    #[error("image {width}x{height} has {len} bytes of pixel data")]
    BadImage {
        /// Image width.
        width: u32,
        /// Image height.
        height: u32,
        /// Actual pixel buffer length.
        len: usize,
    },
    /// File I/O failed.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// PNG encoding failed.
    #[error("PNG encoding failed: {0}")]
    Png(#[from] png::EncodingError),
}
