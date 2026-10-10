//! [`SceneRenderer`]: draws [`MeshData`] with a [`Camera`] into a caller-provided
//! colour target, using a caller-provided `wgpu::Device` and `wgpu::Queue`.

use crate::camera::{Camera, Projection};
use crate::mesh::MeshData;
use bytemuck::{Pod, Zeroable};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use wgpu::util::DeviceExt;

/// Depth buffer format (reverse-Z, cleared to 0).
pub const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

/// MSAA sample count used by [`SceneRenderer::new`].
pub const DEFAULT_SAMPLE_COUNT: u32 = 4;

/// Edge lines are pulled towards the eye by this many pixels' worth of
/// depth so they win the depth test against the faces they bound, including
/// moderately oblique concave corners, without showing through solid walls.
pub const EDGE_DEPTH_BIAS_PIXELS: f64 = 2.0;

/// Slope-scaled depth bias applied to surfaces (in units of the polygon's
/// depth slope per pixel; negative pushes away under reverse-Z).
const FACE_SLOPE_DEPTH_BIAS: f32 = -1.0;

/// Default background gradient, top colour (sRGB).
pub const DEFAULT_BACKGROUND_TOP: [f32; 3] = [0.38, 0.44, 0.54];
/// Default background gradient, bottom colour (sRGB).
pub const DEFAULT_BACKGROUND_BOTTOM: [f32; 3] = [0.90, 0.91, 0.93];
/// Default edge line colour (sRGB + alpha).
pub const DEFAULT_EDGE_COLOR: [f32; 4] = [0.06, 0.07, 0.09, 1.0];

/// Default edge line width in physical pixels.
///
/// Edges are drawn as screen-space quads rather than hardware lines, so they
/// have the same width and weight on every driver (hardware 1-pixel lines
/// under MSAA come out at different strengths on different GPUs).
pub const DEFAULT_EDGE_WIDTH_PX: f32 = 1.5;

/// Upper bound on the size of one GPU vertex or index buffer created by the
/// renderer (further limited by the device's `max_buffer_size`).
const MAX_BATCH_BYTES: u64 = 128 << 20;

/// How surfaces and edges are drawn.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DisplayMode {
    /// Lit surfaces only.
    Shaded,
    /// Lit surfaces with dark B-rep edges on top (the default).
    #[default]
    ShadedWithEdges,
    /// B-rep edges only (no hidden-line removal).
    Wireframe,
}

/// Counts describing what [`SceneRenderer::set_meshes`] uploaded.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SceneStats {
    /// Meshes that contributed triangles or edges.
    pub meshes: usize,
    /// Triangles uploaded (opaque and translucent).
    pub triangles: usize,
    /// Edge line segments uploaded.
    pub edge_segments: usize,
    /// GPU batches (each is one draw call per frame).
    pub batches: usize,
}

/// GPU vertex: position, normal, sRGB colour.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct GpuVertex {
    position: [f32; 3],
    normal: [f32; 3],
    color: [u8; 4],
}

const MESH_ATTRIBUTES: [wgpu::VertexAttribute; 3] =
    wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3, 2 => Unorm8x4];
/// Per-instance attributes of an edge segment: its two end points, read
/// from the packed line-list buffer (two vertices per segment).
const EDGE_ATTRIBUTES: [wgpu::VertexAttribute; 2] =
    wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3];

/// Uniform block; layout must match `Globals` in `shader.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct Globals {
    view: [[f32; 4]; 4],
    proj: [[f32; 4]; 4],
    key_light: [f32; 4],
    sky: [f32; 4],
    ground: [f32; 4],
    edge_color: [f32; 4],
    bg_top: [f32; 4],
    bg_bottom: [f32; 4],
    params: [f32; 4],
    material: [f32; 4],
    viewport: [f32; 4],
}

/// Key light direction in view space (towards the light): from slightly
/// left of and above the eye, so the three faces of an isometric cube get
/// three clearly different shades.
const KEY_LIGHT_DIR: [f32; 3] = [-0.2, 0.28, 0.94];
const KEY_LIGHT_INTENSITY: f32 = 0.85;
/// Hemisphere ambient for normals towards world +Z (sky) and -Z (ground), linear.
const AMBIENT_SKY: [f32; 3] = [0.40, 0.41, 0.44];
const AMBIENT_GROUND: [f32; 3] = [0.24, 0.23, 0.22];
const SPECULAR_STRENGTH: f32 = 0.12;
const SHININESS: f32 = 48.0;

struct MeshBatch {
    vertices: wgpu::Buffer,
    indices: wgpu::Buffer,
    index_count: u32,
}

struct EdgeBatch {
    vertices: wgpu::Buffer,
    vertex_count: u32,
}

struct Attachments {
    size: (u32, u32),
    /// Multisampled colour target, absent when `sample_count == 1`.
    msaa_color: Option<wgpu::TextureView>,
    depth: wgpu::TextureView,
}

/// Renders meshes into a caller-provided colour target.
///
/// The renderer owns its pipelines, mesh buffers, a uniform buffer and its
/// depth / MSAA attachments, but borrows the device and queue on each call,
/// so it can share the device of an egui-wgpu application.
///
/// Typical use: [`SceneRenderer::set_meshes`] when geometry changes, then
/// [`SceneRenderer::render`] every frame. Rendering is fully headless; the
/// target is any `TextureView` of the format given to [`SceneRenderer::new`]
/// with `RENDER_ATTACHMENT` usage and sample count 1.
///
/// Translucent meshes (`color[3] < 1`) are blended after opaque ones without
/// sorting and without depth writes.
pub struct SceneRenderer {
    color_format: wgpu::TextureFormat,
    sample_count: u32,
    globals: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    background_pipeline: wgpu::RenderPipeline,
    opaque_pipeline: wgpu::RenderPipeline,
    transparent_pipeline: wgpu::RenderPipeline,
    edge_pipeline: wgpu::RenderPipeline,
    opaque: Vec<MeshBatch>,
    transparent: Vec<MeshBatch>,
    edges: Vec<EdgeBatch>,
    bounds: Option<([f32; 3], [f32; 3])>,
    stats: SceneStats,
    attachments: Option<Attachments>,
    display_mode: DisplayMode,
    background_top: [f32; 3],
    background_bottom: [f32; 3],
    edge_color: [f32; 4],
    edge_width: f32,
}

impl SceneRenderer {
    /// Creates a renderer drawing into targets of `color_format` with 4x MSAA.
    ///
    /// If `color_format` is not an sRGB format the shader encodes its output
    /// to sRGB itself, so colours look the same either way. Prefer an sRGB
    /// format when possible: with a non-sRGB target the MSAA resolve averages
    /// encoded values, which makes antialiased edge lines look heavier.
    pub fn new(device: &wgpu::Device, color_format: wgpu::TextureFormat) -> Self {
        Self::with_sample_count(device, color_format, DEFAULT_SAMPLE_COUNT)
    }

    /// Like [`SceneRenderer::new`] with an explicit MSAA sample count
    /// (1 disables multisampling; 4 is supported everywhere wgpu runs).
    pub fn with_sample_count(
        device: &wgpu::Device,
        color_format: wgpu::TextureFormat,
        sample_count: u32,
    ) -> Self {
        let sample_count = sample_count.max(1);
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("openshape-render shader"),
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(include_str!("shader.wgsl"))),
        });
        let globals = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("openshape-render globals"),
            size: std::mem::size_of::<Globals>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("openshape-render globals layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: wgpu::BufferSize::new(std::mem::size_of::<Globals>() as u64),
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("openshape-render globals"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: globals.as_entire_binding(),
            }],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("openshape-render pipeline layout"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });

        let builder = PipelineBuilder {
            device,
            layout: &layout,
            module: &module,
            color_format,
            sample_count,
        };
        let background_pipeline = builder.build(PipelineKind::Background);
        let opaque_pipeline = builder.build(PipelineKind::Opaque);
        let transparent_pipeline = builder.build(PipelineKind::Transparent);
        let edge_pipeline = builder.build(PipelineKind::Edges);

        SceneRenderer {
            color_format,
            sample_count,
            globals,
            bind_group,
            background_pipeline,
            opaque_pipeline,
            transparent_pipeline,
            edge_pipeline,
            opaque: Vec::new(),
            transparent: Vec::new(),
            edges: Vec::new(),
            bounds: None,
            stats: SceneStats::default(),
            attachments: None,
            display_mode: DisplayMode::default(),
            background_top: DEFAULT_BACKGROUND_TOP,
            background_bottom: DEFAULT_BACKGROUND_BOTTOM,
            edge_color: DEFAULT_EDGE_COLOR,
            edge_width: DEFAULT_EDGE_WIDTH_PX,
        }
    }

    /// Colour format of the targets this renderer draws into.
    pub fn color_format(&self) -> wgpu::TextureFormat {
        self.color_format
    }

    /// MSAA sample count of the internal attachments.
    pub fn sample_count(&self) -> u32 {
        self.sample_count
    }

    /// Current display mode.
    pub fn display_mode(&self) -> DisplayMode {
        self.display_mode
    }

    /// Sets how surfaces and edges are drawn.
    pub fn set_display_mode(&mut self, mode: DisplayMode) {
        self.display_mode = mode;
    }

    /// Sets the vertical background gradient (sRGB colours, top and bottom).
    pub fn set_background(&mut self, top: [f32; 3], bottom: [f32; 3]) {
        self.background_top = top;
        self.background_bottom = bottom;
    }

    /// Sets the edge line colour (sRGB + alpha).
    pub fn set_edge_color(&mut self, color: [f32; 4]) {
        self.edge_color = color;
    }

    /// Edge line width in physical pixels.
    pub fn edge_width(&self) -> f32 {
        self.edge_width
    }

    /// Sets the edge line width in physical pixels (clamped to 0.5..=16;
    /// non-finite values are ignored). Scale by the display's pixels per
    /// point for a constant apparent width on high-DPI screens.
    pub fn set_edge_width(&mut self, px: f32) {
        if px.is_finite() {
            self.edge_width = px.clamp(0.5, 16.0);
        }
    }

    /// Bounds `(min, max)` of everything uploaded by the last
    /// [`SceneRenderer::set_meshes`], or `None` if the scene is empty.
    /// Use it to [`Camera::fit`] the camera to the scene.
    pub fn bounds(&self) -> Option<([f32; 3], [f32; 3])> {
        self.bounds
    }

    /// What the last [`SceneRenderer::set_meshes`] uploaded.
    pub fn stats(&self) -> SceneStats {
        self.stats
    }

    /// Replaces the displayed geometry. Call only when geometry or colours
    /// change; this packs all meshes into a few large GPU buffers so a frame
    /// costs a handful of draw calls regardless of the mesh count.
    ///
    /// Invalid input is tolerated: missing normals are computed, triangles
    /// with out-of-range indices and incomplete trailing triangles or edge
    /// segments are dropped (with a log warning), fully transparent meshes
    /// only contribute their edges.
    pub fn set_meshes(&mut self, device: &wgpu::Device, meshes: &[MeshData]) {
        let max_bytes = device.limits().max_buffer_size.min(MAX_BATCH_BYTES);
        let mut opaque = MeshBatcher::new(max_bytes);
        let mut transparent = MeshBatcher::new(max_bytes);
        let mut edges = EdgeBatcher::new(max_bytes);
        let mut stats = SceneStats::default();
        let mut bounds: Option<([f32; 3], [f32; 3])> = None;

        for (mesh_index, mesh) in meshes.iter().enumerate() {
            let mesh = sanitize(mesh_index, mesh);
            if mesh.is_empty() {
                continue;
            }
            stats.meshes += 1;
            if let Some((lo, hi)) = mesh.bounds() {
                bounds = Some(match bounds {
                    None => (lo, hi),
                    Some((a, b)) => (
                        std::array::from_fn(|i| a[i].min(lo[i])),
                        std::array::from_fn(|i| b[i].max(hi[i])),
                    ),
                });
            }
            let alpha = mesh.color[3].clamp(0.0, 1.0);
            if alpha > 0.0 && mesh.triangle_count() > 0 {
                stats.triangles += mesh.triangle_count();
                let batcher = if alpha >= 0.999 {
                    &mut opaque
                } else {
                    &mut transparent
                };
                batcher.push(device, &mesh);
            }
            stats.edge_segments += mesh.edge_segment_count();
            edges.push(device, &mesh.edge_lines);
        }

        self.opaque = opaque.finish(device);
        self.transparent = transparent.finish(device);
        self.edges = edges.finish(device);
        stats.batches = self.opaque.len() + self.transparent.len() + self.edges.len();
        self.stats = stats;
        self.bounds = bounds;
        log::debug!("openshape-render: uploaded {stats:?}");
    }

    /// Renders one frame into `target` (a view of a texture of
    /// [`SceneRenderer::color_format`], sample count 1, `size` = its
    /// dimensions in pixels), submitting the work to `queue`.
    ///
    /// The depth buffer and MSAA attachments are (re)created when `size`
    /// changes. A zero-sized `size` renders nothing.
    pub fn render(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        camera: &Camera,
        target: &wgpu::TextureView,
        size: (u32, u32),
    ) {
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("openshape-render frame"),
        });
        self.render_to_encoder(device, queue, &mut encoder, camera, target, size);
        queue.submit(Some(encoder.finish()));
    }

    /// Like [`SceneRenderer::render`] but records into a caller-provided
    /// encoder (e.g. egui-wgpu's `prepare` callback encoder) instead of
    /// submitting. The uniform update is written through `queue` and takes
    /// effect at the next submission, so record at most one frame per
    /// renderer per submission.
    pub fn render_to_encoder(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        encoder: &mut wgpu::CommandEncoder,
        camera: &Camera,
        target: &wgpu::TextureView,
        size: (u32, u32),
    ) {
        if size.0 == 0 || size.1 == 0 {
            return;
        }
        self.ensure_attachments(device, size);
        let globals = self.globals_for(camera, size);
        queue.write_buffer(&self.globals, 0, bytemuck::bytes_of(&globals));

        let Some(attachments) = self.attachments.as_ref() else {
            return;
        };
        let color_attachment = match &attachments.msaa_color {
            Some(msaa) => wgpu::RenderPassColorAttachment {
                view: msaa,
                depth_slice: None,
                resolve_target: Some(target),
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                    store: wgpu::StoreOp::Discard,
                },
            },
            None => wgpu::RenderPassColorAttachment {
                view: target,
                depth_slice: None,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                    store: wgpu::StoreOp::Store,
                },
            },
        };
        let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("openshape-render scene"),
            color_attachments: &[Some(color_attachment)],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: &attachments.depth,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(0.0),
                    store: wgpu::StoreOp::Discard,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
            multiview_mask: None,
        });
        pass.set_bind_group(0, &self.bind_group, &[]);

        pass.set_pipeline(&self.background_pipeline);
        pass.draw(0..3, 0..1);

        let draw_surfaces = self.display_mode != DisplayMode::Wireframe;
        let draw_edges = self.display_mode != DisplayMode::Shaded;
        if draw_surfaces && !self.opaque.is_empty() {
            pass.set_pipeline(&self.opaque_pipeline);
            draw_mesh_batches(&mut pass, &self.opaque);
        }
        if draw_edges && !self.edges.is_empty() {
            pass.set_pipeline(&self.edge_pipeline);
            for batch in &self.edges {
                // One instance per segment, six vertices (two triangles) each.
                pass.set_vertex_buffer(0, batch.vertices.slice(..));
                pass.draw(0..6, 0..batch.vertex_count / 2);
            }
        }
        if draw_surfaces && !self.transparent.is_empty() {
            pass.set_pipeline(&self.transparent_pipeline);
            draw_mesh_batches(&mut pass, &self.transparent);
        }
    }

    fn ensure_attachments(&mut self, device: &wgpu::Device, size: (u32, u32)) {
        if self.attachments.as_ref().is_some_and(|a| a.size == size) {
            return;
        }
        let extent = wgpu::Extent3d {
            width: size.0,
            height: size.1,
            depth_or_array_layers: 1,
        };
        let make = |label: &str, format: wgpu::TextureFormat| {
            device
                .create_texture(&wgpu::TextureDescriptor {
                    label: Some(label),
                    size: extent,
                    mip_level_count: 1,
                    sample_count: self.sample_count,
                    dimension: wgpu::TextureDimension::D2,
                    format,
                    usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
                    view_formats: &[],
                })
                .create_view(&wgpu::TextureViewDescriptor::default())
        };
        let msaa_color =
            (self.sample_count > 1).then(|| make("openshape-render msaa color", self.color_format));
        let depth = make("openshape-render depth", DEPTH_FORMAT);
        self.attachments = Some(Attachments {
            size,
            msaa_color,
            depth,
        });
    }

    fn globals_for(&self, camera: &Camera, size: (u32, u32)) -> Globals {
        let aspect = size.0 as f32 / size.1 as f32;
        let view = camera.view_matrix();
        let proj = match self.bounds {
            Some((lo, hi)) => camera.proj_matrix_with_bounds(aspect, lo, hi),
            None => camera.proj_matrix(aspect),
        };
        let (perspective, edge_scale, edge_offset) = match camera.projection {
            Projection::Perspective { .. } => {
                // Relative pull towards the eye: one pixel spans
                // 2 tan(fov/2) / height world units per unit of depth.
                let pixel = 2.0 * (camera.fov_y() * 0.5).tan() / size.1 as f64;
                let rel = EDGE_DEPTH_BIAS_PIXELS * pixel;
                (1.0, (1.0 - rel).clamp(0.9, 1.0) as f32, 0.0)
            }
            Projection::Orthographic => {
                let bias = EDGE_DEPTH_BIAS_PIXELS * camera.world_units_per_pixel(size.1 as f32);
                (0.0, 1.0, bias as f32)
            }
        };
        let [lx, ly, lz] = KEY_LIGHT_DIR;
        let len = (lx * lx + ly * ly + lz * lz).sqrt();
        let srgb_encode = if self.color_format.is_srgb() {
            0.0
        } else {
            1.0
        };
        let edge = srgb_to_linear3([self.edge_color[0], self.edge_color[1], self.edge_color[2]]);
        Globals {
            view: view.to_cols_array_2d(),
            proj: proj.to_cols_array_2d(),
            key_light: [lx / len, ly / len, lz / len, KEY_LIGHT_INTENSITY],
            sky: [AMBIENT_SKY[0], AMBIENT_SKY[1], AMBIENT_SKY[2], 0.0],
            ground: [AMBIENT_GROUND[0], AMBIENT_GROUND[1], AMBIENT_GROUND[2], 0.0],
            edge_color: [edge[0], edge[1], edge[2], self.edge_color[3]],
            bg_top: with_w(srgb_to_linear3(self.background_top), 1.0),
            bg_bottom: with_w(srgb_to_linear3(self.background_bottom), 1.0),
            params: [edge_scale, edge_offset, srgb_encode, perspective],
            material: [SPECULAR_STRENGTH, SHININESS, 0.0, 0.0],
            viewport: [size.0 as f32, size.1 as f32, 0.5 * self.edge_width, 0.0],
        }
    }
}

fn draw_mesh_batches(pass: &mut wgpu::RenderPass<'_>, batches: &[MeshBatch]) {
    for batch in batches {
        pass.set_vertex_buffer(0, batch.vertices.slice(..));
        pass.set_index_buffer(batch.indices.slice(..), wgpu::IndexFormat::Uint32);
        pass.draw_indexed(0..batch.index_count, 0, 0..1);
    }
}

fn with_w(v: [f32; 3], w: f32) -> [f32; 4] {
    [v[0], v[1], v[2], w]
}

/// sRGB transfer function decode, per channel.
pub(crate) fn srgb_to_linear(c: f32) -> f32 {
    let c = c.clamp(0.0, 1.0);
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

fn srgb_to_linear3(c: [f32; 3]) -> [f32; 3] {
    c.map(srgb_to_linear)
}

/// Returns the mesh with normals present, indices in range and complete
/// primitives only. Borrows when nothing needs fixing.
fn sanitize(mesh_index: usize, mesh: &MeshData) -> Cow<'_, MeshData> {
    let n = mesh.positions.len();
    let bad_indices = mesh.indices.len() % 3 != 0 || mesh.indices.iter().any(|&i| i as usize >= n);
    let odd_edges = mesh.edge_lines.len() % 2 != 0;
    if mesh.has_normals() && !bad_indices && !odd_edges {
        return Cow::Borrowed(mesh);
    }
    let mut fixed = mesh.clone();
    if bad_indices {
        let before = fixed.triangle_count();
        fixed.indices = mesh
            .indices
            .chunks_exact(3)
            .filter(|t| t.iter().all(|&i| (i as usize) < n))
            .flatten()
            .copied()
            .collect();
        log::warn!(
            "openshape-render: mesh {mesh_index}: dropped {} triangle(s) with invalid indices",
            before - fixed.triangle_count()
        );
    }
    if odd_edges {
        fixed.edge_lines.pop();
        log::warn!(
            "openshape-render: mesh {mesh_index}: odd edge_lines length, dropped last vertex"
        );
    }
    if !fixed.has_normals() {
        if !mesh.normals.is_empty() {
            log::warn!(
                "openshape-render: mesh {mesh_index}: {} normals for {n} positions, recomputing",
                mesh.normals.len()
            );
        }
        fixed.compute_smooth_normals();
    }
    Cow::Owned(fixed)
}

fn color_to_unorm8(c: [f32; 4]) -> [u8; 4] {
    c.map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8)
}

/// Accumulates triangle meshes into vertex/index buffers of bounded size.
struct MeshBatcher {
    max_vertices: usize,
    max_indices: usize,
    vertices: Vec<GpuVertex>,
    indices: Vec<u32>,
    done: Vec<MeshBatch>,
}

impl MeshBatcher {
    fn new(max_bytes: u64) -> Self {
        let max_vertices = (max_bytes / std::mem::size_of::<GpuVertex>() as u64) as usize;
        let max_indices = (max_bytes / 4) as usize;
        MeshBatcher {
            max_vertices,
            max_indices,
            vertices: Vec::new(),
            indices: Vec::new(),
            done: Vec::new(),
        }
    }

    fn push(&mut self, device: &wgpu::Device, mesh: &MeshData) {
        let nv = mesh.positions.len();
        let ni = mesh.indices.len();
        if nv > self.max_vertices || ni > self.max_indices {
            log::error!(
                "openshape-render: mesh with {nv} vertices / {ni} indices exceeds the GPU buffer \
                 limit and is not drawn; split it into smaller meshes"
            );
            return;
        }
        if self.vertices.len() + nv > self.max_vertices
            || self.indices.len() + ni > self.max_indices
        {
            self.flush(device);
        }
        let base = self.vertices.len() as u32;
        let color = color_to_unorm8(mesh.color);
        self.vertices.extend(
            mesh.positions
                .iter()
                .zip(&mesh.normals)
                .map(|(p, n)| GpuVertex {
                    position: *p,
                    normal: *n,
                    color,
                }),
        );
        self.indices.extend(mesh.indices.iter().map(|&i| i + base));
    }

    fn flush(&mut self, device: &wgpu::Device) {
        if self.indices.is_empty() {
            self.vertices.clear();
            return;
        }
        let vertices = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("openshape-render mesh vertices"),
            contents: bytemuck::cast_slice(&self.vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });
        let indices = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("openshape-render mesh indices"),
            contents: bytemuck::cast_slice(&self.indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        self.done.push(MeshBatch {
            vertices,
            indices,
            index_count: self.indices.len() as u32,
        });
        self.vertices.clear();
        self.indices.clear();
    }

    fn finish(mut self, device: &wgpu::Device) -> Vec<MeshBatch> {
        self.flush(device);
        self.done
    }
}

/// Accumulates edge line vertices into vertex buffers of bounded size.
struct EdgeBatcher {
    max_vertices: usize,
    vertices: Vec<[f32; 3]>,
    done: Vec<EdgeBatch>,
}

impl EdgeBatcher {
    fn new(max_bytes: u64) -> Self {
        // Keep an even count so segments never straddle two batches.
        let max_vertices = ((max_bytes / 12) as usize) & !1;
        EdgeBatcher {
            max_vertices,
            vertices: Vec::new(),
            done: Vec::new(),
        }
    }

    fn push(&mut self, device: &wgpu::Device, lines: &[[f32; 3]]) {
        for chunk in lines.chunks(self.max_vertices.max(2)) {
            if self.vertices.len() + chunk.len() > self.max_vertices {
                self.flush(device);
            }
            self.vertices.extend_from_slice(chunk);
        }
    }

    fn flush(&mut self, device: &wgpu::Device) {
        if self.vertices.is_empty() {
            return;
        }
        let vertices = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("openshape-render edge vertices"),
            contents: bytemuck::cast_slice(&self.vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });
        self.done.push(EdgeBatch {
            vertices,
            vertex_count: self.vertices.len() as u32,
        });
        self.vertices.clear();
    }

    fn finish(mut self, device: &wgpu::Device) -> Vec<EdgeBatch> {
        self.flush(device);
        self.done
    }
}

#[derive(Clone, Copy)]
enum PipelineKind {
    Background,
    Opaque,
    Transparent,
    Edges,
}

struct PipelineBuilder<'a> {
    device: &'a wgpu::Device,
    layout: &'a wgpu::PipelineLayout,
    module: &'a wgpu::ShaderModule,
    color_format: wgpu::TextureFormat,
    sample_count: u32,
}

impl PipelineBuilder<'_> {
    fn build(&self, kind: PipelineKind) -> wgpu::RenderPipeline {
        let mesh_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<GpuVertex>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &MESH_ATTRIBUTES,
        };
        let edge_layout = wgpu::VertexBufferLayout {
            array_stride: 24,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &EDGE_ATTRIBUTES,
        };
        let (label, vs, fs, buffers, topology, depth_write, depth_compare, blend) = match kind {
            PipelineKind::Background => (
                "openshape-render background",
                "vs_background",
                "fs_background",
                None,
                wgpu::PrimitiveTopology::TriangleList,
                false,
                wgpu::CompareFunction::Always,
                None,
            ),
            PipelineKind::Opaque => (
                "openshape-render opaque",
                "vs_mesh",
                "fs_mesh",
                Some(mesh_layout),
                wgpu::PrimitiveTopology::TriangleList,
                true,
                wgpu::CompareFunction::Greater,
                None,
            ),
            PipelineKind::Transparent => (
                "openshape-render transparent",
                "vs_mesh",
                "fs_mesh",
                Some(mesh_layout),
                wgpu::PrimitiveTopology::TriangleList,
                false,
                wgpu::CompareFunction::Greater,
                Some(wgpu::BlendState::ALPHA_BLENDING),
            ),
            PipelineKind::Edges => (
                "openshape-render edges",
                "vs_edge",
                "fs_edge",
                Some(edge_layout),
                wgpu::PrimitiveTopology::TriangleList,
                true,
                wgpu::CompareFunction::GreaterEqual,
                Some(wgpu::BlendState::ALPHA_BLENDING),
            ),
        };
        // Zero or one vertex buffer slot.
        let vertex_buffers: Vec<Option<wgpu::VertexBufferLayout<'_>>> =
            buffers.into_iter().map(Some).collect();
        self.device
            .create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(self.layout),
                vertex: wgpu::VertexState {
                    module: self.module,
                    entry_point: Some(vs),
                    compilation_options: Default::default(),
                    buffers: &vertex_buffers,
                },
                primitive: wgpu::PrimitiveState {
                    topology,
                    strip_index_format: None,
                    front_face: wgpu::FrontFace::Ccw,
                    // CAD shells may be open or inconsistently wound.
                    cull_mode: None,
                    unclipped_depth: false,
                    polygon_mode: wgpu::PolygonMode::Fill,
                    conservative: false,
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: DEPTH_FORMAT,
                    depth_write_enabled: Some(depth_write),
                    depth_compare: Some(depth_compare),
                    stencil: wgpu::StencilState::default(),
                    bias: match kind {
                        // Push surfaces back in proportion to their depth
                        // slope (reverse-Z: negative is farther), so edges on
                        // nearly edge-on faces still win the depth test.
                        PipelineKind::Opaque | PipelineKind::Transparent => wgpu::DepthBiasState {
                            constant: 0,
                            slope_scale: FACE_SLOPE_DEPTH_BIAS,
                            clamp: 0.0,
                        },
                        PipelineKind::Background | PipelineKind::Edges => {
                            wgpu::DepthBiasState::default()
                        }
                    },
                }),
                multisample: wgpu::MultisampleState {
                    count: self.sample_count,
                    mask: !0,
                    alpha_to_coverage_enabled: false,
                },
                fragment: Some(wgpu::FragmentState {
                    module: self.module,
                    entry_point: Some(fs),
                    compilation_options: Default::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: self.color_format,
                        blend,
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
    }
}
