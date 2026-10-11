//! The texture the 3D view is rendered into and shown through egui.
//!
//! egui-wgpu samples native textures as gamma-encoded (`Rgba8Unorm`) data.
//! The scene renderer should still render into an sRGB target, so that
//! blending and the multisample resolve happen in linear space. Where the
//! device can reinterpret a texture's sRGB-ness (`DownlevelFlags::
//! VIEW_FORMATS`: Vulkan, Metal, DX12), the texture is `Rgba8UnormSrgb`, the
//! scene renders through an sRGB view and egui reads the same bytes through
//! an `Rgba8Unorm` view. Otherwise (GL) the texture is `Rgba8Unorm` and the
//! scene renderer encodes sRGB in its shader.

use openshape_render::wgpu;

/// The format the scene renderer must be created with on this adapter.
pub fn scene_format(adapter: &wgpu::Adapter) -> wgpu::TextureFormat {
    if reinterprets(adapter) {
        wgpu::TextureFormat::Rgba8UnormSrgb
    } else {
        wgpu::TextureFormat::Rgba8Unorm
    }
}

fn reinterprets(adapter: &wgpu::Adapter) -> bool {
    adapter
        .get_downlevel_capabilities()
        .flags
        .contains(wgpu::DownlevelFlags::VIEW_FORMATS)
}

/// A viewport texture with its two views.
#[derive(Debug)]
pub struct ViewportTexture {
    /// The texture.
    pub texture: wgpu::Texture,
    /// The view the scene renders into (format [`scene_format`]).
    pub render_view: wgpu::TextureView,
    /// The view registered with egui (`Rgba8Unorm`).
    pub egui_view: wgpu::TextureView,
    /// Size in pixels.
    pub size: (u32, u32),
}

impl ViewportTexture {
    /// A texture of `size` pixels for a scene renderer created with
    /// `format` (from [`scene_format`]).
    pub fn new(device: &wgpu::Device, format: wgpu::TextureFormat, size: (u32, u32)) -> Self {
        let srgb = format == wgpu::TextureFormat::Rgba8UnormSrgb;
        let view_formats: &[wgpu::TextureFormat] = if srgb {
            &[wgpu::TextureFormat::Rgba8Unorm]
        } else {
            &[]
        };
        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("openshape viewport"),
            size: wgpu::Extent3d {
                width: size.0,
                height: size.1,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                | wgpu::TextureUsages::TEXTURE_BINDING
                | wgpu::TextureUsages::COPY_SRC,
            view_formats,
        });
        let render_view = texture.create_view(&wgpu::TextureViewDescriptor::default());
        let egui_view = texture.create_view(&wgpu::TextureViewDescriptor {
            format: Some(wgpu::TextureFormat::Rgba8Unorm),
            ..Default::default()
        });
        ViewportTexture {
            texture,
            render_view,
            egui_view,
            size,
        }
    }
}
