//! Shared helpers for the GPU integration tests.
#![allow(dead_code)]

use glam::Vec3;
use openshape_render::{write_png, Camera, MeshData, OffscreenRenderer, RgbaImage};
use std::f32::consts::TAU;

/// Creates a headless renderer or panics loudly: these tests must never pass
/// silently on a machine without a usable adapter.
pub fn renderer() -> OffscreenRenderer {
    match OffscreenRenderer::new() {
        Ok(r) => {
            let info = r.adapter_info();
            eprintln!(
                "adapter: {} ({:?}, {:?}, {} {})",
                info.name, info.backend, info.device_type, info.driver, info.driver_info
            );
            r
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

/// Unit cube centred at the origin.
pub fn unit_cube() -> MeshData {
    MeshData::cuboid([-0.5; 3], [0.5; 3], MeshData::DEFAULT_COLOR)
}

/// Pixel coordinates of world point `p` for `camera` at `width` x `height`.
pub fn project(camera: &Camera, p: [f32; 3], width: u32, height: u32) -> (u32, u32) {
    let vp = camera.view_proj(width as f32 / height as f32);
    let ndc = vp.project_point3(Vec3::from_array(p));
    let x = ((ndc.x * 0.5 + 0.5) * width as f32).floor();
    let y = ((0.5 - ndc.y * 0.5) * height as f32).floor();
    (
        (x.max(0.0) as u32).min(width - 1),
        (y.max(0.0) as u32).min(height - 1),
    )
}

/// Rec. 709 luma of an sRGB pixel, 0..=255.
pub fn luma(p: [u8; 4]) -> f32 {
    0.2126 * p[0] as f32 + 0.7152 * p[1] as f32 + 0.0722 * p[2] as f32
}

/// Sum of absolute RGB channel differences.
pub fn rgb_diff(a: [u8; 4], b: [u8; 4]) -> u32 {
    (0..3)
        .map(|i| (a[i] as i32 - b[i] as i32).unsigned_abs())
        .sum()
}

/// Number of pixels whose RGB differs by more than `tolerance` (channel sum)
/// inside the rectangle `[x0, x1) x [y0, y1)`.
pub fn count_different(
    a: &RgbaImage,
    b: &RgbaImage,
    (x0, y0, x1, y1): (u32, u32, u32, u32),
    tolerance: u32,
) -> usize {
    let mut n = 0;
    for y in y0..y1 {
        for x in x0..x1 {
            if rgb_diff(a.pixel(x, y), b.pixel(x, y)) > tolerance {
                n += 1;
            }
        }
    }
    n
}

/// Closed cylinder along Z: smooth side normals, flat caps, and the two cap
/// circles as edges.
pub fn cylinder(center: [f32; 3], radius: f32, height: f32, segments: u32) -> MeshData {
    let c = Vec3::from_array(center);
    let mut m = MeshData {
        color: [0.85, 0.62, 0.35, 1.0],
        ..MeshData::default()
    };
    let ring = |i: u32| {
        let a = TAU * i as f32 / segments as f32;
        Vec3::new(a.cos(), a.sin(), 0.0)
    };
    let (z0, z1) = (-height * 0.5, height * 0.5);
    // Side: two vertices per segment boundary, smooth normals.
    for i in 0..=segments {
        let d = ring(i % segments);
        for z in [z0, z1] {
            m.positions.push((c + d * radius + Vec3::Z * z).to_array());
            m.normals.push(d.to_array());
        }
    }
    for i in 0..segments {
        let b = 2 * i;
        m.indices
            .extend_from_slice(&[b, b + 2, b + 3, b, b + 3, b + 1]);
    }
    // Caps: flat normals, triangle fans.
    for (z, n) in [(z1, Vec3::Z), (z0, -Vec3::Z)] {
        let center_index = m.positions.len() as u32;
        m.positions.push((c + Vec3::Z * z).to_array());
        m.normals.push(n.to_array());
        for i in 0..segments {
            m.positions
                .push((c + ring(i) * radius + Vec3::Z * z).to_array());
            m.normals.push(n.to_array());
        }
        for i in 0..segments {
            let a = center_index + 1 + i;
            let b = center_index + 1 + (i + 1) % segments;
            if n.z > 0.0 {
                m.indices.extend_from_slice(&[center_index, a, b]);
            } else {
                m.indices.extend_from_slice(&[center_index, b, a]);
            }
        }
        for i in 0..segments {
            m.edge_lines
                .push((c + ring(i) * radius + Vec3::Z * z).to_array());
            m.edge_lines
                .push((c + ring(i + 1) * radius + Vec3::Z * z).to_array());
        }
    }
    m
}

/// UV sphere with `2 * slices * (stacks - 1)` triangles, smooth normals and
/// its equator as an edge loop.
pub fn uv_sphere(radius: f32, slices: u32, stacks: u32) -> MeshData {
    let mut m = MeshData {
        color: MeshData::DEFAULT_COLOR,
        ..MeshData::default()
    };
    for j in 0..=stacks {
        let theta = std::f32::consts::PI * j as f32 / stacks as f32;
        for i in 0..=slices {
            let phi = TAU * i as f32 / slices as f32;
            let n = Vec3::new(
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            );
            m.positions.push((n * radius).to_array());
            m.normals.push(n.to_array());
        }
    }
    let row = slices + 1;
    for j in 0..stacks {
        for i in 0..slices {
            let a = j * row + i;
            let b = a + row;
            if j != 0 {
                m.indices.extend_from_slice(&[a, b, a + 1]);
            }
            if j != stacks - 1 {
                m.indices.extend_from_slice(&[a + 1, b, b + 1]);
            }
        }
    }
    for i in 0..slices {
        for k in [i, i + 1] {
            let phi = TAU * k as f32 / slices as f32;
            m.edge_lines
                .push([radius * phi.cos(), radius * phi.sin(), 0.0]);
        }
    }
    m
}
