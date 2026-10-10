//! GPU acceptance tests: run headlessly (lavapipe in CI). They fail loudly
//! if no adapter is available.

mod common;

use common::*;
use glam::Mat4;
use openshape_render::{Camera, DisplayMode, MeshData, Projection, RgbaImage};
use std::f64::consts::FRAC_PI_2;
use std::time::Instant;

const W: u32 = 320;
const H: u32 = 240;

fn background(r: &mut openshape_render::OffscreenRenderer, w: u32, h: u32) -> RgbaImage {
    r.set_meshes(&[]).unwrap();
    r.render_rgba(&Camera::isometric(), w, h).unwrap()
}

/// Test 1: shaded unit cube, isometric, fitted.
#[test]
fn unit_cube_isometric() {
    let mut r = renderer();
    let bg = background(&mut r, W, H);

    let cube = unit_cube();
    assert_eq!(cube.positions.len(), 24);
    assert_eq!(cube.normals.len(), 24);
    assert_eq!(cube.triangle_count(), 12);
    assert_eq!(cube.edge_segment_count(), 12);
    let (min, max) = cube.bounds().unwrap();
    let mut camera = Camera::isometric();
    camera.fit(min, max);

    let img = r.render_meshes_rgba(&[cube], &camera, W, H).unwrap();
    maybe_dump("unit_cube_isometric", &img);
    assert_eq!(img.pixels.len(), (W * H * 4) as usize);

    let (cx, cy) = (W / 2, H / 2);
    assert_ne!(img.pixel(cx, cy), bg.pixel(cx, cy), "centre shows the cube");
    for (x, y) in [(0, 0), (W - 1, 0), (0, H - 1), (W - 1, H - 1)] {
        assert_eq!(
            img.pixel(x, y),
            bg.pixel(x, y),
            "corner ({x},{y}) is background"
        );
    }

    // Sample the centre of each visible face.
    let top = img.pixel_at(project(&camera, [0.0, 0.0, 0.5], W, H));
    let front = img.pixel_at(project(&camera, [0.0, -0.5, 0.0], W, H));
    let right = img.pixel_at(project(&camera, [0.5, 0.0, 0.0], W, H));
    eprintln!("top {top:?} front {front:?} right {right:?}");
    for (name, face, (x, y)) in [
        ("top", top, project(&camera, [0.0, 0.0, 0.5], W, H)),
        ("front", front, project(&camera, [0.0, -0.5, 0.0], W, H)),
        ("right", right, project(&camera, [0.5, 0.0, 0.0], W, H)),
    ] {
        assert!(
            rgb_diff(face, bg.pixel(x, y)) > 30,
            "{name} face differs from background"
        );
    }
    assert!(rgb_diff(top, front) > 20, "top vs front: {top:?} {front:?}");
    assert!(rgb_diff(top, right) > 20, "top vs right: {top:?} {right:?}");
    assert!(
        rgb_diff(front, right) > 20,
        "front vs right: {front:?} {right:?}"
    );
    // The lighting design: top brightest, then front, then right.
    assert!(luma(top) > luma(front) && luma(front) > luma(right));

    // A visible edge (top/front) is dark: some pixel near its projected
    // midpoint is much darker than both faces.
    let (ex, ey) = project(&camera, [0.0, -0.5, 0.5], W, H);
    let darkest = (ey - 1..=ey + 1)
        .flat_map(|y| (ex - 1..=ex + 1).map(move |x| (x, y)))
        .map(|(x, y)| luma(img.pixel(x, y)))
        .fold(f32::MAX, f32::min);
    assert!(
        darkest < 0.6 * luma(right),
        "edge is drawn: darkest {darkest}"
    );
}

/// Display modes: shaded has no edge lines, wireframe has no faces.
#[test]
fn display_modes() {
    let mut r = renderer();
    let bg = background(&mut r, W, H);
    let cube = unit_cube();
    let camera = Camera::isometric().fitted([-0.5; 3], [0.5; 3]);
    r.set_meshes(&[cube]).unwrap();

    let edge_px = project(&camera, [0.0, -0.5, 0.5], W, H);
    // Off the symmetry lines: in an isometric view every face centre lies on
    // the projection of a hidden edge, which wireframe mode draws.
    let face_px = project(&camera, [0.2, -0.1, 0.5], W, H);

    r.scene_mut().set_display_mode(DisplayMode::ShadedWithEdges);
    let with_edges = r.render_rgba(&camera, W, H).unwrap();
    r.scene_mut().set_display_mode(DisplayMode::Shaded);
    let shaded = r.render_rgba(&camera, W, H).unwrap();
    r.scene_mut().set_display_mode(DisplayMode::Wireframe);
    let wire = r.render_rgba(&camera, W, H).unwrap();
    maybe_dump("display_shaded", &shaded);
    maybe_dump("display_wireframe", &wire);

    // Shaded: the edge location shows the face colour, not a dark line.
    let near_edge = |img: &RgbaImage| {
        let (ex, ey) = edge_px;
        (ey - 1..=ey + 1)
            .flat_map(|y| (ex - 1..=ex + 1).map(move |x| (x, y)))
            .map(|(x, y)| luma(img.pixel(x, y)))
            .fold(f32::MAX, f32::min)
    };
    assert!(near_edge(&shaded) > near_edge(&with_edges) + 40.0);
    // Faces identical away from edges.
    assert_eq!(shaded.pixel_at(face_px), with_edges.pixel_at(face_px));
    // Hidden edges stay hidden: the top face centre lies on the projection
    // of the hidden back vertical edge.
    let hidden_px = project(&camera, [0.0, 0.0, 0.5], W, H);
    let around_hidden = || {
        let (hx, hy) = hidden_px;
        (hy - 1..=hy + 1).flat_map(move |y| (hx - 1..=hx + 1).map(move |x| (x, y)))
    };
    assert!(around_hidden().all(|(x, y)| shaded.pixel(x, y) == with_edges.pixel(x, y)));
    // Wireframe draws it: some pixel next to the projected point (the line
    // need not cover that exact pixel's centre) differs from the background.
    let strongest = around_hidden()
        .map(|(x, y)| rgb_diff(wire.pixel(x, y), bg.pixel(x, y)))
        .max()
        .unwrap();
    assert!(strongest > 30, "hidden edge in wireframe: diff {strongest}");
    // Wireframe: face centre is background, edges still drawn.
    assert_eq!(wire.pixel_at(face_px), bg.pixel_at(face_px));
    assert!(near_edge(&wire) < luma(bg.pixel_at(edge_px)) - 40.0);
}

/// Test 2: orbiting changes the image; orbiting back restores it exactly.
#[test]
fn orbit_changes_and_restores_image() {
    let mut r = renderer();
    // A unit cube is symmetric under a 90 degree turn about Z (and the
    // lighting is camera- and Z-relative), so use an asymmetric scene: a
    // flat box plus a smaller cube on one end.
    let meshes = [
        MeshData::cuboid([-1.0, -0.5, -0.3], [1.0, 0.5, 0.3], MeshData::DEFAULT_COLOR),
        MeshData::cuboid([0.4, -0.3, 0.3], [0.9, 0.2, 0.8], [0.8, 0.5, 0.4, 1.0]),
    ];
    r.set_meshes(&meshes).unwrap();
    let (min, max) = r.scene().bounds().unwrap();
    let mut camera = Camera::isometric();
    camera.fit(min, max);

    let before = r.render_rgba(&camera, W, H).unwrap();
    let center = (W / 4, H / 4, 3 * W / 4, 3 * H / 4);
    let area = ((W / 2) * (H / 2)) as usize;

    camera.orbit_radians(FRAC_PI_2, 0.0);
    let turned = r.render_rgba(&camera, W, H).unwrap();
    maybe_dump("orbit_before", &before);
    maybe_dump("orbit_turned", &turned);
    let changed = count_different(&before, &turned, center, 6);
    assert!(
        changed > area / 20,
        "90 degree orbit changed only {changed} of {area} centre pixels"
    );

    camera.orbit_radians(-FRAC_PI_2, 0.0);
    let back = r.render_rgba(&camera, W, H).unwrap();
    assert!(back == before, "orbiting back restores an identical image");

    // Same with pixel-based orbiting, including pitch.
    camera.orbit(150.0, -40.0);
    let moved = r.render_rgba(&camera, W, H).unwrap();
    assert!(count_different(&before, &moved, center, 6) > area / 20);
    camera.orbit(-150.0, 40.0);
    let back = r.render_rgba(&camera, W, H).unwrap();
    assert!(back == before, "pixel orbit round trip restores the image");
}

/// Test 3: orthographic and perspective both render the cube.
#[test]
fn orthographic_and_perspective() {
    let mut r = renderer();
    let bg = background(&mut r, W, H);
    r.set_meshes(&[unit_cube()]).unwrap();

    let mut images = Vec::new();
    for projection in [Projection::default(), Projection::Orthographic] {
        let camera = Camera::isometric()
            .with_projection(projection)
            .fitted([-0.5; 3], [0.5; 3]);
        let img = r.render_rgba(&camera, W, H).unwrap();
        maybe_dump(&format!("projection_{projection:?}"), &img);
        assert_ne!(
            img.pixel(W / 2, H / 2),
            bg.pixel(W / 2, H / 2),
            "{projection:?}"
        );
        for (x, y) in [(0, 0), (W - 1, H - 1)] {
            assert_eq!(img.pixel(x, y), bg.pixel(x, y), "{projection:?} corner");
        }
        let faces: Vec<[u8; 4]> = [[0.0, 0.0, 0.5], [0.0, -0.5, 0.0], [0.5, 0.0, 0.0]]
            .iter()
            .map(|p| img.pixel_at(project(&camera, *p, W, H)))
            .collect();
        for i in 0..3 {
            assert!(
                rgb_diff(faces[i], faces[(i + 1) % 3]) > 20,
                "{projection:?} {faces:?}"
            );
        }
        // The whole cube is in view: count cube pixels.
        let covered = count_different(&img, &bg, (0, 0, W, H), 6);
        assert!(
            covered > (W * H / 10) as usize,
            "{projection:?}: {covered} pixels"
        );
        images.push(img);
    }
    assert!(images[0] != images[1], "projections differ");

    // In an orthographic front view the cube is an exact centred square.
    let camera = Camera::front()
        .with_projection(Projection::Orthographic)
        .fitted([-0.5; 3], [0.5; 3]);
    let img = r.render_rgba(&camera, W, H).unwrap();
    let (l, rr) = row_extent(&img, &bg, H / 2);
    let (t, b) = col_extent(&img, &bg, W / 2);
    let (w, h) = (rr - l, b - t);
    assert!(w.abs_diff(h) <= 2, "square: {w}x{h}");
}

/// Test 4: readback of sizes whose row pitch needs padding.
#[test]
fn readback_handles_unaligned_widths() {
    let mut r = renderer();
    for (w, h) in [(333, 217), (1, 1), (65, 3), (257, 64)] {
        let bg = background(&mut r, w, h);
        assert_eq!((bg.width, bg.height), (w, h));
        assert_eq!(bg.pixels.len(), (w * h * 4) as usize);
        // The background is a vertical gradient: every row is uniform and
        // brightness grows towards the bottom. Row-padding mistakes shear
        // the gradient and break both properties.
        for y in 0..h {
            let first = bg.pixel(0, y);
            for x in 0..w {
                assert!(
                    rgb_diff(bg.pixel(x, y), first) <= 3,
                    "{w}x{h} row {y} not uniform"
                );
            }
            assert_eq!(first[3], 255);
        }
        if h > 1 {
            assert!(luma(bg.pixel(0, h - 1)) > luma(bg.pixel(0, 0)));
        }
    }

    // A front-view cube must be left-right symmetric on every row.
    let (w, h) = (333, 217);
    let bg = background(&mut r, w, h);
    let camera = Camera::front().fitted([-0.5; 3], [0.5; 3]);
    let img = r.render_meshes_rgba(&[unit_cube()], &camera, w, h).unwrap();
    maybe_dump("readback_333x217", &img);
    assert_ne!(img.pixel(w / 2, h / 2), bg.pixel(w / 2, h / 2));
    assert_eq!(img.pixel(0, 0), bg.pixel(0, 0));
    assert_eq!(img.pixel(w - 1, h - 1), bg.pixel(w - 1, h - 1));
    for y in [h / 3, h / 2, 2 * h / 3] {
        let (l, rr) = row_extent(&img, &bg, y);
        let mid2 = l + rr; // twice the midpoint
        assert!(
            mid2.abs_diff(w - 1) <= 2,
            "row {y}: extent {l}..{rr} not centred"
        );
    }
}

/// Larger scene with concave edges and smooth surfaces (also a visual check
/// via OPENSHAPE_RENDER_TEST_OUT).
#[test]
fn mixed_scene_renders() {
    let mut r = renderer();
    let meshes = vec![
        MeshData::cuboid([0.0, 0.0, 0.0], [2.0, 1.0, 0.3], MeshData::DEFAULT_COLOR),
        MeshData::cuboid([0.0, 0.0, 0.3], [0.3, 1.0, 1.5], MeshData::DEFAULT_COLOR),
        cylinder([1.3, 0.5, 0.65], 0.3, 0.7, 48),
        uv_sphere(0.25, 24, 16)
            .transformed(&Mat4::from_translation(glam::Vec3::new(1.3, 0.5, 1.3))),
        MeshData {
            color: [0.3, 0.7, 0.4, 0.45],
            ..MeshData::cuboid([1.8, -0.4, 0.0], [2.4, 0.2, 0.6], MeshData::DEFAULT_COLOR)
        },
    ];
    r.set_meshes(&meshes).unwrap();
    let stats = r.scene().stats();
    assert_eq!(stats.meshes, 5);
    let (min, max) = r.scene().bounds().unwrap();
    for (name, camera) in [
        ("mixed_iso", Camera::isometric()),
        (
            "mixed_iso_ortho",
            Camera::isometric().with_projection(Projection::Orthographic),
        ),
        ("mixed_top", Camera::top()),
    ] {
        let mut camera = camera;
        camera.fit_with_aspect(min, max, 640.0 / 480.0);
        let img = r.render_rgba(&camera, 640, 480).unwrap();
        maybe_dump(name, &img);
        assert_ne!(img.pixel(320, 240), img.pixel(0, 0));
    }
    // Scale invariance: the same scene at 1 mm and 10 km looks the same.
    let base = {
        let camera = Camera::isometric().fitted(min, max);
        r.render_rgba(&camera, 200, 150).unwrap()
    };
    for scale in [1e-3f32, 1e4] {
        let scaled: Vec<MeshData> = meshes
            .iter()
            .map(|m| m.transformed(&Mat4::from_scale(glam::Vec3::splat(scale))))
            .collect();
        r.set_meshes(&scaled).unwrap();
        let (min, max) = r.scene().bounds().unwrap();
        let camera = Camera::isometric().fitted(min, max);
        let img = r.render_rgba(&camera, 200, 150).unwrap();
        maybe_dump(&format!("mixed_scale_{scale}"), &img);
        let differing = count_different(&base, &img, (0, 0, 200, 150), 30);
        assert!(differing < 200, "scale {scale}: {differing} pixels differ");
    }
}

/// Performance smoke test: 1,000 instances of a ~200-triangle mesh at
/// 640x480. Prints timings (run with `--nocapture`); asserts nothing about
/// speed since software rasterisers vary widely.
#[test]
fn thousand_instances_timing() {
    let mut r = renderer();
    let sphere = uv_sphere(1.0, 12, 9);
    assert_eq!(sphere.triangle_count(), 192);
    let mut meshes = Vec::with_capacity(1000);
    for i in 0..10 {
        for j in 0..10 {
            for k in 0..10 {
                let t = glam::Vec3::new(i as f32, j as f32, k as f32) * 3.0;
                let mut m = sphere.transformed(&Mat4::from_translation(t));
                m.color = [
                    0.4 + 0.06 * i as f32,
                    0.5 + 0.04 * j as f32,
                    0.8 - 0.05 * k as f32,
                    1.0,
                ];
                meshes.push(m);
            }
        }
    }
    let t0 = Instant::now();
    r.set_meshes(&meshes).unwrap();
    let upload = t0.elapsed();
    let stats = r.scene().stats();
    assert_eq!(stats.meshes, 1000);
    assert_eq!(stats.triangles, 192_000);

    let (min, max) = r.scene().bounds().unwrap();
    let camera = Camera::isometric().fitted(min, max);
    let t1 = Instant::now();
    let first = r.render_rgba(&camera, 640, 480).unwrap();
    let first_time = t1.elapsed();
    maybe_dump("thousand_instances", &first);

    let frames = 5;
    let t2 = Instant::now();
    for _ in 0..frames {
        r.render_rgba(&camera, 640, 480).unwrap();
    }
    let per_frame = t2.elapsed() / frames;
    let info = r.adapter_info();
    println!(
        "1000 x 192-triangle instances @ 640x480 on {} ({:?}): upload {:.1} ms, \
         first frame {:.1} ms, steady {:.1} ms/frame incl. readback ({} GPU batches)",
        info.name,
        info.backend,
        upload.as_secs_f64() * 1e3,
        first_time.as_secs_f64() * 1e3,
        per_frame.as_secs_f64() * 1e3,
        stats.batches,
    );
    assert_ne!(first.pixel(320, 240), first.pixel(0, 0));
}

/// Invalid sizes are reported, not panicked on.
#[test]
fn zero_size_is_an_error() {
    let mut r = renderer();
    assert!(r.render_rgba(&Camera::default(), 0, 10).is_err());
    assert!(r.render_rgba(&Camera::default(), 10, 0).is_err());
}

/// Malformed meshes are sanitised instead of crashing the GPU.
#[test]
fn malformed_meshes_are_tolerated() {
    let mut r = renderer();
    let mut bad = unit_cube();
    bad.normals.clear(); // recomputed
    bad.indices.push(9999); // incomplete + out of range
    bad.indices.extend_from_slice(&[0, 1, 5000]); // out of range triangle
    bad.edge_lines.push([0.0; 3]); // odd count
    r.set_meshes(&[bad, MeshData::default()]).unwrap();
    let stats = r.scene().stats();
    assert_eq!(stats.meshes, 1);
    assert_eq!(stats.triangles, 12);
    assert_eq!(stats.edge_segments, 12);
    let img = r
        .render_rgba(&Camera::isometric().fitted([-0.5; 3], [0.5; 3]), 64, 64)
        .unwrap();
    assert_ne!(img.pixel(32, 32), img.pixel(0, 0));
}

trait PixelAt {
    fn pixel_at(&self, p: (u32, u32)) -> [u8; 4];
}

impl PixelAt for RgbaImage {
    fn pixel_at(&self, (x, y): (u32, u32)) -> [u8; 4] {
        self.pixel(x, y)
    }
}

/// First and last x on row `y` that differ from the background.
fn row_extent(img: &RgbaImage, bg: &RgbaImage, y: u32) -> (u32, u32) {
    let xs: Vec<u32> = (0..img.width)
        .filter(|&x| rgb_diff(img.pixel(x, y), bg.pixel(x, y)) > 6)
        .collect();
    (*xs.first().expect("row hits"), *xs.last().unwrap())
}

/// First and last y in column `x` that differ from the background.
fn col_extent(img: &RgbaImage, bg: &RgbaImage, x: u32) -> (u32, u32) {
    let ys: Vec<u32> = (0..img.height)
        .filter(|&y| rgb_diff(img.pixel(x, y), bg.pixel(x, y)) > 6)
        .collect();
    (*ys.first().expect("column hits"), *ys.last().unwrap())
}

/// Concave edges (the inner corner of an L bracket) stay visible down to
/// shallow viewing angles, where the adjacent faces approach edge-on.
#[test]
fn concave_edges_visible_at_oblique_angles() {
    let mut r = renderer();
    let meshes = [
        MeshData::cuboid([0.0, 0.0, 0.0], [2.0, 1.0, 0.3], MeshData::DEFAULT_COLOR),
        MeshData::cuboid([0.0, 0.0, 0.3], [0.3, 1.0, 1.5], MeshData::DEFAULT_COLOR),
    ];
    r.set_meshes(&meshes).unwrap();
    let (min, max) = r.scene().bounds().unwrap();
    let (w, h) = (480, 360);
    for pitch_deg in [35.0f64, 15.0, 8.0, 4.0] {
        for projection in [Projection::default(), Projection::Orthographic] {
            let mut camera = Camera::front().with_projection(projection);
            camera.yaw = 30f64.to_radians();
            camera.pitch = pitch_deg.to_radians();
            camera.fit_with_aspect(min, max, w as f32 / h as f32);
            camera.zoom(1.6);
            let img = r.render_rgba(&camera, w, h).unwrap();
            maybe_dump(&format!("concave_{pitch_deg}_{projection:?}"), &img);
            // Reference: the shade of the base's top face away from edges.
            let face = luma(img.pixel_at(project(&camera, [1.2, 0.6, 0.3], w, h)));
            let samples = 21;
            let dark = (1..samples)
                .filter(|i| {
                    let y = *i as f32 / samples as f32;
                    let (ex, ey) = project(&camera, [0.3, y, 0.3], w, h);
                    let darkest = (ey.saturating_sub(1)..=(ey + 1).min(h - 1))
                        .flat_map(|yy| {
                            (ex.saturating_sub(1)..=(ex + 1).min(w - 1)).map(move |xx| (xx, yy))
                        })
                        .map(|(xx, yy)| luma(img.pixel(xx, yy)))
                        .fold(f32::MAX, f32::min);
                    darkest < 0.6 * face
                })
                .count();
            eprintln!(
                "pitch {pitch_deg:>4} {projection:?}: {dark}/{} concave edge samples dark",
                samples - 1
            );
            assert!(
                dark >= 18,
                "pitch {pitch_deg} {projection:?}: only {dark} dark"
            );
        }
    }
}

/// The UI path: a SceneRenderer created on someone else's device, drawing
/// into a non-sRGB BGRA texture (as egui-wgpu commonly uses), matches the
/// offscreen sRGB render.
#[test]
fn shared_device_non_srgb_target_matches() {
    use openshape_render::{read_texture_rgba8, wgpu, SceneRenderer};
    let mut r = renderer();
    let camera = Camera::isometric().fitted([-0.5; 3], [0.5; 3]);
    let reference = r.render_meshes_rgba(&[unit_cube()], &camera, W, H).unwrap();
    // Where edge lines land: pixels that change when edges are drawn.
    r.scene_mut().set_display_mode(DisplayMode::Shaded);
    let no_edges = r.render_rgba(&camera, W, H).unwrap();
    r.scene_mut().set_display_mode(DisplayMode::ShadedWithEdges);
    let on_line = |x: u32, y: u32| rgb_diff(reference.pixel(x, y), no_edges.pixel(x, y)) > 6;

    let device = r.device().clone();
    let queue = r.queue().clone();
    let format = wgpu::TextureFormat::Bgra8Unorm;
    let mut scene = SceneRenderer::new(&device, format);
    scene.set_meshes(&device, &[unit_cube()]);
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("test viewport"),
        size: wgpu::Extent3d {
            width: W,
            height: H,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let view = texture.create_view(&Default::default());
    scene.render(&device, &queue, &camera, &view, (W, H));
    let img = read_texture_rgba8(&device, &queue, &texture).unwrap();
    maybe_dump("shared_bgra", &img);
    // Everything matches to rounding except antialiased pixels at edge
    // lines: a non-sRGB target resolves MSAA on encoded values, not linear.
    // Every silhouette of the cube carries an edge line, so this covers all
    // partially covered pixels.
    let mut unexplained = 0;
    for y in 0..H {
        for x in 0..W {
            if rgb_diff(reference.pixel(x, y), img.pixel(x, y)) <= 6 {
                continue;
            }
            let near_line = (y.saturating_sub(2)..(y + 3).min(H))
                .flat_map(|yy| (x.saturating_sub(2)..(x + 3).min(W)).map(move |xx| (xx, yy)))
                .any(|(xx, yy)| on_line(xx, yy));
            if !near_line {
                unexplained += 1;
            }
        }
    }
    assert_eq!(unexplained, 0, "pixels differ away from edge lines");
    assert_eq!(img.pixel(W / 2, H / 4), reference.pixel(W / 2, H / 4));

    // Resizing recreates attachments.
    let small = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("test viewport small"),
        size: wgpu::Extent3d {
            width: 100,
            height: 70,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    scene.render(
        &device,
        &queue,
        &camera,
        &small.create_view(&Default::default()),
        (100, 70),
    );
    let img = read_texture_rgba8(&device, &queue, &small).unwrap();
    assert_ne!(img.pixel(50, 35), img.pixel(0, 0));
}

#[test]
fn renderer_types_are_send_and_sync() {
    fn check<T: Send + Sync + 'static>() {}
    check::<openshape_render::SceneRenderer>();
    check::<openshape_render::OffscreenRenderer>();
    check::<Camera>();
    check::<MeshData>();
}

/// Edge lines are screen-space quads of the configured width, so the inked
/// area grows in proportion to the width on every driver.
#[test]
fn edge_width_is_respected() {
    let mut r = renderer();
    let bg = background(&mut r, W, H);
    let camera = Camera::isometric().fitted([-0.5; 3], [0.5; 3]);
    r.set_meshes(&[unit_cube()]).unwrap();
    r.scene_mut().set_display_mode(DisplayMode::Wireframe);
    let ink = |r: &mut openshape_render::OffscreenRenderer, width: f32| {
        r.scene_mut().set_edge_width(width);
        let img = r.render_rgba(&camera, W, H).unwrap();
        maybe_dump(&format!("edge_width_{width}"), &img);
        // Coverage-weighted: each pixel counts by how far it moved from the
        // background towards the edge colour.
        let mut total = 0.0;
        for y in 0..H {
            for x in 0..W {
                let b = luma(bg.pixel(x, y));
                total += ((b - luma(img.pixel(x, y))) / b).clamp(0.0, 1.0) as f64;
            }
        }
        total
    };
    let thin = ink(&mut r, 1.5);
    let thick = ink(&mut r, 4.5);
    let ratio = thick / thin;
    assert!(
        thin > 0.0 && (2.4..3.6).contains(&ratio),
        "ink ratio {ratio}"
    );
    assert_eq!(r.scene().edge_width(), 4.5);
    r.scene_mut().set_edge_width(f32::NAN);
    assert_eq!(r.scene().edge_width(), 4.5);
    r.scene_mut().set_edge_width(100.0);
    assert_eq!(r.scene().edge_width(), 16.0);
}
