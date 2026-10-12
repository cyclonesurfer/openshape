//! Test 5: camera unit tests (no GPU needed).

use glam::{DVec3, Vec3, Vec4};
use opendatum_render::{Camera, Projection, ViewPreset};
use std::f64::consts::FRAC_PI_2;

fn box_corners(min: [f32; 3], max: [f32; 3]) -> Vec<Vec3> {
    (0..8)
        .map(|i| {
            Vec3::new(
                if i & 1 == 0 { min[0] } else { max[0] },
                if i & 2 == 0 { min[1] } else { max[1] },
                if i & 4 == 0 { min[2] } else { max[2] },
            )
        })
        .collect()
}

fn presets() -> [ViewPreset; 7] {
    [
        ViewPreset::Isometric,
        ViewPreset::Front,
        ViewPreset::Back,
        ViewPreset::Top,
        ViewPreset::Bottom,
        ViewPreset::Right,
        ViewPreset::Left,
    ]
}

fn projections() -> [Projection; 3] {
    [
        Projection::default(),
        Projection::Perspective { fov_y_rad: 1.2 },
        Projection::Orthographic,
    ]
}

#[test]
fn fit_puts_all_corners_inside_clip_volume() {
    // 1 mm to 10 m (model units are metres), off-origin, non-cubic boxes.
    let boxes = [
        ([0.0, 0.0, 0.0], [0.001, 0.001, 0.001]),
        ([10.0, -3.0, 2.0], [10.002, -2.999, 2.0005]),
        ([-5.0, -5.0, 0.0], [5.0, 5.0, 2.0]),
        ([-1.0, 100.0, -1.0], [1.0, 102.0, 9.0]),
    ];
    for (min, max) in boxes {
        for preset in presets() {
            for projection in projections() {
                for aspect in [1.0f32, 16.0 / 9.0, 0.5] {
                    let mut camera = Camera::from_preset(preset).with_projection(projection);
                    if aspect == 1.0 {
                        camera.fit(min, max);
                    } else {
                        camera.fit_with_aspect(min, max, aspect);
                    }
                    let vp = camera.view_proj(aspect);
                    for c in box_corners(min, max) {
                        let clip = vp * Vec4::new(c.x, c.y, c.z, 1.0);
                        assert!(clip.w > 0.0, "{preset:?} {projection:?}: w {}", clip.w);
                        let ndc = clip.truncate() / clip.w;
                        let ctx = format!("{preset:?} {projection:?} aspect {aspect} box {min:?}..{max:?}: {ndc:?}");
                        assert!(ndc.x.abs() <= 1.0 && ndc.y.abs() <= 1.0, "{ctx}");
                        // Reverse-Z: depth in (0, 1], never clipped by near.
                        assert!(ndc.z > 0.0 && ndc.z <= 1.0, "{ctx}");
                    }
                    // The target (box centre) projects to the screen centre.
                    // Tolerance: f32 matrices for a 2 mm box 10 m from the
                    // origin carry ~1e-6 m of error, i.e. ~5e-4 NDC (0.1 px).
                    let t = camera.target.as_vec3();
                    let ndc = vp.project_point3(t);
                    assert!(ndc.x.abs() < 2e-3 && ndc.y.abs() < 2e-3, "{ndc:?}");
                }
            }
        }
    }
}

#[test]
fn fit_frames_tightly() {
    // The fitted sphere should fill a reasonable part of the view, not a speck.
    let (min, max) = ([-0.5f32; 3], [0.5f32; 3]);
    for projection in projections() {
        let camera = Camera::front().with_projection(projection).fitted(min, max);
        let vp = camera.view_proj(1.0);
        let max_extent = box_corners(min, max)
            .into_iter()
            .map(|c| vp.project_point3(c).x.abs())
            .fold(0.0f32, f32::max);
        assert!(max_extent > 0.4, "{projection:?}: {max_extent}");
    }
}

#[test]
fn pitch_is_clamped() {
    let mut c = Camera::isometric();
    c.orbit(0.0, 1e6);
    assert_eq!(c.pitch, FRAC_PI_2);
    assert!((c.forward() - (-DVec3::Z)).length() < 1e-12);
    c.orbit(0.0, -1e7);
    assert_eq!(c.pitch, -FRAC_PI_2);
    c.orbit_radians(0.0, -10.0);
    assert_eq!(c.pitch, -FRAC_PI_2);
    c.orbit_radians(0.0, 3.0 * FRAC_PI_2);
    assert_eq!(c.pitch, FRAC_PI_2);

    // Even a hand-edited out-of-range pitch never flips the view: the screen
    // up vector keeps a non-negative world Z component.
    c.pitch = 2.5;
    let up = c.up_dir();
    assert!(up.z >= -1e-12, "{up:?}");
    assert!(c.view_matrix().is_finite());
    // Sweeping pitch through the full range keeps up.z >= 0 and up continuous.
    let mut prev: Option<DVec3> = None;
    for i in 0..=200 {
        c.pitch = -FRAC_PI_2 + std::f64::consts::PI * i as f64 / 200.0;
        let up = c.up_dir();
        assert!(up.z >= -1e-12);
        if let Some(p) = prev {
            assert!((up - p).length() < 0.05, "up jumps at step {i}");
        }
        prev = Some(up);
    }
}

#[test]
fn matrices_are_finite() {
    for preset in presets() {
        for projection in projections() {
            for distance in [1e-6, 1e-3, 1.0, 1e3, 1e6] {
                for aspect in [0.01f32, 1.0, 3.0, 0.0, f32::NAN] {
                    let mut c = Camera::from_preset(preset).with_projection(projection);
                    c.distance = distance;
                    c.target = DVec3::new(1.0, -2.0, 3.0);
                    assert!(c.view_matrix().is_finite());
                    assert!(
                        c.proj_matrix(aspect).is_finite(),
                        "{preset:?} {projection:?}"
                    );
                    assert!(c.view_proj(aspect).is_finite());
                    c.fit([0.0; 3], [1.0, 2.0, 3.0]);
                    assert!(c.view_proj(aspect).is_finite());
                }
            }
        }
    }
    // Degenerate bounds (a point) still give a usable camera.
    let c = Camera::isometric().fitted([1.0; 3], [1.0; 3]);
    assert!(c.distance.is_finite() && c.distance > 0.0);
    assert!(c.view_proj(1.0).is_finite());
}

#[test]
fn reverse_z_depth_ordering() {
    for projection in projections() {
        let c = Camera::front()
            .with_projection(projection)
            .fitted([-1.0; 3], [1.0; 3]);
        let vp = c.view_proj(1.0);
        // Front view looks along +Y: smaller y is nearer, so deeper (larger) z.
        let near = vp.project_point3(Vec3::new(0.0, -1.0, 0.0)).z;
        let mid = vp.project_point3(Vec3::ZERO).z;
        let far = vp.project_point3(Vec3::new(0.0, 1.0, 0.0)).z;
        assert!(
            near > mid && mid > far,
            "{projection:?}: {near} {mid} {far}"
        );
        assert!(far > 0.0 && near <= 1.0);
    }
}

#[test]
fn zoom_is_multiplicative() {
    let mut c = Camera::isometric();
    c.distance = 8.0;
    c.zoom(2.0);
    assert!((c.distance - 4.0).abs() < 1e-12);
    c.zoom(0.25);
    assert!((c.distance - 16.0).abs() < 1e-12);
    let d = c.distance;
    c.zoom(1.1);
    c.zoom(1.3);
    assert!((c.distance - d / (1.1 * 1.3)).abs() < 1e-9);
    // Invalid factors are ignored.
    let d = c.distance;
    for f in [0.0, -2.0, f64::NAN, f64::INFINITY] {
        c.zoom(f);
        assert_eq!(c.distance, d);
    }
    // Zoom keeps target and direction.
    let (t, fwd) = (c.target, c.forward());
    c.zoom(3.0);
    assert_eq!(c.target, t);
    assert_eq!(c.forward(), fwd);
}

#[test]
fn orbit_keeps_target_and_distance() {
    let mut c = Camera::isometric().fitted([0.0; 3], [2.0, 1.0, 1.0]);
    let (t, d) = (c.target, c.distance);
    let eye_dist = (c.eye() - t).length();
    c.orbit(37.0, -12.0);
    assert_eq!(c.target, t);
    assert_eq!(c.distance, d);
    assert!(((c.eye() - t).length() - eye_dist).abs() < 1e-9);
    // Dragging right turns the model right, i.e. the eye moves left (yaw decreases).
    let yaw = c.yaw;
    c.orbit(10.0, 0.0);
    assert!(c.yaw < yaw);
}

#[test]
fn pan_follows_the_cursor() {
    for projection in projections() {
        let (w, h) = (800.0f32, 600.0f32);
        let mut c = Camera::isometric()
            .with_projection(projection)
            .fitted([-1.0; 3], [1.0; 3]);
        let p = c.target.as_vec3();
        let to_px = |c: &Camera, p: Vec3| {
            let ndc = c.view_proj(w / h).project_point3(p);
            ((ndc.x * 0.5 + 0.5) * w, (0.5 - ndc.y * 0.5) * h)
        };
        let (x0, y0) = to_px(&c, p);
        c.pan(30.0, -20.0, h);
        let (x1, y1) = to_px(&c, p);
        assert!(
            (x1 - x0 - 30.0).abs() < 0.01,
            "{projection:?}: dx {}",
            x1 - x0
        );
        assert!(
            (y1 - y0 + 20.0).abs() < 0.01,
            "{projection:?}: dy {}",
            y1 - y0
        );
    }
}

#[test]
fn presets_keep_zoom_when_switching_views() {
    let mut c = Camera::isometric().fitted([0.0; 3], [1.0; 3]);
    let (t, d) = (c.target, c.distance);
    c.set_view(ViewPreset::Top);
    assert_eq!((c.target, c.distance), (t, d));
    assert!((c.forward() + DVec3::Z).length() < 1e-12);
}
