//! View-cube geometry: target angles, picking from the standard views and
//! from the isometric view, and points off the cube. No GPU needed.

use egui::{Pos2, Rect, Vec2};
use glam::DVec3;
use openshape_render::{Camera, ViewPreset};
use openshape_ui_kit::view_cube::{self, CubeCorner, CubeEdge, CubeFace, CubeTarget, CENTER_HALF};
use proptest::prelude::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

fn rect() -> Rect {
    Rect::from_min_size(Pos2::new(500.0, 40.0), Vec2::splat(120.0))
}

fn camera_at((yaw, pitch): (f64, f64)) -> Camera {
    Camera {
        yaw,
        pitch,
        ..Camera::default()
    }
}

fn dir(t: CubeTarget) -> DVec3 {
    let d = t.direction();
    DVec3::new(d[0] as f64, d[1] as f64, d[2] as f64)
}

/// Centre of the cell of `target` on `face` (which must be one of the
/// target's faces), in cube coordinates.
fn cell_center(face: CubeFace, target: CubeTarget) -> [f64; 3] {
    let n = {
        let n = face.normal();
        DVec3::new(n[0] as f64, n[1] as f64, n[2] as f64)
    };
    let d = dir(target);
    assert!(
        (d.dot(n) - 1.0).abs() < 1e-12,
        "{face:?} is not a face of {target:?}"
    );
    // Edge and corner bands span CENTER_HALF..1, so their middle is here.
    let band_mid = (1.0 + CENTER_HALF) * 0.5;
    let c = n + (d - n) * band_mid;
    [c.x, c.y, c.z]
}

fn faces_of(t: CubeTarget) -> Vec<CubeFace> {
    match t {
        CubeTarget::Face(f) => vec![f],
        CubeTarget::Edge(e) => e.faces().to_vec(),
        CubeTarget::Corner(c) => c.faces().to_vec(),
    }
}

#[test]
fn there_are_26_distinct_targets_that_round_trip_through_directions() {
    let all = CubeTarget::all();
    assert_eq!(all.len(), 26);
    let dirs: std::collections::HashSet<[i8; 3]> = all.iter().map(|t| t.direction()).collect();
    assert_eq!(dirs.len(), 26);
    for t in all {
        assert_eq!(CubeTarget::from_direction(t.direction()), Some(t));
    }
    assert_eq!(CubeTarget::from_direction([0, 0, 0]), None);
    assert_eq!(CubeTarget::from_direction([2, 0, 0]), None);
}

#[test]
fn face_angles_match_the_view_presets_exactly() {
    let pairs = [
        (CubeFace::Front, ViewPreset::Front),
        (CubeFace::Back, ViewPreset::Back),
        (CubeFace::Top, ViewPreset::Top),
        (CubeFace::Bottom, ViewPreset::Bottom),
        (CubeFace::Left, ViewPreset::Left),
        (CubeFace::Right, ViewPreset::Right),
    ];
    for (face, preset) in pairs {
        assert_eq!(CubeTarget::Face(face).yaw_pitch(), preset.yaw_pitch());
    }
    assert_eq!(
        CubeTarget::HOME.yaw_pitch(),
        ViewPreset::Isometric.yaw_pitch()
    );
    assert_eq!(
        CubeTarget::Corner(CubeCorner::FrontTopRight),
        CubeTarget::HOME
    );
}

#[test]
fn edge_and_corner_angles_match_the_expected_table() {
    let iso = (1.0f64 / 2.0f64.sqrt()).atan();
    let q = FRAC_PI_4;
    let expected: [(CubeTarget, (f64, f64)); 20] = [
        (CubeTarget::Edge(CubeEdge::FrontTop), (0.0, q)),
        (CubeTarget::Edge(CubeEdge::FrontBottom), (0.0, -q)),
        (CubeTarget::Edge(CubeEdge::FrontLeft), (-q, 0.0)),
        (CubeTarget::Edge(CubeEdge::FrontRight), (q, 0.0)),
        (CubeTarget::Edge(CubeEdge::BackTop), (PI, q)),
        (CubeTarget::Edge(CubeEdge::BackBottom), (PI, -q)),
        (CubeTarget::Edge(CubeEdge::BackLeft), (-3.0 * q, 0.0)),
        (CubeTarget::Edge(CubeEdge::BackRight), (3.0 * q, 0.0)),
        (CubeTarget::Edge(CubeEdge::TopLeft), (-FRAC_PI_2, q)),
        (CubeTarget::Edge(CubeEdge::TopRight), (FRAC_PI_2, q)),
        (CubeTarget::Edge(CubeEdge::BottomLeft), (-FRAC_PI_2, -q)),
        (CubeTarget::Edge(CubeEdge::BottomRight), (FRAC_PI_2, -q)),
        (CubeTarget::Corner(CubeCorner::FrontTopLeft), (-q, iso)),
        (CubeTarget::Corner(CubeCorner::FrontTopRight), (q, iso)),
        (CubeTarget::Corner(CubeCorner::FrontBottomLeft), (-q, -iso)),
        (CubeTarget::Corner(CubeCorner::FrontBottomRight), (q, -iso)),
        (CubeTarget::Corner(CubeCorner::BackTopLeft), (-3.0 * q, iso)),
        (CubeTarget::Corner(CubeCorner::BackTopRight), (3.0 * q, iso)),
        (
            CubeTarget::Corner(CubeCorner::BackBottomLeft),
            (-3.0 * q, -iso),
        ),
        (
            CubeTarget::Corner(CubeCorner::BackBottomRight),
            (3.0 * q, -iso),
        ),
    ];
    for (target, (yaw, pitch)) in expected {
        let (y, p) = target.yaw_pitch();
        assert!(
            (y - yaw).abs() < 1e-12 && (p - pitch).abs() < 1e-12,
            "{target:?}: got ({y}, {p}), expected ({yaw}, {pitch})"
        );
    }
}

#[test]
fn every_target_angle_looks_along_its_direction() {
    for t in CubeTarget::all() {
        let back = camera_at(t.yaw_pitch()).back_dir();
        let d = dir(t).normalize();
        assert!((back - d).length() < 1e-12, "{t:?}: {back:?} vs {d:?}");
    }
}

#[test]
fn picking_the_centre_from_each_face_view_returns_that_face() {
    for face in CubeFace::ALL {
        let camera = Camera::from_preset(face.view_preset());
        assert_eq!(
            view_cube::pick(&camera, rect(), rect().center()),
            Some(CubeTarget::Face(face)),
            "{face:?}"
        );
    }
}

#[test]
fn picking_the_centre_from_each_target_view_returns_that_target() {
    // Looking along an edge or a corner puts it at the centre of the cube.
    for t in CubeTarget::all() {
        let camera = camera_at(t.yaw_pitch());
        assert_eq!(view_cube::pick(&camera, rect(), rect().center()), Some(t));
    }
}

#[test]
fn from_isometric_every_visible_cell_picks_its_target() {
    let camera = Camera::isometric();
    let visible = [CubeFace::Front, CubeFace::Right, CubeFace::Top];
    let mut checked = 0;
    for t in CubeTarget::all() {
        for face in faces_of(t) {
            if !visible.contains(&face) {
                continue;
            }
            let pos = view_cube::project(&camera, rect(), cell_center(face, t));
            assert_eq!(
                view_cube::pick(&camera, rect(), pos),
                Some(t),
                "cell of {t:?} on {face:?} at {pos:?}"
            );
            checked += 1;
        }
    }
    // 3 visible faces x 9 cells each.
    assert_eq!(checked, 27);
}

#[test]
fn hidden_faces_are_never_picked_from_isometric() {
    let camera = Camera::isometric();
    let r = rect();
    let hidden = [CubeFace::Back, CubeFace::Left, CubeFace::Bottom];
    for i in 0..=60 {
        for j in 0..=60 {
            let pos = Pos2::new(
                r.left() + r.width() * i as f32 / 60.0,
                r.top() + r.height() * j as f32 / 60.0,
            );
            if let Some(t) = view_cube::pick(&camera, r, pos) {
                assert!(
                    faces_of(t).iter().any(|f| !hidden.contains(f)),
                    "{t:?} picked at {pos:?}"
                );
            }
        }
    }
}

#[test]
fn points_off_the_cube_pick_nothing() {
    let r = rect();
    for camera in [Camera::isometric(), Camera::front(), Camera::top()] {
        for pos in [
            r.left_top() + Vec2::splat(2.0),
            r.right_bottom() - Vec2::splat(2.0),
            r.right_top() + Vec2::new(-2.0, 2.0),
            Pos2::new(0.0, 0.0),
            r.center() + Vec2::new(0.0, r.height()),
        ] {
            assert_eq!(view_cube::pick(&camera, r, pos), None, "{pos:?}");
        }
    }
    // A degenerate rectangle has no cube, and a non-finite pointer hits
    // nothing.
    let empty = Rect::from_min_size(r.min, Vec2::ZERO);
    assert_eq!(
        view_cube::pick(&Camera::isometric(), empty, empty.min),
        None
    );
    let nan = Pos2::new(f32::NAN, r.center().y);
    assert_eq!(view_cube::pick(&Camera::isometric(), r, nan), None);
}

#[test]
fn the_cube_stays_inside_its_rectangle_at_any_orientation() {
    let r = rect();
    for k in 0..72 {
        let camera = camera_at((k as f64 * 0.37, (k as f64 * 0.23).sin() * FRAC_PI_2));
        for x in [-1.0, 1.0] {
            for y in [-1.0, 1.0] {
                for z in [-1.0, 1.0] {
                    let p = view_cube::project(&camera, r, [x, y, z]);
                    assert!(r.contains(p), "{p:?} outside {r:?}");
                }
            }
        }
    }
}

proptest! {
    /// For any camera, the projected centre of every face turned clearly
    /// towards the viewer picks that face, and the screen centre picks one
    /// of the faces that face the viewer.
    #[test]
    fn visible_face_centres_pick_their_face(
        yaw in -PI..PI,
        pitch in -FRAC_PI_2..FRAC_PI_2,
    ) {
        let camera = camera_at((yaw, pitch));
        let back = camera.back_dir();
        for face in CubeFace::ALL {
            let n = face.normal();
            let n = DVec3::new(n[0] as f64, n[1] as f64, n[2] as f64);
            if n.dot(back) > 0.2 {
                let pos = view_cube::project(&camera, rect(), [n.x, n.y, n.z]);
                prop_assert_eq!(
                    view_cube::pick(&camera, rect(), pos),
                    Some(CubeTarget::Face(face))
                );
            }
        }
        let centre = view_cube::pick(&camera, rect(), rect().center());
        prop_assert!(centre.is_some());
        let t = centre.unwrap();
        prop_assert!(dir(t).dot(back) > 0.0, "{:?} faces away", t);
    }
}
