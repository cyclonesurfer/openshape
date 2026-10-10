//! Orbit (turntable) camera with scene-adaptive clip planes.
//!
//! # Conventions
//!
//! - World space is right-handed with **Z up** (the usual CAD convention).
//! - The camera orbits [`Camera::target`] at [`Camera::distance`]. `yaw` turns
//!   around world Z and `pitch` tilts above (positive) or below the XY plane.
//!   `yaw = 0, pitch = 0` is the front view: the eye is on the -Y side looking
//!   along +Y, with +X to the right and +Z up.
//! - The screen basis is built directly from yaw and pitch (the right vector
//!   is always horizontal), so straight top and bottom views are well defined.
//!   Pitch is clamped to `[-pi/2, pi/2]`, so the view never flips upside down.
//! - Matrices follow wgpu conventions (clip-space depth in `0..=1`) with
//!   **reverse-Z**: the near plane maps to depth 1 and far to depth 0, so the
//!   depth test is `Greater` and the depth buffer is cleared to 0. Perspective
//!   projections use an infinite far plane; together with a 32-bit float depth
//!   buffer this gives nearly uniform relative depth precision over the whole
//!   range. The near plane comes from the scene bounds (see
//!   [`Camera::clip_planes`]), so models from 1 mm to 10 m look the same.
//! - All state is `f64`; matrices are handed out as `f32` ([`glam::Mat4`]).

use glam::{DMat4, DVec3, DVec4, Mat4};
use serde::{Deserialize, Serialize};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

/// Vertical field of view of the default perspective projection (30 degrees).
pub const DEFAULT_FOV_Y_RAD: f32 = std::f32::consts::FRAC_PI_6;

/// Orbit sensitivity of [`Camera::orbit`], in radians per pixel of drag.
pub const ORBIT_RADIANS_PER_PIXEL: f64 = 0.008;

const MIN_DISTANCE: f64 = 1e-9;
const MAX_DISTANCE: f64 = 1e12;
/// Extra room around the bounding sphere when framing it.
const FIT_MARGIN: f64 = 1.15;
/// Lower bound of `near / far-extent` for perspective projections.
const MIN_NEAR_RATIO: f64 = 1e-5;

/// How the camera projects onto the screen.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Projection {
    /// Perspective projection with the given vertical field of view.
    Perspective {
        /// Vertical field of view in radians (clamped to 1..=170 degrees when used).
        fov_y_rad: f32,
    },
    /// Orthographic (parallel) projection.
    ///
    /// The visible height at the target is the same as a perspective camera
    /// with [`DEFAULT_FOV_Y_RAD`] at the same distance, so switching between
    /// the default perspective and orthographic keeps the model's size.
    /// Zooming changes [`Camera::distance`] in both modes.
    Orthographic,
}

impl Default for Projection {
    fn default() -> Self {
        Projection::Perspective {
            fov_y_rad: DEFAULT_FOV_Y_RAD,
        }
    }
}

/// Named standard view directions.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ViewPreset {
    /// Looking from front-right-top (eye towards +X, -Y, +Z): shows the
    /// front, right and top faces equally.
    Isometric,
    /// Looking along +Y (eye on the -Y side), +X right, +Z up.
    Front,
    /// Looking along -Y (eye on the +Y side), -X right, +Z up.
    Back,
    /// Looking down along -Z, +X right, +Y up.
    Top,
    /// Looking up along +Z, +X right, -Y up.
    Bottom,
    /// Looking along -X (eye on the +X side), +Y right, +Z up.
    Right,
    /// Looking along +X (eye on the -X side), -Y right, +Z up.
    Left,
}

impl ViewPreset {
    /// `(yaw, pitch)` in radians for this view.
    pub fn yaw_pitch(self) -> (f64, f64) {
        match self {
            // atan(1/sqrt(2)): the classic isometric elevation, 35.26 degrees.
            ViewPreset::Isometric => (FRAC_PI_4, (1.0f64 / 2.0f64.sqrt()).atan()),
            ViewPreset::Front => (0.0, 0.0),
            ViewPreset::Back => (std::f64::consts::PI, 0.0),
            ViewPreset::Top => (0.0, FRAC_PI_2),
            ViewPreset::Bottom => (0.0, -FRAC_PI_2),
            ViewPreset::Right => (FRAC_PI_2, 0.0),
            ViewPreset::Left => (-FRAC_PI_2, 0.0),
        }
    }
}

/// Orbit camera. See the [module documentation](self) for conventions.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Camera {
    /// Point the camera orbits around and looks at.
    pub target: DVec3,
    /// Distance from the eye to [`Camera::target`] (kept positive).
    pub distance: f64,
    /// Rotation around world +Z in radians; 0 looks along +Y, positive values
    /// move the eye towards +X.
    pub yaw: f64,
    /// Elevation in radians, clamped to `[-pi/2, pi/2]`; positive looks down.
    pub pitch: f64,
    /// Perspective or orthographic projection.
    pub projection: Projection,
    /// Axis-aligned bounds `(min, max)` of the scene, used to choose the clip
    /// planes. Set by [`Camera::fit`] / [`Camera::set_scene_bounds`].
    /// `None` uses generous planes derived from the distance alone.
    #[serde(default)]
    pub scene_bounds: Option<(DVec3, DVec3)>,
}

impl Default for Camera {
    fn default() -> Self {
        Camera::isometric()
    }
}

impl Camera {
    /// Camera at a preset view direction looking at the origin from distance 10,
    /// with the default perspective projection. Call [`Camera::fit`] to frame a model.
    pub fn from_preset(preset: ViewPreset) -> Camera {
        let (yaw, pitch) = preset.yaw_pitch();
        Camera {
            target: DVec3::ZERO,
            distance: 10.0,
            yaw,
            pitch,
            projection: Projection::default(),
            scene_bounds: None,
        }
    }

    /// Isometric view (front, right and top faces visible).
    pub fn isometric() -> Camera {
        Camera::from_preset(ViewPreset::Isometric)
    }

    /// Front view: looking along +Y, +Z up.
    pub fn front() -> Camera {
        Camera::from_preset(ViewPreset::Front)
    }

    /// Top view: looking down along -Z, +Y up on screen.
    pub fn top() -> Camera {
        Camera::from_preset(ViewPreset::Top)
    }

    /// Right view: looking along -X, +Z up.
    pub fn right() -> Camera {
        Camera::from_preset(ViewPreset::Right)
    }

    /// Returns `self` with the given projection (builder style).
    pub fn with_projection(mut self, projection: Projection) -> Camera {
        self.projection = projection;
        self
    }

    /// Returns `self` fitted to the bounds (builder style); see [`Camera::fit`].
    pub fn fitted(mut self, min: [f32; 3], max: [f32; 3]) -> Camera {
        self.fit(min, max);
        self
    }

    /// Switches to a preset view direction, keeping target, distance and projection.
    pub fn set_view(&mut self, preset: ViewPreset) {
        let (yaw, pitch) = preset.yaw_pitch();
        self.yaw = yaw;
        self.pitch = pitch;
    }

    /// Pitch actually used for the matrices (clamped, NaN-safe).
    fn clamped_pitch(&self) -> f64 {
        if self.pitch.is_nan() {
            0.0
        } else {
            self.pitch.clamp(-FRAC_PI_2, FRAC_PI_2)
        }
    }

    /// Distance actually used for the matrices (positive, finite).
    fn clamped_distance(&self) -> f64 {
        if self.distance.is_finite() {
            self.distance.clamp(MIN_DISTANCE, MAX_DISTANCE)
        } else {
            MAX_DISTANCE
        }
    }

    /// Vertical field of view used for perspective, or the reference field
    /// of view that sets the orthographic view height.
    pub fn fov_y(&self) -> f64 {
        let fov = match self.projection {
            Projection::Perspective { fov_y_rad } => fov_y_rad as f64,
            Projection::Orthographic => DEFAULT_FOV_Y_RAD as f64,
        };
        if fov.is_finite() {
            fov.clamp(1f64.to_radians(), 170f64.to_radians())
        } else {
            DEFAULT_FOV_Y_RAD as f64
        }
    }

    /// Unit vector from the target towards the eye.
    pub fn back_dir(&self) -> DVec3 {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.clamped_pitch().sin_cos();
        DVec3::new(sy * cp, -cy * cp, sp)
    }

    /// Unit viewing direction (from the eye towards the target).
    pub fn forward(&self) -> DVec3 {
        -self.back_dir()
    }

    /// Unit screen-right direction in world space (always horizontal).
    pub fn right_dir(&self) -> DVec3 {
        let (sy, cy) = self.yaw.sin_cos();
        DVec3::new(cy, sy, 0.0)
    }

    /// Unit screen-up direction in world space.
    pub fn up_dir(&self) -> DVec3 {
        self.right_dir().cross(self.forward())
    }

    /// Eye position in world space.
    pub fn eye(&self) -> DVec3 {
        self.target + self.back_dir() * self.clamped_distance()
    }

    /// Rotates around the target as if dragging the model: dragging right
    /// (positive `dx_pixels`) turns the model to the right, dragging down
    /// (positive `dy_pixels`, screen y grows downwards) tilts its top
    /// towards the viewer. Uses [`ORBIT_RADIANS_PER_PIXEL`].
    pub fn orbit(&mut self, dx_pixels: f32, dy_pixels: f32) {
        self.orbit_radians(
            -(dx_pixels as f64) * ORBIT_RADIANS_PER_PIXEL,
            dy_pixels as f64 * ORBIT_RADIANS_PER_PIXEL,
        );
    }

    /// Adds `d_yaw` and `d_pitch` (radians) to the camera angles, clamping pitch.
    pub fn orbit_radians(&mut self, d_yaw: f64, d_pitch: f64) {
        if d_yaw.is_finite() {
            self.yaw += d_yaw;
        }
        if d_pitch.is_finite() {
            self.pitch = (self.clamped_pitch() + d_pitch).clamp(-FRAC_PI_2, FRAC_PI_2);
        }
    }

    /// Magnifies the view by `factor`: the distance is divided by it, so
    /// `zoom(2.0)` makes the model appear twice as large and `zoom(0.5)`
    /// half as large. Non-positive or non-finite factors are ignored.
    pub fn zoom(&mut self, factor: f64) {
        if factor.is_finite() && factor > 0.0 {
            self.distance = (self.clamped_distance() / factor).clamp(MIN_DISTANCE, MAX_DISTANCE);
        }
    }

    /// World-space size of one pixel at the target's depth for a viewport
    /// `viewport_height` pixels tall.
    pub fn world_units_per_pixel(&self, viewport_height: f32) -> f64 {
        let h = 2.0 * self.clamped_distance() * (self.fov_y() * 0.5).tan();
        h / (viewport_height.max(1.0) as f64)
    }

    /// Moves the target in the view plane so the scene follows the cursor:
    /// a drag of `(dx, dy)` pixels (y down) moves points at the target depth
    /// by exactly that many pixels.
    pub fn pan(&mut self, dx_pixels: f32, dy_pixels: f32, viewport_height: f32) {
        let upp = self.world_units_per_pixel(viewport_height);
        let delta = self.right_dir() * (-(dx_pixels as f64) * upp)
            + self.up_dir() * (dy_pixels as f64 * upp);
        if delta.is_finite() {
            self.target += delta;
        }
    }

    /// Frames the axis-aligned box `min..max` for a square viewport (see
    /// [`Camera::fit_with_aspect`]); keeps the view direction and projection.
    pub fn fit(&mut self, min: [f32; 3], max: [f32; 3]) {
        self.fit_with_aspect(min, max, 1.0);
    }

    /// Frames the box `min..max`: the target moves to the box centre, the
    /// distance is chosen so the bounding sphere (plus a small margin) fits in
    /// both the vertical and horizontal field of view for a viewport with the
    /// given `aspect` (width / height), and the box is stored as the scene
    /// bounds for clip-plane selection. Degenerate boxes get a radius of 1.
    pub fn fit_with_aspect(&mut self, min: [f32; 3], max: [f32; 3], aspect: f32) {
        let lo = DVec3::from_array(min.map(f64::from));
        let hi = DVec3::from_array(max.map(f64::from));
        if !lo.is_finite() || !hi.is_finite() {
            return;
        }
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        let center = (lo + hi) * 0.5;
        let mut radius = (hi - lo).length() * 0.5;
        if radius <= 1e-12 * center.length().max(1.0) {
            radius = 1.0;
        }
        let aspect = if aspect.is_finite() && aspect > 0.0 {
            aspect as f64
        } else {
            1.0
        };
        let half_v = self.fov_y() * 0.5;
        let half_h = (half_v.tan() * aspect).atan();
        let half = half_v.min(half_h);
        let r = radius * FIT_MARGIN;
        self.distance = match self.projection {
            // Sphere tangent to the narrower frustum half-angle.
            Projection::Perspective { .. } => r / half.sin(),
            // Orthographic half-extent (distance * tan(half_v)) along the
            // narrower axis must cover the radius.
            Projection::Orthographic => r / (half_v.tan() * aspect.min(1.0)),
        }
        .clamp(MIN_DISTANCE, MAX_DISTANCE);
        self.target = center;
        self.scene_bounds = Some((lo, hi));
    }

    /// Records the scene bounds used to choose the clip planes without moving
    /// the camera (call when geometry changes).
    pub fn set_scene_bounds(&mut self, min: [f32; 3], max: [f32; 3]) {
        let lo = DVec3::from_array(min.map(f64::from));
        let hi = DVec3::from_array(max.map(f64::from));
        if lo.is_finite() && hi.is_finite() {
            self.scene_bounds = Some((lo.min(hi), lo.max(hi)));
        }
    }

    /// World-to-view matrix (right-handed, camera looks down -Z).
    pub fn view_matrix(&self) -> Mat4 {
        self.view_matrix_f64().as_mat4()
    }

    /// World-to-view matrix in double precision.
    pub fn view_matrix_f64(&self) -> DMat4 {
        let r = self.right_dir();
        let u = self.up_dir();
        let b = self.back_dir();
        let e = self.eye();
        DMat4::from_cols(
            DVec4::new(r.x, u.x, b.x, 0.0),
            DVec4::new(r.y, u.y, b.y, 0.0),
            DVec4::new(r.z, u.z, b.z, 0.0),
            DVec4::new(-r.dot(e), -u.dot(e), -b.dot(e), 1.0),
        )
    }

    /// Near and far clip distances (along the view direction) for the given
    /// scene bounds, or for [`Camera::scene_bounds`] when `bounds` is `None`.
    ///
    /// The planes bracket the scene's bounding sphere. For perspective the
    /// far plane is informational only (the projection uses an infinite far
    /// plane) and near never drops below `1e-5` of the far extent. For
    /// orthographic the near plane may be negative (behind the eye), so the
    /// scene is never clipped when zoomed in.
    pub fn clip_planes(&self, bounds: Option<(DVec3, DVec3)>) -> (f64, f64) {
        let d = self.clamped_distance();
        let eye = self.eye();
        let forward = self.forward();
        let (center_depth, radius) = match bounds.or(self.scene_bounds) {
            Some((lo, hi)) if lo.is_finite() && hi.is_finite() => {
                let c = (lo + hi) * 0.5;
                let r = ((hi - lo).length() * 0.5).max(1e-12 * c.length().max(1.0));
                ((c - eye).dot(forward), r * 1.02)
            }
            // Unknown scene: be generous around the target.
            _ => (d, d * 1e3),
        };
        let far = (center_depth + radius).max(d * 1e-3);
        match self.projection {
            Projection::Perspective { .. } => {
                let near = (center_depth - radius).max(far * MIN_NEAR_RATIO);
                (near, far)
            }
            Projection::Orthographic => {
                let near = center_depth - radius;
                (near, far.max(near + 1e-9 * d.max(1.0)))
            }
        }
    }

    /// Projection matrix for a viewport of the given `aspect` (width / height),
    /// with clip planes from [`Camera::scene_bounds`]. Reverse-Z: see the
    /// [module documentation](self).
    pub fn proj_matrix(&self, aspect: f32) -> Mat4 {
        self.proj_matrix_f64(aspect, None).as_mat4()
    }

    /// Like [`Camera::proj_matrix`] but with explicit scene bounds for the
    /// clip planes (the renderer passes the bounds of its uploaded meshes).
    pub fn proj_matrix_with_bounds(&self, aspect: f32, min: [f32; 3], max: [f32; 3]) -> Mat4 {
        let lo = DVec3::from_array(min.map(f64::from));
        let hi = DVec3::from_array(max.map(f64::from));
        self.proj_matrix_f64(aspect, Some((lo, hi))).as_mat4()
    }

    /// Double-precision projection matrix; `bounds` overrides the stored scene bounds.
    pub fn proj_matrix_f64(&self, aspect: f32, bounds: Option<(DVec3, DVec3)>) -> DMat4 {
        let aspect = if aspect.is_finite() && aspect > 0.0 {
            aspect as f64
        } else {
            1.0
        };
        let (near, far) = self.clip_planes(bounds);
        let half_v = self.fov_y() * 0.5;
        match self.projection {
            Projection::Perspective { .. } => {
                // Infinite reverse-Z: z_ndc = near / -z_view.
                let f = 1.0 / half_v.tan();
                DMat4::from_cols(
                    DVec4::new(f / aspect, 0.0, 0.0, 0.0),
                    DVec4::new(0.0, f, 0.0, 0.0),
                    DVec4::new(0.0, 0.0, 0.0, -1.0),
                    DVec4::new(0.0, 0.0, near, 0.0),
                )
            }
            Projection::Orthographic => {
                // Reverse-Z: z_view = -near -> 1, z_view = -far -> 0.
                let hh = self.clamped_distance() * half_v.tan();
                let hw = hh * aspect;
                let inv = 1.0 / (far - near);
                DMat4::from_cols(
                    DVec4::new(1.0 / hw, 0.0, 0.0, 0.0),
                    DVec4::new(0.0, 1.0 / hh, 0.0, 0.0),
                    DVec4::new(0.0, 0.0, inv, 0.0),
                    DVec4::new(0.0, 0.0, far * inv, 1.0),
                )
            }
        }
    }

    /// `proj_matrix(aspect) * view_matrix()`, computed in double precision.
    pub fn view_proj(&self, aspect: f32) -> Mat4 {
        (self.proj_matrix_f64(aspect, None) * self.view_matrix_f64()).as_mat4()
    }

    /// World-space ray `(origin, unit direction)` through pixel `(x, y)` of a
    /// `width` x `height` viewport (origin top-left, y down), for picking.
    pub fn pixel_ray(&self, x: f32, y: f32, width: f32, height: f32) -> (DVec3, DVec3) {
        let (w, h) = (width.max(1.0) as f64, height.max(1.0) as f64);
        let ndc_x = 2.0 * (x as f64) / w - 1.0;
        let ndc_y = 1.0 - 2.0 * (y as f64) / h;
        let half_h = (self.fov_y() * 0.5).tan();
        let half_w = half_h * w / h;
        let (r, u, f) = (self.right_dir(), self.up_dir(), self.forward());
        match self.projection {
            Projection::Perspective { .. } => {
                let dir = (f + r * (ndc_x * half_w) + u * (ndc_y * half_h)).normalize();
                (self.eye(), dir)
            }
            Projection::Orthographic => {
                let d = self.clamped_distance();
                let origin = self.eye() + r * (ndc_x * half_w * d) + u * (ndc_y * half_h * d);
                (origin, f)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn corners(min: [f32; 3], max: [f32; 3]) -> Vec<DVec3> {
        let mut v = Vec::new();
        for i in 0..8 {
            v.push(DVec3::new(
                if i & 1 == 0 { min[0] } else { max[0] } as f64,
                if i & 2 == 0 { min[1] } else { max[1] } as f64,
                if i & 4 == 0 { min[2] } else { max[2] } as f64,
            ));
        }
        v
    }

    #[test]
    fn presets_look_where_documented() {
        let close = |a: DVec3, b: DVec3| (a - b).length() < 1e-12;
        assert!(close(Camera::front().forward(), DVec3::Y));
        assert!(close(Camera::front().up_dir(), DVec3::Z));
        assert!(close(Camera::top().forward(), -DVec3::Z));
        assert!(close(Camera::top().up_dir(), DVec3::Y));
        assert!(close(Camera::right().forward(), -DVec3::X));
        assert!(close(Camera::right().right_dir(), DVec3::Y));
        let iso = Camera::isometric().back_dir();
        assert!(close(iso, DVec3::new(1.0, -1.0, 1.0).normalize()));
    }

    #[test]
    fn view_matrix_maps_eye_to_origin_and_target_forward() {
        let mut c = Camera::isometric();
        c.target = DVec3::new(3.0, -2.0, 7.0);
        c.distance = 5.0;
        let v = c.view_matrix_f64();
        assert!(v.transform_point3(c.eye()).length() < 1e-9);
        let t = v.transform_point3(c.target);
        assert!((t - DVec3::new(0.0, 0.0, -5.0)).length() < 1e-9);
    }

    #[test]
    fn pixel_ray_through_center_hits_target() {
        for proj in [Projection::default(), Projection::Orthographic] {
            let c = Camera::isometric().with_projection(proj);
            let (o, d) = c.pixel_ray(50.0, 40.0, 100.0, 80.0);
            let closest = o + d * (c.target - o).dot(d);
            assert!((closest - c.target).length() < 1e-9);
        }
    }

    #[test]
    fn serde_round_trip() {
        let mut c = Camera::top().with_projection(Projection::Orthographic);
        c.fit([0.0; 3], [1.0, 2.0, 3.0]);
        let json = serde_json::to_string(&c).unwrap();
        let back: Camera = serde_json::from_str(&json).unwrap();
        assert_eq!(c, back);
    }

    #[test]
    fn fit_contains_corners_for_many_scales() {
        for scale in [1e-3f32, 1.0, 1e4] {
            let min = [-scale, 2.0 * scale, -0.5 * scale];
            let max = [3.0 * scale, 2.5 * scale, 4.0 * scale];
            for preset in [ViewPreset::Isometric, ViewPreset::Top, ViewPreset::Left] {
                let c = Camera::from_preset(preset).fitted(min, max);
                let vp = c.proj_matrix_f64(1.0, None) * c.view_matrix_f64();
                for p in corners(min, max) {
                    let clip = vp * p.extend(1.0);
                    assert!(clip.w > 0.0);
                    let n = clip.truncate() / clip.w;
                    assert!(n.x.abs() <= 1.0 && n.y.abs() <= 1.0, "{n:?}");
                    assert!(n.z > 0.0 && n.z <= 1.0, "{n:?}");
                }
            }
        }
    }
}
