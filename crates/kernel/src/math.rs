//! Small math helpers shared by the kernel: frames and tolerances.

use glam::{DMat4, DVec2, DVec3};
use serde::{Deserialize, Serialize};

/// Model tolerances. All lengths are in metres, angles in radians.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Tolerance {
    /// Two points closer than this are the same point.
    pub linear: f64,
    /// Two directions closer than this angle are parallel.
    pub angular: f64,
}

impl Default for Tolerance {
    fn default() -> Self {
        // 0.1 µm, comparable to commercial kernels' 1e-7 m resolution.
        Self {
            linear: 1e-7,
            angular: 1e-9,
        }
    }
}

/// A right-handed orthonormal coordinate frame. Used for sketch planes,
/// primitive placement and surface parameterization.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    /// Origin.
    pub origin: DVec3,
    /// Unit x axis.
    pub x: DVec3,
    /// Unit y axis.
    pub y: DVec3,
    /// Unit z axis (`x × y`).
    pub z: DVec3,
}

impl Default for Frame {
    fn default() -> Self {
        Self::WORLD
    }
}

impl Frame {
    /// The world XY frame.
    pub const WORLD: Frame = Frame {
        origin: DVec3::ZERO,
        x: DVec3::X,
        y: DVec3::Y,
        z: DVec3::Z,
    };

    /// Builds a frame from an origin, a normal and an x direction hint.
    /// The hint is projected to be perpendicular to the normal; if it is
    /// parallel to the normal an arbitrary perpendicular is chosen.
    pub fn from_normal(origin: DVec3, normal: DVec3, x_hint: DVec3) -> Frame {
        let z = normal.normalize();
        let mut x = x_hint - z * x_hint.dot(z);
        if x.length_squared() < 1e-24 {
            x = z.any_orthonormal_vector();
        }
        let x = x.normalize();
        let y = z.cross(x);
        Frame { origin, x, y, z }
    }

    /// The world XY plane (top), with normal +Z.
    pub fn top() -> Frame {
        Self::WORLD
    }

    /// The world XZ plane (front), with normal −Y.
    pub fn front() -> Frame {
        Frame {
            origin: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Z,
            z: -DVec3::Y,
        }
    }

    /// The world YZ plane (right), with normal +X.
    pub fn right() -> Frame {
        Frame {
            origin: DVec3::ZERO,
            x: DVec3::Y,
            y: DVec3::Z,
            z: DVec3::X,
        }
    }

    /// Maps local 2D plane coordinates to world space.
    pub fn point2(&self, p: DVec2) -> DVec3 {
        self.origin + self.x * p.x + self.y * p.y
    }

    /// Maps a local 2D direction to world space.
    pub fn vector2(&self, v: DVec2) -> DVec3 {
        self.x * v.x + self.y * v.y
    }

    /// Maps local 3D coordinates to world space.
    pub fn point(&self, p: DVec3) -> DVec3 {
        self.origin + self.x * p.x + self.y * p.y + self.z * p.z
    }

    /// Expresses a world point in local coordinates.
    pub fn to_local(&self, p: DVec3) -> DVec3 {
        let d = p - self.origin;
        DVec3::new(d.dot(self.x), d.dot(self.y), d.dot(self.z))
    }

    /// The local-to-world affine matrix.
    pub fn matrix(&self) -> DMat4 {
        DMat4::from_cols(
            self.x.extend(0.0),
            self.y.extend(0.0),
            self.z.extend(0.0),
            self.origin.extend(1.0),
        )
    }

    /// The frame translated by `d`.
    pub fn translated(&self, d: DVec3) -> Frame {
        Frame {
            origin: self.origin + d,
            ..*self
        }
    }

    /// True if the axes are orthonormal and right-handed within `eps`.
    pub fn is_valid(&self, eps: f64) -> bool {
        (self.x.length() - 1.0).abs() < eps
            && (self.y.length() - 1.0).abs() < eps
            && self.x.dot(self.y).abs() < eps
            && (self.x.cross(self.y) - self.z).length() < eps
            && self.origin.is_finite()
    }
}

/// Rotates `v` about the unit `axis` by `angle` (Rodrigues' formula).
pub(crate) fn rotate_vec(v: DVec3, axis: DVec3, angle: f64) -> DVec3 {
    let (s, c) = angle.sin_cos();
    v * c + axis.cross(v) * s + axis * axis.dot(v) * (1.0 - c)
}

/// Rotates the point `p` about the axis through `origin` along unit `axis`.
#[cfg(test)]
pub(crate) fn rotate_point(p: DVec3, origin: DVec3, axis: DVec3, angle: f64) -> DVec3 {
    origin + rotate_vec(p - origin, axis, angle)
}

/// Affine rotation matrix about an axis through a point.
pub(crate) fn rotation_about(origin: DVec3, axis: DVec3, angle: f64) -> DMat4 {
    DMat4::from_translation(origin)
        * DMat4::from_axis_angle(axis, angle)
        * DMat4::from_translation(-origin)
}

/// Wraps an angle into `[-π/2, 3π/2)`. Periodic parameters in the kernel use
/// this range; faces never span more than π, and every face's frame starts
/// at its own boundary, so a face's domain `[0, ≤π]` never straddles the cut.
pub(crate) fn wrap_angle(a: f64) -> f64 {
    use std::f64::consts::{FRAC_PI_2, TAU};
    let mut x = (a + FRAC_PI_2).rem_euclid(TAU) - FRAC_PI_2;
    if x >= 3.0 * FRAC_PI_2 {
        x -= TAU;
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_angle_range() {
        use std::f64::consts::PI;
        for i in -100..100 {
            let a = i as f64 * 0.173;
            let w = wrap_angle(a);
            assert!((-PI / 2.0..1.5 * PI).contains(&w));
            assert!(((w - a) / (2.0 * PI)).round() * 2.0 * PI - (w - a) < 1e-12);
        }
        assert!((wrap_angle(PI) - PI).abs() < 1e-15);
        assert!(wrap_angle(0.0).abs() < 1e-15);
    }

    #[test]
    fn frames_are_valid() {
        for f in [Frame::top(), Frame::front(), Frame::right()] {
            assert!(f.is_valid(1e-12));
        }
        let f = Frame::from_normal(DVec3::ONE, DVec3::new(1.0, 2.0, 3.0), DVec3::X);
        assert!(f.is_valid(1e-12));
        let p = DVec3::new(0.3, -2.0, 5.0);
        assert!((f.point(f.to_local(p)) - p).length() < 1e-12);
    }

    #[test]
    fn rotation_matches_matrix() {
        let o = DVec3::new(1.0, 2.0, 3.0);
        let axis = DVec3::new(0.0, 1.0, 1.0).normalize();
        let p = DVec3::new(-1.0, 0.5, 2.0);
        let a = rotate_point(p, o, axis, 0.7);
        let b = rotation_about(o, axis, 0.7).transform_point3(p);
        assert!((a - b).length() < 1e-12);
    }
}
