//! Helpers for the boolean module's unit tests.

use glam::{DVec2, DVec3};

use crate::math::Frame;
use crate::profile::{Profile, ProfileCurve};
use crate::sweep::{self, RevolveSpec};
use crate::topo::Brep;

/// A small deterministic generator (xorshift), for repeatable tests.
pub(crate) struct Rng(pub u64);

impl Rng {
    /// Uniform in `[0, 1)`.
    pub fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[a, b)`.
    pub fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }

    /// Uniform in `0..n`.
    pub fn below(&mut self, n: usize) -> usize {
        ((self.next() * n as f64) as usize).min(n - 1)
    }

    /// Uniform in the cube `[−r, r]³`.
    pub fn point(&mut self, r: f64) -> DVec3 {
        DVec3::new(self.range(-r, r), self.range(-r, r), self.range(-r, r))
    }

    /// A random unit vector.
    pub fn unit(&mut self) -> DVec3 {
        loop {
            let v = self.point(1.0);
            let l = v.length();
            if l > 0.1 && l <= 1.0 {
                return v / l;
            }
        }
    }
}

/// A solid turned a whole turn about the `y` axis of `plane` from the
/// profile `curves` in it, its faces of half a turn each (ADR 0008).
pub(super) fn turned(plane: Frame, curves: Vec<ProfileCurve>) -> Brep {
    let spec = RevolveSpec {
        axis_origin: DVec2::ZERO,
        axis_dir: DVec2::Y,
        angle: std::f64::consts::TAU,
    };
    sweep::revolve(
        &Profile { plane, curves },
        &spec,
        &crate::Tolerance::default(),
    )
    .expect("a turned solid")
    .0
}

/// A ball of radius `r` about the origin of `plane`: a half disc turned
/// about the plane's `y` axis, in two faces on the sphere.
pub(super) fn ball(plane: Frame, r: f64) -> Brep {
    turned(
        plane,
        vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: r,
                start_angle: -std::f64::consts::FRAC_PI_2,
                sweep_angle: std::f64::consts::PI,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, r),
                end: DVec2::new(0.0, -r),
            },
        ],
    )
}

/// A pointed cone on the plane's origin: its base a disc of radius `r`
/// square to the plane's `y` axis, its apex `h` along that axis.
pub(super) fn cone(plane: Frame, r: f64, h: f64) -> Brep {
    let corners = [DVec2::ZERO, DVec2::new(r, 0.0), DVec2::new(0.0, h)];
    turned(
        plane,
        (0..3)
            .map(|k| ProfileCurve::Line {
                start: corners[k],
                end: corners[(k + 1) % 3],
            })
            .collect(),
    )
}
