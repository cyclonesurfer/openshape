//! Helpers for the boolean module's unit tests.

use glam::DVec3;

/// A small deterministic generator (xorshift), for repeatable tests.
pub(super) struct Rng(pub u64);

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
