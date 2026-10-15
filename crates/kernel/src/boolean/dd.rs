//! Double-double numbers (Dekker 1971, "A floating-point technique for
//! extending the available precision"), for the boolean module's tests:
//! about 106 bits, for references that the rounding of `f64` does not
//! reach.

/// A number as the unevaluated sum of two `f64`, the second within half a
/// unit in the last place of the first.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Dd(pub f64, pub f64);

impl Dd {
    pub const ZERO: Dd = Dd(0.0, 0.0);
    pub const ONE: Dd = Dd(1.0, 0.0);

    pub fn of(a: f64) -> Dd {
        Dd(a, 0.0)
    }

    fn fast(a: f64, b: f64) -> Dd {
        let s = a + b;
        Dd(s, b - (s - a))
    }

    pub fn add(self, o: Dd) -> Dd {
        // Knuth's two-sum of the leading parts, then the rest.
        let s = self.0 + o.0;
        let v = s - self.0;
        let e = (self.0 - (s - v)) + (o.0 - v);
        Dd::fast(s, e + self.1 + o.1)
    }

    pub fn neg(self) -> Dd {
        Dd(-self.0, -self.1)
    }

    pub fn sub(self, o: Dd) -> Dd {
        self.add(o.neg())
    }

    pub fn mul(self, o: Dd) -> Dd {
        let p = self.0 * o.0;
        let e = self.0.mul_add(o.0, -p);
        Dd::fast(p, e + (self.0 * o.1 + self.1 * o.0))
    }

    pub fn mulf(self, k: f64) -> Dd {
        self.mul(Dd::of(k))
    }

    pub fn div(self, o: Dd) -> Dd {
        let q1 = self.0 / o.0;
        let r = self.sub(o.mul(Dd::of(q1)));
        let q2 = r.0 / o.0;
        let r = r.sub(o.mul(Dd::of(q2)));
        Dd::fast(q1, q2).add(Dd::of(r.0 / o.0))
    }

    pub fn sqrt(self) -> Dd {
        if self.0 <= 0.0 {
            return Dd::ZERO;
        }
        let x = Dd::of(self.0.sqrt());
        x.add(self.sub(x.mul(x)).div(x.add(x)))
    }

    /// The nearest `f64` (within a unit in its last place).
    pub fn f(self) -> f64 {
        self.0 + self.1
    }

    /// Whether this is less than `o`.
    pub fn lt(self, o: Dd) -> bool {
        self.sub(o).0 < 0.0
    }

    /// Whether this is at most `b`.
    pub fn le(self, b: f64) -> bool {
        self.sub(Dd::of(b)).0 <= 0.0
    }

    /// The larger of this and `o`.
    pub fn max(self, o: Dd) -> Dd {
        if self.lt(o) {
            o
        } else {
            self
        }
    }

    /// Sine and cosine by their Taylor series about zero, for `|x| ≤ 4`:
    /// the terms are summed until they fall under `10⁻⁴⁰`; each is at most
    /// `sinh 4 < 28`, so the sum is within about `10⁻³⁰` absolute.
    pub fn sin_cos(self) -> (Dd, Dd) {
        assert!(self.0.abs() <= 4.0, "{self:?}");
        let x2 = self.mul(self);
        let series = |first: Dd, mut n: f64| {
            let (mut sum, mut term) = (first, first);
            while term.0 != 0.0 {
                term = term.mul(x2).div(Dd::of((n + 1.0) * (n + 2.0))).neg();
                sum = sum.add(term);
                n += 2.0;
                if term.0.abs() < 1e-40 {
                    break;
                }
            }
            sum
        };
        (series(self, 1.0), series(Dd::ONE, 0.0))
    }

    /// `atan2(y, x)`: the `f64` angle of the leading parts, corrected by
    /// one Newton step, `θ = θ₀ + (y cos θ₀ − x sin θ₀) / (x cos θ₀ + y sin
    /// θ₀)` (the tangent of the remaining angle, whose cube, about 10⁻⁴⁸,
    /// is dropped).
    pub fn atan2(y: Dd, x: Dd) -> Dd {
        let a = y.0.atan2(x.0);
        let (s, c) = Dd::of(a).sin_cos();
        let num = y.mul(c).sub(x.mul(s));
        let den = x.mul(c).add(y.mul(s));
        Dd::of(a).add(num.div(den))
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::Dd;
    use crate::boolean::testing::Rng;

    /// The references hold to about 10⁻³⁰: sine and cosine, atan2 and the
    /// square root against identities, and sin 1 and cos 1 against their
    /// decimal values.
    #[test]
    fn references_hold_to_double_double() {
        let mut rng = Rng(0xD0D0_0001);
        for _ in 0..1000 {
            let a = Dd::of(rng.range(-PI, PI)).add(Dd::of(rng.range(-1e-17, 1e-17)));
            let (s, c) = a.sin_cos();
            let one = s.mul(s).add(c.mul(c)).sub(Dd::ONE).f();
            assert!(one.abs() < 1e-30, "{one}");
            let back = Dd::atan2(s.mulf(3.7), c.mulf(3.7)).sub(a).f();
            assert!(back.abs() < 1e-30, "{back}");
            let r = Dd::of(10f64.powf(rng.range(-12.0, 12.0)));
            let root = r.sqrt();
            assert!((root.mul(root).sub(r).f() / r.f()).abs() < 1e-30);
            let q = r.div(Dd::of(rng.range(0.5, 2.0)));
            assert!(q.0.is_finite());
        }
        // sin 1 and cos 1 as double-doubles (their series summed in exact
        // rationals, sixty terms).
        let (s, c) = Dd::of(1.0).sin_cos();
        let near = |x: Dd, hi: f64, lo: f64| x.sub(Dd(hi, 0.0)).sub(Dd::of(lo)).f().abs();
        assert!(near(s, 0.8414709848078965, 1.776845092935536e-18) < 1e-30);
        assert!(near(c, 0.5403023058681398, -4.760954612604417e-17) < 1e-30);
    }
}
