//! Polynomials on `[0, 1]` in Bernstein form and their real roots, for the
//! crossings of curved edges with surfaces (design note 0012 §10.2).
//!
//! A polynomial of degree `n` is `Σ cᵢ Bᵢⁿ(s)` with `Bᵢⁿ(s) = C(n, i) sⁱ
//! (1 − s)ⁿ⁻ⁱ`. The form is the numerically stable one on the interval
//! (Farouki & Rajan, "On the numerical condition of polynomials in Bernstein
//! form", *CAGD* 4, 1987), and its coefficients bound the polynomial: the
//! number of real roots in `(0, 1)` is at most the number of sign changes
//! of the coefficients, and of the same parity (the variation-diminishing
//! property; Lane & Riesenfeld, "Bounds on a polynomial", *BIT* 21, 1981).
//! Roots are isolated by halving the interval until each part has at most
//! one sign change, then refined by bisection on the polynomial's value.

/// Most coefficients a polynomial here has (degree 15): enough for the
/// distance to a quadric along a rational Bézier piece of degree five.
pub(super) const MAX: usize = 16;

/// Depth of halving after which an interval that still has more than one
/// sign change is taken to hold a cluster of roots (a near-multiple root),
/// reported at its middle: an interval of `2⁻⁴⁸` of the piece.
const DEPTH: u32 = 48;

/// A polynomial on `[0, 1]` in Bernstein form.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Bernstein {
    c: [f64; MAX],
    /// Number of coefficients: degree plus one (at least one).
    n: usize,
}

/// `C(n, k)` for `n < 2 · MAX`, exact in `f64`.
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut b = 1.0;
    for i in 0..k {
        b = b * (n - i) as f64 / (i + 1) as f64;
    }
    b.round()
}

impl Bernstein {
    /// The polynomial with these coefficients, or `None` with none or more
    /// than [`MAX`].
    pub fn new(coefficients: &[f64]) -> Option<Bernstein> {
        if coefficients.is_empty() || coefficients.len() > MAX {
            return None;
        }
        let mut c = [0.0; MAX];
        c[..coefficients.len()].copy_from_slice(coefficients);
        Some(Bernstein {
            c,
            n: coefficients.len(),
        })
    }

    /// The constant `v`.
    pub fn constant(v: f64) -> Bernstein {
        Bernstein::new(&[v]).expect("one coefficient")
    }

    /// The coefficients.
    pub fn coefficients(&self) -> &[f64] {
        &self.c[..self.n]
    }

    /// The degree.
    pub fn degree(&self) -> usize {
        self.n - 1
    }

    /// The value at `s`, by de Casteljau's algorithm.
    pub fn eval(&self, s: f64) -> f64 {
        let mut b = self.c;
        let t = 1.0 - s;
        for r in 1..self.n {
            for i in 0..self.n - r {
                b[i] = b[i] * t + b[i + 1] * s;
            }
        }
        b[0]
    }

    /// The two halves on `[0, ½]` and `[½, 1]`, each reparameterised to
    /// `[0, 1]`.
    pub fn halves(&self) -> (Bernstein, Bernstein) {
        let n = self.n;
        let mut left = [0.0; MAX];
        let mut right = [0.0; MAX];
        let mut b = self.c;
        left[0] = b[0];
        right[n - 1] = b[n - 1];
        for r in 1..n {
            for i in 0..n - r {
                b[i] = 0.5 * (b[i] + b[i + 1]);
            }
            left[r] = b[0];
            right[n - 1 - r] = b[n - 1 - r];
        }
        (Bernstein { c: left, n }, Bernstein { c: right, n })
    }

    /// The derivative (degree one less; the zero constant for a constant).
    pub fn deriv(&self) -> Bernstein {
        if self.n == 1 {
            return Bernstein::constant(0.0);
        }
        let d = self.degree() as f64;
        let mut c = [0.0; MAX];
        for (i, ci) in c.iter_mut().enumerate().take(self.n - 1) {
            *ci = d * (self.c[i + 1] - self.c[i]);
        }
        Bernstein { c, n: self.n - 1 }
    }

    /// The same polynomial written with degree `degree` (at least its own
    /// and less than [`MAX`]).
    pub fn elevated(&self, degree: usize) -> Bernstein {
        debug_assert!(degree >= self.degree() && degree < MAX);
        let mut p = *self;
        while p.degree() < degree {
            // One step (Farin §5.1): c'ᵢ = i/(n+1) cᵢ₋₁ + (1 − i/(n+1)) cᵢ.
            let m = p.n as f64;
            let mut c = [0.0; MAX];
            c[0] = p.c[0];
            c[p.n] = p.c[p.n - 1];
            for (i, ci) in c.iter_mut().enumerate().take(p.n).skip(1) {
                let a = i as f64 / m;
                *ci = a * p.c[i - 1] + (1.0 - a) * p.c[i];
            }
            p = Bernstein { c, n: p.n + 1 };
        }
        p
    }

    /// The product, or `None` if its degree would reach [`MAX`].
    pub fn mul(&self, other: &Bernstein) -> Option<Bernstein> {
        let (m, n) = (self.degree(), other.degree());
        if m + n >= MAX {
            return None;
        }
        let mut c = [0.0; MAX];
        for i in 0..=m {
            let bi = binomial(m, i) * self.c[i];
            for j in 0..=n {
                c[i + j] += bi * binomial(n, j) * other.c[j];
            }
        }
        for (k, ck) in c.iter_mut().enumerate().take(m + n + 1) {
            *ck /= binomial(m + n, k);
        }
        Some(Bernstein { c, n: m + n + 1 })
    }

    /// `self + k · other`, at the higher of the two degrees.
    pub fn add_scaled(&self, k: f64, other: &Bernstein) -> Bernstein {
        let d = self.degree().max(other.degree());
        let a = self.elevated(d);
        let b = other.elevated(d);
        let mut c = [0.0; MAX];
        for (i, ci) in c.iter_mut().enumerate().take(d + 1) {
            *ci = a.c[i] + k * b.c[i];
        }
        Bernstein { c, n: d + 1 }
    }

    /// `k` times the polynomial.
    pub fn scaled(&self, k: f64) -> Bernstein {
        let mut p = *self;
        for c in &mut p.c[..p.n] {
            *c *= k;
        }
        p
    }

    /// Sets to zero every coefficient of magnitude at most `floor`: those
    /// that are rounding noise of a computation whose terms were about
    /// `floor / 1e-12` in size, so that a polynomial that is zero but for
    /// rounding (the distance to a surface along an edge lying on it) has
    /// no spurious roots.
    pub fn cleaned(&self, floor: f64) -> Bernstein {
        let mut p = *self;
        for c in &mut p.c[..p.n] {
            if c.abs() <= floor {
                *c = 0.0;
            }
        }
        p
    }

    /// The largest coefficient magnitude: a bound on `|p(s)|` over `[0, 1]`
    /// (the Bernstein polynomials are non-negative and sum to one).
    pub fn bound(&self) -> f64 {
        self.c[..self.n].iter().fold(0.0, |m, c| m.max(c.abs()))
    }

    /// The real roots in the open interval `(0, 1)`, increasing. A simple
    /// root is found to within a few units of rounding of `s`; a cluster of
    /// roots closer together than `2⁻⁴⁸` (a near-multiple root) is
    /// reported once, at the middle of the interval holding it. The zero
    /// polynomial has none.
    pub fn roots(&self) -> Vec<f64> {
        let mut out = Vec::new();
        self.isolate(self, 0.0, 1.0, 0, &mut out);
        out.sort_by(f64::total_cmp);
        out.dedup_by(|a, b| (*a - *b).abs() <= 1e-15);
        out
    }

    /// Sign changes of the nonzero coefficients, with the sign of the first
    /// nonzero one (the polynomial's sign just right of 0).
    fn variations(&self) -> (usize, f64) {
        let mut count = 0;
        let mut last = 0.0;
        let mut first = 0.0;
        for &c in &self.c[..self.n] {
            if c == 0.0 {
                continue;
            }
            let s = c.signum();
            if last == 0.0 {
                first = s;
            } else if s != last {
                count += 1;
            }
            last = s;
        }
        (count, first)
    }

    /// Adds the roots of `whole` in `(lo, hi)` to `out`, `self` being
    /// `whole` on `[lo, hi]` reparameterised to `[0, 1]`.
    fn isolate(&self, whole: &Bernstein, lo: f64, hi: f64, depth: u32, out: &mut Vec<f64>) {
        let (v, first) = self.variations();
        if v == 0 {
            return;
        }
        if v == 1 {
            out.push(whole.bisect(lo, hi, first));
            return;
        }
        let mid = 0.5 * (lo + hi);
        if depth == DEPTH || mid <= lo || mid >= hi {
            out.push(mid);
            return;
        }
        let (left, right) = self.halves();
        left.isolate(whole, lo, mid, depth + 1, out);
        if left.c[left.n - 1] == 0.0 {
            out.push(mid);
        }
        right.isolate(whole, mid, hi, depth + 1, out);
    }

    /// The root of `self` in `(lo, hi)`, where it has exactly one and its
    /// sign just right of `lo` is `first`, by bisection to the last bit.
    fn bisect(&self, mut lo: f64, mut hi: f64, first: f64) -> f64 {
        loop {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                return mid;
            }
            let v = self.eval(mid);
            if v == 0.0 {
                return mid;
            }
            if v.signum() == first {
                lo = mid;
            } else {
                hi = mid;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::testing::Rng;
    use super::*;

    /// `∏ (s − rᵢ)` times `k`, in Bernstein form.
    fn from_roots(roots: &[f64], k: f64) -> Bernstein {
        let mut p = Bernstein::constant(k);
        for &r in roots {
            // s − r = −r · B₀¹ + (1 − r) · B₁¹.
            p = p.mul(&Bernstein::new(&[-r, 1.0 - r]).unwrap()).unwrap();
        }
        p
    }

    /// The same polynomial evaluated from its power form, for comparison.
    fn power_eval(roots: &[f64], k: f64, s: f64) -> f64 {
        roots.iter().fold(k, |v, r| v * (s - r))
    }

    #[test]
    fn arithmetic_matches_the_power_form() {
        let mut rng = Rng(0x5EED_0001);
        for _ in 0..200 {
            let n = 1 + (rng.next() * 7.0) as usize;
            let roots: Vec<f64> = (0..n).map(|_| rng.next() * 3.0 - 1.0).collect();
            let p = from_roots(&roots, 2.5);
            assert_eq!(p.degree(), n);
            for i in 0..=10 {
                let s = i as f64 / 10.0;
                let want = power_eval(&roots, 2.5, s);
                assert!((p.eval(s) - want).abs() <= 1e-12 * (1.0 + want.abs()));
                // The derivative, against a central difference.
                let h = 1e-6;
                let fd =
                    (power_eval(&roots, 2.5, s + h) - power_eval(&roots, 2.5, s - h)) / (2.0 * h);
                assert!((p.deriv().eval(s) - fd).abs() <= 1e-6 * (1.0 + fd.abs()));
                // Elevation and halving keep the polynomial.
                assert!((p.elevated(12).eval(s) - want).abs() <= 1e-12 * (1.0 + want.abs()));
                let (l, r) = p.halves();
                assert!((l.eval(s) - p.eval(0.5 * s)).abs() <= 1e-12 * (1.0 + want.abs()));
                assert!((r.eval(s) - p.eval(0.5 + 0.5 * s)).abs() <= 1e-12 * (1.0 + want.abs()));
            }
            let q = from_roots(&roots[..1], -1.0);
            let sum = p.add_scaled(3.0, &q);
            for i in 0..=10 {
                let s = i as f64 / 10.0;
                let want = power_eval(&roots, 2.5, s) + 3.0 * power_eval(&roots[..1], -1.0, s);
                assert!((sum.eval(s) - want).abs() <= 1e-12 * (1.0 + want.abs()));
            }
        }
        assert!(Bernstein::new(&[0.0; MAX + 1]).is_none());
        let big = Bernstein::new(&[1.0; 9]).unwrap();
        assert!(big.mul(&big).is_none());
    }

    #[test]
    fn simple_roots_are_found_to_rounding() {
        let mut rng = Rng(0x5EED_0002);
        for case in 0..2000 {
            let n = 1 + case % 10;
            // Roots inside and outside the interval, at least 1e-3 apart.
            let mut roots: Vec<f64> = Vec::new();
            while roots.len() < n {
                let r = rng.next() * 1.6 - 0.3;
                let apart = roots.iter().all(|q| (q - r).abs() > 1e-3);
                if apart && r.abs() > 1e-6 && (r - 1.0).abs() > 1e-6 {
                    roots.push(r);
                }
            }
            let k = if case % 2 == 0 { 1.0 } else { -7.0 };
            let p = from_roots(&roots, k);
            let mut want: Vec<f64> = roots
                .iter()
                .copied()
                .filter(|r| *r > 0.0 && *r < 1.0)
                .collect();
            want.sort_by(f64::total_cmp);
            let got = p.roots();
            assert_eq!(got.len(), want.len(), "case {case}: {roots:?} gave {got:?}");
            for (g, w) in got.iter().zip(&want) {
                // Within rounding of the coefficients over the slope there
                // (the root's condition number).
                let slope: f64 = k * roots
                    .iter()
                    .filter(|r| *r != w)
                    .map(|r| w - r)
                    .product::<f64>();
                let bound = 1e-14 * p.bound() / slope.abs() + 1e-15;
                assert!(
                    (g - w).abs() <= bound,
                    "case {case}: {g} vs {w}, bound {bound}"
                );
            }
        }
    }

    #[test]
    fn multiple_and_clustered_roots_are_reported_once() {
        // A double root: the polynomial touches zero without changing sign.
        let p = from_roots(&[0.3, 0.3], 1.0);
        let r = p.roots();
        assert!(r.iter().all(|x| (x - 0.3).abs() < 1e-7), "{r:?}");
        assert!(r.len() <= 2);
        // Two roots 1e-16 apart are one cluster.
        let p = from_roots(&[0.6, 0.6 + 1e-16, 0.2], 1.0);
        let r = p.roots();
        assert!((r[0] - 0.2).abs() < 1e-12);
        assert!(r[1..].iter().all(|x| (x - 0.6).abs() < 1e-7), "{r:?}");
        // A root exactly at a halving point.
        let p = from_roots(&[0.5, 0.25, 0.75], 1.0);
        let r = p.roots();
        assert_eq!(r.len(), 3, "{r:?}");
        for (g, w) in r.iter().zip([0.25, 0.5, 0.75]) {
            assert!((g - w).abs() < 1e-14);
        }
        // Roots at the ends are not reported.
        assert!(from_roots(&[0.0, 1.0], 1.0).roots().is_empty());
        // The zero polynomial, and a polynomial that is zero up to rounding
        // once cleaned, have none.
        assert!(Bernstein::new(&[0.0; 5]).unwrap().roots().is_empty());
        let noise = Bernstein::new(&[1e-17, -2e-17, 3e-17, -1e-17]).unwrap();
        assert!(noise.cleaned(1e-12 * 1e-4).roots().is_empty());
    }

    #[test]
    fn noise_takes_bounded_work() {
        // Random noise of degree 15: every root is reported at most once per
        // sign change, so at most fifteen.
        let mut rng = Rng(0x5EED_0003);
        for _ in 0..100 {
            let c: Vec<f64> = (0..MAX).map(|_| rng.next() - 0.5).collect();
            let p = Bernstein::new(&c).unwrap();
            let r = p.roots();
            assert!(r.len() < MAX);
            for &s in &r {
                // Each is a root up to rounding, or the middle of a cluster.
                assert!(p.eval(s).abs() < 1e-9, "{s}: {}", p.eval(s));
            }
        }
    }
}
