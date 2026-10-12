//! Mass properties by boundary integration on the exact geometry.
//!
//! Volume integrals are converted to surface integrals with the divergence
//! theorem, and each face's surface integral is converted to an integral
//! around the face's boundary in parameter space with Green's theorem
//! (`docs/design/0003-mass-properties.md`):
//!
//! ```text
//! ∫_V g dV = ∮_∂V F·n dA,          div F = g
//! ∬_D f du dv = ∮_∂D (∫_{u0}^{u} f(s, v) ds) dv        (inner along u)
//!             = −∮_∂D (∫_{v0}^{v} f(u, s) ds) du       (inner along v)
//! ```
//!
//! Both the inner integrals and the boundary integrals use Gauss–Legendre
//! quadrature on the exact surfaces and edge curves, so the result does not
//! depend on a tessellation. The boundary integrals are adaptive (see
//! [`Adaptive`]) to the requested accuracy, which for smooth geometry is
//! reached with few evaluations and can be asked down to rounding.

use glam::{DMat3, DVec2, DVec3};

use crate::geom::Surface;
use crate::topo::{Brep, Face};
use crate::{KernelError, MassProperties};

/// Gauss–Legendre nodes and weights on [-1, 1] (Newton iteration on P_n).
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

/// The eleven accumulated integrals: area, volume, first moments (x, y, z)
/// and second moments (xx, yy, zz, xy, yz, zx), all about `origin`.
type Acc = [f64; 11];

fn add(a: &mut Acc, b: &Acc, s: f64) {
    for i in 0..11 {
        a[i] += b[i] * s;
    }
}

/// Integrand vector at parameter point `uv`, with coordinates relative to `o`.
fn integrand(s: &Surface, uv: DVec2, o: DVec3) -> Acc {
    let (su, sv) = s.derivs(uv.x, uv.y);
    integrand_at(s.point(uv.x, uv.y) - o, su, sv)
}

/// Integrand vector at a surface point `p` (relative to the origin) with
/// partial derivatives `su`, `sv`.
fn integrand_at(p: DVec3, su: DVec3, sv: DVec3) -> Acc {
    let n = su.cross(sv);
    [
        n.length(),
        p.dot(n) / 3.0,
        p.x * p.x / 2.0 * n.x,
        p.y * p.y / 2.0 * n.y,
        p.z * p.z / 2.0 * n.z,
        p.x * p.x * p.x / 3.0 * n.x,
        p.y * p.y * p.y / 3.0 * n.y,
        p.z * p.z * p.z / 3.0 * n.z,
        p.x * p.x * p.y / 2.0 * n.x,
        p.y * p.y * p.z / 2.0 * n.y,
        p.z * p.z * p.x / 2.0 * n.z,
    ]
}

struct Ctx<'a> {
    surface: &'a Surface,
    origin: DVec3,
    /// True: inner integral along u (accumulate `F dv`); false: along v
    /// (accumulate `−G du`).
    inner_u: bool,
    reference: f64,
    inner: Vec<(f64, f64)>,
}

impl Ctx<'_> {
    /// Inner integral from the reference line to `uv`.
    fn inner_integral(&self, uv: DVec2) -> Acc {
        let mut acc = [0.0; 11];
        let (a, b) = if self.inner_u {
            (self.reference, uv.x)
        } else {
            (self.reference, uv.y)
        };
        let half = 0.5 * (b - a);
        if half == 0.0 {
            return acc;
        }
        // Along an extrusion's v the base curve's point and tangent are
        // fixed: evaluate them once.
        if let (Surface::Extrusion { curve, dir, .. }, false) = (self.surface, self.inner_u) {
            let base = curve.point(uv.x) - self.origin;
            let su = curve.deriv(uv.x);
            for &(x, w) in &self.inner {
                let t = a + half * (x + 1.0);
                add(&mut acc, &integrand_at(base + *dir * t, su, *dir), w * half);
            }
            return acc;
        }
        for &(x, w) in &self.inner {
            let t = a + half * (x + 1.0);
            let q = if self.inner_u {
                DVec2::new(t, uv.y)
            } else {
                DVec2::new(uv.x, t)
            };
            add(&mut acc, &integrand(self.surface, q, self.origin), w * half);
        }
        acc
    }

    /// Contribution of an infinitesimal boundary step `duv` at `uv`.
    fn boundary(&self, uv: DVec2, duv: DVec2) -> (Acc, f64) {
        if self.inner_u {
            (self.inner_integral(uv), duv.y)
        } else {
            (self.inner_integral(uv), -duv.x)
        }
    }
}

/// Parameter-space velocity of a curve point moving with 3D velocity `c1`
/// on the surface: least-squares solution of `[S_u S_v] · duv = c1`.
fn param_velocity(s: &Surface, uv: DVec2, c1: DVec3) -> DVec2 {
    let (su, sv) = s.derivs(uv.x, uv.y);
    let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let det = a * c - b * b;
    let (r0, r1) = (su.dot(c1), sv.dot(c1));
    if det.abs() > 1e-14 * a.max(c).max(1e-300).powi(2) {
        DVec2::new((c * r0 - b * r1) / det, (a * r1 - b * r0) / det)
    } else if a < c {
        DVec2::new(0.0, r1 / c)
    } else {
        DVec2::new(r0 / a, 0.0)
    }
}

/// Parameter value of a coedge end, resolving the ambiguous coordinate at a
/// singular point from a point slightly inside the coedge.
fn end_uv(s: &Surface, p_end: DVec3, p_inside: DVec3) -> DVec2 {
    let (u, v) = s.project(p_end);
    let uv = DVec2::new(u, v);
    let (su, sv) = s.derivs(u, v);
    let (a, b) = (su.length(), sv.length());
    let (iu, iv) = s.project(p_inside);
    if a < 1e-9 * b.max(1e-300) {
        DVec2::new(iu, v)
    } else if b < 1e-9 * a.max(1e-300) {
        DVec2::new(u, iv)
    } else {
        uv
    }
}

/// The deepest bisection of a boundary piece. A piece that still misses
/// its tolerance there makes the whole computation fail (at once) rather
/// than return a value it cannot vouch for.
const MAX_DEPTH: u32 = 30;

/// Differences below this fraction of a piece's absolute integral are
/// rounding noise: the integrand is evaluated in floating point, so no
/// rule resolves finer differences (the termination criterion of Gander &
/// Gautschi, below). Requested accuracies finer than this mean this.
const NOISE: f64 = 1e-13;

/// Adaptive Gauss–Legendre quadrature over one smooth piece of an edge
/// (Gander & Gautschi, "Adaptive Quadrature — Revisited", BIT 40, 2000).
struct Adaptive<'a> {
    /// The integrand.
    g: &'a dyn Fn(f64) -> Acc,
    /// The rule (8 points: exact to degree 15).
    gl: &'a [(f64, f64)],
    /// Allowed error as a fraction of a piece's absolute integral.
    relative: f64,
    /// Allowed error per component per unit of parameter length regardless
    /// of the integral (the rounding floor of a component that vanishes on
    /// the piece up to rounding).
    floor: Acc,
}

/// The result of [`Adaptive::run`] on a piece: its integral, and whether
/// every piece met its tolerance.
struct Quad {
    value: Acc,
    converged: bool,
}

impl Adaptive<'_> {
    /// Integrates over `[a, b]`, given the rule's value `whole` there.
    ///
    /// The rule on the two halves is compared with `whole`; the halves are
    /// accepted when every component differs by at most the relative
    /// tolerance times the integral of the component's magnitude over the
    /// piece, or by at most its floor times the piece's length (the
    /// halves' own error is far smaller for smooth integrands). Otherwise
    /// each half is refined in turn.
    ///
    /// The estimate assumes the integrand is smooth on the piece, which
    /// the splitting into pieces at knots ensures: a jump between samples
    /// is invisible to any sampling rule.
    fn run(&self, a: f64, b: f64, whole: &Acc, depth: u32, out: &mut Quad) {
        if !out.converged {
            // The result is already known to be refused.
            return;
        }
        let m = 0.5 * (a + b);
        let (left, left_abs) = rule(self.g, a, m, self.gl);
        let (right, right_abs) = rule(self.g, m, b, self.gl);
        let len = (b - a).abs();
        let ok = (0..11).all(|k| {
            let d = (left[k] + right[k] - whole[k]).abs();
            d <= self.relative * (left_abs[k] + right_abs[k]) || d <= self.floor[k] * len
        });
        if ok {
            add(&mut out.value, &left, 1.0);
            add(&mut out.value, &right, 1.0);
        } else if depth >= MAX_DEPTH || m == a || m == b {
            add(&mut out.value, &left, 1.0);
            add(&mut out.value, &right, 1.0);
            out.converged = false;
        } else {
            self.run(a, m, &left, depth + 1, out);
            self.run(m, b, &right, depth + 1, out);
        }
    }
}

/// The Gauss–Legendre rule `gl` for `g` on `[a, b]`, and the same rule
/// for `|g|` (the scale of the rounding noise in the first).
fn rule(g: &dyn Fn(f64) -> Acc, a: f64, b: f64, gl: &[(f64, f64)]) -> (Acc, Acc) {
    let half = 0.5 * (b - a);
    let mut acc = [0.0; 11];
    let mut abs = [0.0; 11];
    for &(x, w) in gl {
        let v = g(a + half * (x + 1.0));
        add(&mut acc, &v, w * half);
        add(&mut abs, &v.map(f64::abs), w * half.abs());
    }
    (acc, abs)
}

fn face_integrals(
    b: &Brep,
    f: &Face,
    origin: DVec3,
    relative: f64,
    floor: &Acc,
    gl_edge: &[(f64, f64)],
    gl_inner: &[(f64, f64)],
) -> Quad {
    let s = &f.surface;
    let inner_u = !matches!(s, Surface::Extrusion { .. });
    // Reference line: middle of the boundary's parameter range.
    let (mut lo, mut hi) = (DVec2::splat(f64::MAX), DVec2::splat(f64::MIN));
    for l in &f.loops {
        for c in &l.coedges {
            let e = b.edge(c.edge);
            for j in 1..8 {
                let (u, v) = s.project(e.curve.point(e.t0 + (e.t1 - e.t0) * j as f64 / 8.0));
                lo = lo.min(DVec2::new(u, v));
                hi = hi.max(DVec2::new(u, v));
            }
        }
    }
    let mid = (lo + hi) * 0.5;
    let ctx = Ctx {
        surface: s,
        origin,
        inner_u,
        reference: if inner_u { mid.x } else { mid.y },
        inner: gl_inner.to_vec(),
    };
    let mut acc = Quad {
        value: [0.0; 11],
        converged: true,
    };
    for l in &f.loops {
        let n = l.coedges.len();
        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        for c in &l.coedges {
            let e = b.edge(c.edge);
            let eps = (e.t1 - e.t0) * 1e-6;
            let (ta, tb, ia, ib) = if c.forward {
                (e.t0, e.t1, e.t0 + eps, e.t1 - eps)
            } else {
                (e.t1, e.t0, e.t1 - eps, e.t0 + eps)
            };
            starts.push(end_uv(s, e.curve.point(ta), e.curve.point(ia)));
            ends.push(end_uv(s, e.curve.point(tb), e.curve.point(ib)));
            // Boundary integral along the edge, adaptively on each smooth
            // piece; the floor is shared out by parameter length.
            let dir = if c.forward { 1.0 } else { -1.0 };
            let g = |t: f64| -> Acc {
                let (u, v) = s.project(e.curve.point(t));
                let uv = DVec2::new(u, v);
                let duv = param_velocity(s, uv, e.curve.deriv(t)) * dir;
                let (mut val, d) = ctx.boundary(uv, duv);
                for x in &mut val {
                    *x *= d;
                }
                val
            };
            let span = (e.t1 - e.t0).abs();
            let quad = Adaptive {
                g: &g,
                gl: gl_edge,
                relative,
                floor: floor.map(|x| x / span),
            };
            for w in e.curve.smooth_breaks(e.t0, e.t1).windows(2) {
                let (whole, _) = rule(&g, w[0], w[1], gl_edge);
                quad.run(w[0], w[1], &whole, 0, &mut acc);
            }
        }
        // Parameter-space connectors at singular points (zero length in 3D).
        for k in 0..n {
            let a = ends[k];
            let bb = starts[(k + 1) % n];
            let d = bb - a;
            if d.length() > 1e-12 {
                for &(x, wt) in gl_edge {
                    let uv = a + d * (0.5 * (x + 1.0));
                    let (val, dd) = ctx.boundary(uv, d);
                    add(&mut acc.value, &val, dd * wt * 0.5);
                }
            }
        }
    }
    // The loop traversal already carries the face sense for flux terms; the
    // area integrand is unsigned, so restore its sign.
    if !f.forward {
        acc.value[0] = -acc.value[0];
    }
    acc
}

/// Flux integrals through the thin strip between an edge's projections onto
/// the surfaces of its two faces (design note 0003, "Tolerant edges").
///
/// Each face integral follows the edge as projected onto its own surface,
/// so where an edge lies on its surfaces only within the tolerance (a
/// vertex merged from end points that miss each other, a tolerant
/// intersection curve), the two faces trace two slightly different
/// curves and the boundary is open along a strip between them. The flux
/// through an open surface depends on the origin of the integrands, which
/// made the volume of such a body wrong by about the strip's width times
/// its length times the distance to the origin. Closing it with the strip
/// leaves only errors of the order of the width squared.
///
/// `forward` is the direction in which `face`'s loop runs along the edge;
/// `other` uses it the other way. The strip is ruled between the two
/// projections; its flux is taken at the edge curve with the edge's
/// tangent, which differs from the exact ruled strip's by terms of the
/// order of the width squared. The area integrand is left out: the strip
/// is no part of the surface.
#[allow(clippy::too_many_arguments)]
fn strip_integrals(
    b: &Brep,
    edge: usize,
    face: &Face,
    forward: bool,
    other: &Face,
    origin: DVec3,
    relative: f64,
    floor: &Acc,
    gl_edge: &[(f64, f64)],
) -> Quad {
    let e = &b.edges[edge];
    let (s1, s2) = (&face.surface, &other.surface);
    let on = |s: &Surface, p: DVec3| {
        let (u, v) = s.project(p);
        s.point(u, v)
    };
    let width = |t: f64| {
        let p = e.curve.point(t);
        on(s2, p) - on(s1, p)
    };
    let mut acc = Quad {
        value: [0.0; 11],
        converged: true,
    };
    // Edges on both surfaces to rounding (all exact geometry) need nothing.
    let scale = 1.0 + (e.curve.point(e.t0) - origin).length();
    if (0..=8).all(|j| {
        let t = e.t0 + (e.t1 - e.t0) * j as f64 / 8.0;
        width(t).length() <= 1e-14 * scale
    }) {
        return acc;
    }
    let sign = if forward { -1.0 } else { 1.0 };
    let g = |t: f64| -> Acc {
        let mut val = integrand_at(e.curve.point(t) - origin, e.curve.deriv(t), width(t));
        val[0] = 0.0;
        val.map(|x| x * sign)
    };
    let span = (e.t1 - e.t0).abs();
    let quad = Adaptive {
        g: &g,
        gl: gl_edge,
        relative,
        floor: floor.map(|x| x / span),
    };
    for w in e.curve.smooth_breaks(e.t0, e.t1).windows(2) {
        let (whole, _) = rule(&g, w[0], w[1], gl_edge);
        quad.run(w[0], w[1], &whole, 0, &mut acc);
    }
    acc
}

/// Mass properties of a body with unit density: each integral to within
/// `relative_accuracy` (no finer than [`NOISE`]) of the boundary integral
/// of its integrand's magnitude, as estimated by [`Adaptive`], plus a
/// rounding floor of [`NOISE`] times the body's `size` to the integral's
/// dimension (area L², volume L³, first moments L⁴, second moments L⁵).
/// Fails if some piece of the boundary cannot be integrated that
/// accurately.
pub(crate) fn mass_properties(
    b: &Brep,
    size: f64,
    relative_accuracy: f64,
) -> Result<MassProperties, KernelError> {
    let gl_edge = gauss_legendre(8);
    let gl_inner = gauss_legendre(16);
    // Integrate about the vertex centroid for conditioning.
    let origin = if b.vertices.is_empty() {
        DVec3::ZERO
    } else {
        b.vertices.iter().map(|v| v.point).sum::<DVec3>() / b.vertices.len() as f64
    };
    let coedges: usize = b
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .map(|l| l.coedges.len())
        .sum();
    let share = NOISE / coedges.max(1) as f64;
    let dims = [2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5];
    let floor: Acc = dims.map(|d| share * size.powi(d));
    let relative = relative_accuracy.max(NOISE);
    let mut acc = [0.0; 11];
    let refuse = || {
        KernelError::MassProperties(format!(
            "the boundary integrals did not converge to a relative accuracy of \
             {relative_accuracy:e}"
        ))
    };
    // The faces using each edge, with the direction of their loops.
    let mut uses: Vec<Vec<(usize, bool)>> = vec![Vec::new(); b.edges.len()];
    for (fi, f) in b.faces.iter().enumerate() {
        let fa = face_integrals(b, f, origin, relative, &floor, &gl_edge, &gl_inner);
        if !fa.converged {
            return Err(refuse());
        }
        add(&mut acc, &fa.value, 1.0);
        for l in &f.loops {
            for c in &l.coedges {
                uses[c.edge.idx()].push((fi, c.forward));
            }
        }
    }
    for (e, u) in uses.iter().enumerate() {
        if let [(f1, fwd), (f2, _)] = u[..] {
            if f1 == f2 {
                continue;
            }
            let st = strip_integrals(
                b,
                e,
                &b.faces[f1],
                fwd,
                &b.faces[f2],
                origin,
                relative,
                &floor,
                &gl_edge,
            );
            if !st.converged {
                return Err(refuse());
            }
            add(&mut acc, &st.value, 1.0);
        }
    }
    let volume = acc[1];
    let first = DVec3::new(acc[2], acc[3], acc[4]);
    let d = if volume != 0.0 {
        first / volume
    } else {
        DVec3::ZERO
    };
    let cov_o = DMat3::from_cols(
        DVec3::new(acc[5], acc[8], acc[10]),
        DVec3::new(acc[8], acc[6], acc[9]),
        DVec3::new(acc[10], acc[9], acc[7]),
    );
    // Parallel-axis shift of the second moments to the centroid.
    let cov_c = cov_o - DMat3::from_cols(d * d.x, d * d.y, d * d.z) * volume;
    let tr = cov_c.x_axis.x + cov_c.y_axis.y + cov_c.z_axis.z;
    let inertia = DMat3::from_diagonal(DVec3::splat(tr)) - cov_c;
    Ok(MassProperties {
        volume,
        area: acc[0],
        centroid: origin + d,
        inertia: inertia.to_cols_array_2d(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_legendre_integrates_polynomials_exactly() {
        let g = gauss_legendre(8);
        // ∫_{-1}^{1} x^14 dx = 2/15
        let s: f64 = g.iter().map(|(x, w)| w * x.powi(14)).sum();
        assert!((s - 2.0 / 15.0).abs() < 1e-14);
        let w: f64 = g.iter().map(|(_, w)| w).sum();
        assert!((w - 2.0).abs() < 1e-14);
    }

    fn integrate(g: &dyn Fn(f64) -> Acc, relative: f64) -> Quad {
        let gl = gauss_legendre(8);
        let quad = Adaptive {
            g,
            gl: &gl,
            relative,
            floor: [0.0; 11],
        };
        let mut out = Quad {
            value: [0.0; 11],
            converged: true,
        };
        let (whole, _) = rule(g, 0.0, 1.0, &gl);
        quad.run(0.0, 1.0, &whole, 0, &mut out);
        out
    }

    #[test]
    fn a_sharp_peak_is_resolved_however_the_first_rules_agree() {
        // ∫₀¹ w/(1 + w²(t − ½)²) dt = 2·atan(w/2): almost all of it within
        // 1/w of the middle, where the first rules barely look.
        for w in [10.0, 1e3, 1e5] {
            let g = |t: f64| {
                let mut v = [0.0; 11];
                let x = w * (t - 0.5);
                v[1] = w / (1.0 + x * x);
                v
            };
            let q = integrate(&g, 1e-12);
            assert!(q.converged);
            let exact = 2.0 * (w / 2.0).atan();
            assert!(
                (q.value[1] - exact).abs() < 1e-12 * exact,
                "{w}: {}",
                q.value[1]
            );
        }
    }

    #[test]
    fn an_integrand_it_cannot_resolve_is_reported() {
        // Evaluation noise of 10 % (a hash of the parameter's bits): no
        // piece ever passes a relative test, so the result is refused
        // rather than returned.
        let g = |t: f64| {
            let mut v = [0.0; 11];
            let h = t.to_bits().wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 11;
            v[1] = 1.0 + 0.1 * (h as f64 / (1u64 << 53) as f64 - 0.5);
            v
        };
        assert!(!integrate(&g, 1e-9).converged);
    }
}
