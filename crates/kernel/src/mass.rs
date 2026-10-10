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
//! depend on a tessellation; accuracy is limited only by quadrature (close to
//! machine precision for analytic geometry).

use glam::{DMat3, DVec2, DVec3};

use crate::geom::Surface;
use crate::topo::{Brep, Face};
use crate::MassProperties;

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
    let p = s.point(uv.x, uv.y) - o;
    let (su, sv) = s.derivs(uv.x, uv.y);
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

fn face_integrals(
    b: &Brep,
    f: &Face,
    origin: DVec3,
    chord: f64,
    gl_edge: &[(f64, f64)],
    gl_inner: &[(f64, f64)],
) -> Acc {
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
    let mut acc = [0.0; 11];
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
            // Boundary integral along the edge.
            let ts = e.curve.discretize(e.t0, e.t1, chord, 0.3);
            let dir = if c.forward { 1.0 } else { -1.0 };
            for w in ts.windows(2) {
                let half = 0.5 * (w[1] - w[0]);
                for &(x, wt) in gl_edge {
                    let t = w[0] + half * (x + 1.0);
                    let (u, v) = s.project(e.curve.point(t));
                    let uv = DVec2::new(u, v);
                    let duv = param_velocity(s, uv, e.curve.deriv(t)) * dir;
                    let (val, d) = ctx.boundary(uv, duv);
                    add(&mut acc, &val, d * wt * half);
                }
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
                    add(&mut acc, &val, dd * wt * 0.5);
                }
            }
        }
    }
    // The loop traversal already carries the face sense for flux terms; the
    // area integrand is unsigned, so restore its sign.
    if !f.forward {
        acc[0] = -acc[0];
    }
    acc
}

/// Mass properties of a body with unit density.
pub(crate) fn mass_properties(b: &Brep, chord: f64) -> MassProperties {
    let gl_edge = gauss_legendre(8);
    let gl_inner = gauss_legendre(16);
    // Integrate about the vertex centroid for conditioning.
    let origin = if b.vertices.is_empty() {
        DVec3::ZERO
    } else {
        b.vertices.iter().map(|v| v.point).sum::<DVec3>() / b.vertices.len() as f64
    };
    let mut acc = [0.0; 11];
    for f in &b.faces {
        let fa = face_integrals(b, f, origin, chord, &gl_edge, &gl_inner);
        add(&mut acc, &fa, 1.0);
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
    MassProperties {
        volume,
        area: acc[0],
        centroid: origin + d,
        inertia: inertia.to_cols_array_2d(),
    }
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
}
