//! Tests of the rounding allowances in courses' bends ([`Trace::bend`],
//! design note 0012 §16.2) against references in double-double arithmetic
//! ([`Dd`]): the chart points of conic tracks in planes' and cones'
//! charts, the hulls of NURBS curves' clipped pieces in a plane's chart,
//! and NURBS courses in cylinders' and cones' charts measured from a point
//! far along the axis. Each checks the claim the planar graph rests on: the exact
//! course of a stretch lies within the band, the sag's bound plus the
//! trace's slack, of the line through the stretch's computed ends (the
//! computed ends are shared by neighbouring chords, so the course's part
//! beyond an end, displaced along the course by rounding, lies in the
//! next chord's band).
//!
//! The first two are adapted from the measurements that set the
//! allowances (2026-10-15: 1.4 million conic chart points, ten million
//! NURBS stretches). Larger samples than the defaults run with
//! `CONIC_ROUNDING_CASES`, `HULL_ROUNDING_CURVES` and
//! `DEVELOPED_ROUNDING_CURVES`.

use std::f64::consts::FRAC_PI_2;

use glam::{DVec2, DVec3, DVec4};
use opendatum_nurbs::NurbsCurve;

use super::super::dd::Dd;
use super::super::meet::{meet, Meet, Region, Track};
use super::super::planar::Bend;
use super::super::quadric::Quadric;
use super::super::testing::Rng;
use super::{trace, Chart, Trace, TraceOf};
use crate::geom::Curve3;
use crate::math::Frame;

const EPS: f64 = f64::EPSILON;

/// The number of cases from the environment variable `name`, or the
/// default: `release` in release builds, a thirtieth of it in debug ones.
fn cases(name: &str, release: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(if cfg!(debug_assertions) {
            release / 30
        } else {
            release
        })
}

/// `10^x`, `x` uniform in `[a, b)`.
fn log(rng: &mut Rng, a: f64, b: f64) -> f64 {
    10f64.powf(rng.range(a, b))
}

/// `1` or `−1`, evenly.
fn sign(rng: &mut Rng) -> f64 {
    if rng.next() < 0.5 {
        -1.0
    } else {
        1.0
    }
}

/// A random unit vector square to the unit `z`.
fn across(rng: &mut Rng, z: DVec3) -> DVec3 {
    loop {
        let v = rng.unit();
        let w = v - z * v.dot(z);
        if w.length() > 0.1 {
            return w.normalize();
        }
    }
}

type V = [Dd; 3];

fn v(a: DVec3) -> V {
    [Dd::of(a.x), Dd::of(a.y), Dd::of(a.z)]
}

fn v_add(a: V, b: V) -> V {
    [a[0].add(b[0]), a[1].add(b[1]), a[2].add(b[2])]
}

fn v_sub(a: V, b: V) -> V {
    [a[0].sub(b[0]), a[1].sub(b[1]), a[2].sub(b[2])]
}

fn v_scale(a: V, k: Dd) -> V {
    [a[0].mul(k), a[1].mul(k), a[2].mul(k)]
}

fn v_dot(a: V, b: V) -> Dd {
    a[0].mul(b[0]).add(a[1].mul(b[1])).add(a[2].mul(b[2]))
}

fn v_f(a: V) -> DVec3 {
    DVec3::new(a[0].f(), a[1].f(), a[2].f())
}

/// [`Chart::map`] of a plane's, a cylinder's or a cone's chart in
/// double-double: the same formula, with the chart's stored numbers taken
/// as exact.
fn exact_map(chart: &Chart, q: V) -> [Dd; 2] {
    match *chart {
        Chart::Plane(f) => {
            let d = v_sub(q, v(f.origin));
            [v_dot(d, v(f.x)), v_dot(d, v(f.y))]
        }
        Chart::Cylinder {
            point,
            axis,
            x,
            y,
            radius,
        } => {
            let w = v_sub(q, v(point));
            [
                Dd::atan2(v_dot(w, v(y)), v_dot(w, v(x))).mulf(radius),
                v_dot(w, v(axis)),
            ]
        }
        Chart::Cone {
            point,
            axis,
            x,
            y,
            radius,
            sin,
            cos,
        } => {
            let w = v_sub(q, v(point));
            let h = v_dot(w, v(axis));
            let across = v_sub(w, v_scale(v(axis), h));
            let l = v_dot(across, across).sqrt();
            let sigma = l.sub(Dd::of(radius)).mulf(sin).add(h.mulf(cos));
            let psi = Dd::atan2(v_dot(across, v(y)), v_dot(across, v(x))).mulf(sin);
            let s0 = Dd::of(radius).div(Dd::of(sin));
            let (sp, cp) = psi.sin_cos();
            let (half, _) = psi.mulf(0.5).sin_cos();
            [
                sigma.mul(cp).sub(s0.mulf(2.0).mul(half).mul(half)),
                s0.add(sigma).mul(sp),
            ]
        }
        _ => unreachable!(),
    }
}

/// The distance of `y` from the line through `a` and `b`.
fn from_line(a: DVec2, b: DVec2, y: [Dd; 2]) -> f64 {
    let (dx, dy) = (Dd::of(b.x).sub(Dd::of(a.x)), Dd::of(b.y).sub(Dd::of(a.y)));
    let len = dx.mul(dx).add(dy.mul(dy)).sqrt();
    let (ux, uy) = (y[0].sub(Dd::of(a.x)), y[1].sub(Dd::of(a.y)));
    dx.mul(uy).sub(dy.mul(ux)).div(len).f().abs()
}

/// The most of a ratio, with a description of the case attaining it.
#[derive(Default)]
struct Most {
    value: f64,
    case: String,
}

impl Most {
    fn see(&mut self, value: f64, case: impl FnOnce() -> String) {
        if value > self.value {
            self.value = value;
            self.case = case();
        }
    }
}

// Conic tracks.

/// The point of a conic track at `t`, exactly as stored (to about 10⁻³⁰
/// relative), and its `w = √(1 − c t²)`; `None` where `1 − c t² ≤ 0`.
fn conic_point(track: &Track, t: f64) -> Option<(V, Dd)> {
    let Track::Conic { vertex, x, y, p, c } = *track else {
        unreachable!()
    };
    let t = Dd::of(t);
    let tt = t.mul(t);
    let w2 = Dd::ONE.sub(Dd::of(c).mul(tt));
    if w2.0 <= 0.0 {
        return None;
    }
    let w = w2.sqrt();
    let xi = Dd::of(p).mul(tt).div(Dd::ONE.add(w));
    let eta = Dd::of(p).mul(t);
    Some((
        v_add(v_add(v(vertex), v_scale(v(x), xi)), v_scale(v(y), eta)),
        w,
    ))
}

/// The angle about a cone chart's axis of point `q`, from its `x` toward
/// its `y`; zero in a plane's chart.
fn chart_angle(chart: &Chart, q: V) -> f64 {
    let Chart::Cone {
        point, axis, x, y, ..
    } = *chart
    else {
        return 0.0;
    };
    let w = v_sub(q, v(point));
    let h = v_dot(w, v(axis));
    let across = v_f(v_sub(w, v_scale(v(axis), h)));
    across.dot(y).atan2(across.dot(x))
}

/// The conic kinds sampled.
const KINDS: [&str; 3] = ["hyperbola", "parabola", "half ellipse"];

/// A plane across a cone and their conic, from [`meet`].
struct Section {
    plane: Quadric,
    cone: Quadric,
    /// The cone again, with the reference circle a chart is anchored on.
    chart_cone: Quadric,
    track: Track,
    /// The cone's frame: its axis `a3`, and `a1` toward the vertex.
    a1: DVec3,
    a2: DVec3,
}

/// A cone and a plane meeting in a conic of kind `kind` with its vertex
/// at `at`, the cone's reference circle for charts `extra` along the
/// cone's ruling beyond the vertex's parallel; a half ellipse reaching
/// about `scale` from its vertex.
///
/// About the apex, with the axis `a3` and `a1` toward the vertex: the
/// cone of half-angle α and the plane with normal `a1 sin β + a3 cos β`,
/// `d0 = p / tan α` from the apex, meet in the conic with semi-latus
/// rectum `p`, `c = 1 − e² = cos(α + β) cos(α − β) / cos² α` (`e = sin β /
/// cos α`), its vertex `s_v = d0 / cos(β − α)` from the apex along the
/// ruling toward `a1`, opening along `−a1 cos β + a3 sin β` (§14).
fn section(rng: &mut Rng, kind: usize, scale: f64, at: DVec3, extra: f64) -> Section {
    let alpha = [0.01, 0.1, 0.4, 0.8, 1.2, 1.45][rng.below(6)] * rng.range(0.9, 1.03);
    let (sa, ca) = alpha.sin_cos();
    let k = sa / ca;
    let p = log(rng, -3.0, 1.0);
    let c = match kind {
        0 => -k * k * log(rng, -8.0, 0.0),
        1 => 0.0,
        _ => (p / scale * rng.range(0.3, 1.0)).min(0.9),
    };
    // cos² β = sin² α + c cos² α.
    let beta = (sa * sa + c * ca * ca).max(0.0).sqrt().min(1.0).acos();
    let (sb, cb) = beta.sin_cos();
    let d0 = p / k;
    let s_v = d0 / (beta - alpha).cos();
    let a3 = rng.unit();
    let a1 = across(rng, a3);
    let a2 = a3.cross(a1);
    let apex = at - (a1 * (s_v * sa) + a3 * (s_v * ca));
    let normal = a1 * sb + a3 * cb;
    let plane = Quadric::Plane {
        point: apex + normal * d0,
        normal,
    };
    // The cone with its reference circle through the vertex, for meet
    // (whose vertex loses precision as the square of the reference's
    // radius when it is far from the section; the kernel keeps it near the
    // faces), and the same cone with its reference `extra` farther along
    // the ruling, for the chart.
    let at_height = |h0: f64| Quadric::Cone {
        point: apex + a3 * h0,
        axis: a3,
        radius: h0 * k,
        slope: k,
    };
    let cone = at_height(s_v * ca);
    let chart_cone = at_height((s_v + extra) * ca);
    // A region at the vertex, small beside a half ellipse (whose centre,
    // p / c along the axis, must lie well outside it for a half).
    let r = if kind == 2 {
        (0.01 * p).min(p / c / 64.0)
    } else {
        0.01 * p
    };
    let region = Region::of_box(at - DVec3::splat(r), at + DVec3::splat(r));
    let m = meet(&plane, &cone, &region, 1e-7);
    let Meet::Curves(tracks) = &m else {
        panic!("{} α {alpha} p {p} c {c}: {m:?}", KINDS[kind]);
    };
    let near = |t: &Track| match *t {
        Track::Conic { vertex, .. } => (vertex - at).length(),
        _ => f64::INFINITY,
    };
    let mut track = *tracks
        .iter()
        .min_by(|a, b| near(a).total_cmp(&near(b)))
        .expect("a conic");
    let Track::Conic {
        vertex,
        x,
        p: pm,
        c: cm,
        ..
    } = &mut track
    else {
        panic!("{} α {alpha} p {p} c {c}: {m:?}", KINDS[kind]);
    };
    let x_expected = -a1 * cb + a3 * sb;
    let size = 1.0 + s_v + at.length();
    assert!(
        (*vertex - at).length() <= 1e-9 * size
            && (*pm - p).abs() <= 1e-8 * p
            && (*cm - c).abs() <= 1e-9 * c.abs().max(1.0)
            && (*x - x_expected).length() <= 1e-6,
        "{}: α {alpha} p {p} c {c}: {track:?} vertex {at} x {x_expected}",
        KINDS[kind]
    );
    if kind == 1 {
        // The parabola exactly (the section's c is 0 within rounding).
        *cm = 0.0;
    }
    Section {
        plane,
        cone,
        chart_cone,
        track,
        a1,
        a2,
    }
}

/// A chart of the section's cone from its reference circle for charts,
/// its middle at angle `theta` about the axis from `a1` toward `a2`.
fn cone_chart(rng: &mut Rng, s: &Section, theta: f64) -> Chart {
    let Quadric::Cone {
        point,
        axis,
        radius,
        slope,
    } = s.chart_cone
    else {
        unreachable!()
    };
    let mid = s.a1 * theta.cos() + s.a2 * theta.sin();
    let x = (mid - axis * mid.dot(axis)).normalize();
    let y = x.cross(axis);
    let y = if rng.next() < 0.5 { y } else { -y };
    let cos = 1.0 / (1.0 + slope * slope).sqrt();
    Chart::Cone {
        point,
        axis,
        x,
        y,
        radius,
        sin: slope * cos,
        cos,
    }
}

/// The distance of exact point `q` from the section's plane and cone,
/// over the size of its offset from the cone's reference circle.
fn off_surfaces(s: &Section, q: V) -> f64 {
    let (
        Quadric::Plane { point: pp, normal },
        Quadric::Cone {
            point,
            axis,
            radius,
            slope,
        },
    ) = (s.plane, s.cone)
    else {
        unreachable!()
    };
    let to_plane = v_dot(v_sub(q, v(pp)), v(normal)).f().abs();
    let w = v_sub(q, v(point));
    let h = v_dot(w, v(axis));
    let across = v_sub(w, v_scale(v(axis), h));
    let l = v_dot(across, across).sqrt();
    let to_cone =
        l.sub(Dd::of(radius).add(Dd::of(slope).mul(h))).f().abs() / (1.0 + slope * slope).sqrt();
    let size = v_f(w).length() + radius;
    to_plane.max(to_cone) / size
}

/// The parameter of the point of the conic at distance `reach` from its
/// vertex (on the side of positive `t`), or `None` beyond half an
/// ellipse: `ξ² + η² = reach²` with `η² = 2pξ − cξ²`.
fn t_of_reach(p: f64, c: f64, reach: f64) -> Option<f64> {
    let xi = reach * reach / (p + (p * p + (1.0 - c) * reach * reach).sqrt());
    if c > 0.0 && xi >= p / c {
        return None;
    }
    let eta2 = (2.0 * p * xi - c * xi * xi).max(0.0);
    Some(eta2.sqrt() / p)
}

/// Conic tracks' courses in planes' and cones' charts lie within their
/// bands: hyperbolas, parabolas and half ellipses (to `w` under 10⁻⁸ at
/// a half's ends, taken by the parameter's last bits; checked below) cut
/// from real cones by [`meet`], stretches from 1 m to 10⁷ m from their
/// vertices whose bend is from a thousandth of `ε` times that reach to a
/// hundred times it (where the rounding is most of the band), the vertex
/// near the world's origin and the chart's anchor (mode 0), the anchor
/// up to 10⁴ m from the vertex (mode 1), or the vertex up to 10⁴ m from
/// the world's origin (mode 2). The exact course at seventeen points of
/// each stretch is checked against the band about the line through the
/// stretch's computed ends.
///
/// Without the bend's allowance for the rounding that grows with the
/// reach ([`super::CONIC_ROUNDING`]), or without the chart's size in the slack,
/// courses leave their bands (both checked below); the allowance is at
/// least the first-order bound derived for a cone's chart, 27.2 ε reach,
/// and the measured rounding is at most 1.96 ε (reach + V), V the
/// vertex's part (its coordinates, its distance from the anchor and a
/// cone's radius there).
#[test]
fn conic_courses_lie_within_their_bands() {
    let mut rng = Rng(0xC071_C50A);
    let cases = cases("CONIC_ROUNDING_CASES", 36_000);
    let mut chords = 0usize;
    let mut off = 0.0f64;
    // The least `w` of a checked stretch on a half ellipse.
    let mut least_w = f64::INFINITY;
    // The most a course's distance from its chord is of the band; of what
    // the band would be without the bend's allowance; and without the
    // chart's part of the slack.
    let mut most = [Most::default(), Most::default(), Most::default()];
    for case in 0..cases {
        let chart_kind = case % 2;
        let kind = (case / 2) % 3;
        let mode = (case / 6) % 3;
        let scale = log(&mut rng, 0.0, 7.0);
        let at = match mode {
            2 => rng.unit() * log(&mut rng, 0.0, 4.0),
            _ => rng.unit() * log(&mut rng, -3.0, 0.0),
        };
        let extra = if mode == 1 {
            log(&mut rng, 0.0, 4.0)
        } else {
            0.0
        };
        let s = section(&mut rng, kind, scale, at, extra);
        let Track::Conic { vertex, p, c, .. } = s.track else {
            unreachable!()
        };
        // The chart: a plane's frame with its origin near the vertex (or
        // up to 10⁴ m from it in mode 1), or a cone's chart from the
        // reference circle for charts, its middle near the vertex's
        // ruling.
        let chart = if chart_kind == 0 {
            let Quadric::Plane { normal, .. } = s.plane else {
                unreachable!()
            };
            let z = normal * sign(&mut rng);
            let fx = across(&mut rng, z);
            let fy = z.cross(fx);
            let offset = if mode == 1 {
                log(&mut rng, 0.0, 4.0)
            } else {
                log(&mut rng, -6.0, -3.0)
            };
            let dir = (fx * rng.range(-1.0, 1.0) + fy * rng.range(-1.0, 1.0)).normalize();
            Chart::Plane(Frame {
                origin: vertex + dir * offset,
                x: fx,
                y: fy,
                z,
            })
        } else {
            let theta = rng.range(-0.3, 0.3);
            cone_chart(&mut rng, &s, theta)
        };
        let tr = trace(&chart, TraceOf::Track(s.track))
            .expect("traced")
            .expect("curved");
        // Points at reaches from 1 m to 10⁷ m (a half ellipse to its
        // end's) and near the ends of a half ellipse (w from 10⁻⁸ to 1).
        let mut ts = Vec::new();
        let most_reach = if c > 0.0 {
            (p / c) * (1.0 + c).sqrt()
        } else {
            f64::INFINITY
        };
        for _ in 0..8 {
            let reach = log(&mut rng, 0.0, most_reach.min(1e7).log10().max(0.0));
            if let Some(t) = t_of_reach(p, c, reach) {
                ts.push(t * sign(&mut rng));
            }
        }
        if c > 0.0 {
            for _ in 0..4 {
                let w = log(&mut rng, -8.0, 0.0);
                // Back to where the conic is defined, by the parameter's
                // last bits at the half's end.
                let mut t = ((1.0 - w * w) / c).sqrt();
                while conic_point(&s.track, t).is_none() {
                    t = t.next_down();
                }
                ts.push(t * sign(&mut rng));
            }
        }
        for _ in 0..6 {
            let t = ts[rng.below(ts.len())];
            let Some((q, w)) = conic_point(&s.track, t) else {
                continue;
            };
            let reach = v_f(v_sub(q, v(vertex))).length();
            let w = w.f();
            // The stretch from t toward the vertex on a half ellipse, away
            // from it otherwise, its bend p / w³ (t1 − t0)² / 8 about
            // `bend` (w least at t).
            let bend = log(&mut rng, -3.0, 2.0) * EPS * reach.max(1.0);
            let dt = (8.0 * bend * w * w * w / p).sqrt() * t.signum();
            // Near a half's end that stretch is under the parameter's last
            // bit: one to eight units of it instead.
            let ulp = t.abs().next_up() - t.abs();
            let dt = if dt.abs() < ulp {
                (1 + rng.below(8)) as f64 * ulp * t.signum()
            } else {
                dt
            };
            let (t0, t1) = if c > 0.0 { (t - dt, t) } else { (t, t + dt) };
            let (Some((e0, _)), Some((e1, _))) =
                (conic_point(&s.track, t0), conic_point(&s.track, t1))
            else {
                continue;
            };
            if chart_angle(&chart, e0).abs() > FRAC_PI_2
                || chart_angle(&chart, e1).abs() > FRAC_PI_2
            {
                continue;
            }
            let (a, b) = (tr.at(t0), tr.at(t1));
            if t0 == t1 || a == b {
                continue;
            }
            chords += 1;
            if c > 0.0 {
                least_w = least_w.min(w);
            }
            let sag = tr.sag(t0, t1);
            let band = sag.bound + tr.slack();
            let unallowed = sag.bound - sag.floor + tr.slack();
            let unanchored = sag.bound + 1e-13 * (1.0 + vertex.abs().max_element());
            for j in 0..=16 {
                let tj = t0 + (t1 - t0) * j as f64 / 16.0;
                let Some((qj, _)) = conic_point(&s.track, tj) else {
                    continue;
                };
                off = off.max(off_surfaces(&s, qj));
                let d = from_line(a, b, exact_map(&chart, qj));
                let describe = || {
                    format!(
                        "case {case}: {} {} mode {mode}, p {p:.3e} c {c:.3e} t [{t0:?}, {t1:?}] \
                         w {w:.2e} reach {reach:.3e} vertex {vertex} chart {chart:?}",
                        ["plane", "cone"][chart_kind],
                        KINDS[kind]
                    )
                };
                assert!(d <= band, "{}: {d:e} beyond the band {band:e}", describe());
                most[0].see(d / band, describe);
                most[1].see(d / unallowed, describe);
                most[2].see(d / unanchored, describe);
            }
        }
    }
    eprintln!(
        "{chords} chords; the course over its band, at most {:.3} ({}); over the band without \
         the allowance {:.3} ({}); without the chart's part of the slack {:.3} ({})",
        most[0].value, most[0].case, most[1].value, most[1].case, most[2].value, most[2].case
    );
    assert!(off <= 1e-9, "points off the cone or the plane: {off:e}");
    // Stretches reach the ends of half ellipses, where `w` is least
    // (2.9·10⁻⁹ in the default run of a debug build, 9.3·10⁻¹⁰ in a
    // release one).
    eprintln!("least w on a half ellipse: {least_w:e}");
    assert!(least_w < 1e-8, "{least_w:e}");
    // Both parts of the band are needed: without either, some course
    // leaves its band.
    assert!(
        most[1].value > 1.0 && most[2].value > 1.0,
        "{:?}",
        [most[1].value, most[2].value]
    );
}

// NURBS curves in a plane's chart.

/// A homogeneous point in double-double.
type H = [Dd; 4];

fn h(q: DVec4) -> H {
    [Dd::of(q.x), Dd::of(q.y), Dd::of(q.z), Dd::of(q.w)]
}

fn comb(a: &H, b: &H, u: Dd, s: Dd) -> H {
    std::array::from_fn(|k| a[k].mul(u).add(b[k].mul(s)))
}

/// De Casteljau's split at `s`.
fn split_dd(pw: &[H], s: Dd) -> (Vec<H>, Vec<H>) {
    let n = pw.len();
    let u = Dd::ONE.sub(s);
    let mut b = pw.to_vec();
    let mut left = vec![b[0]];
    let mut right = vec![[Dd::ZERO; 4]; n];
    right[n - 1] = b[n - 1];
    for r in 1..n {
        for i in 0..n - r {
            b[i] = comb(&b[i], &b[i + 1], u, s);
        }
        left.push(b[0]);
        right[n - 1 - r] = b[n - 1 - r];
    }
    (left, right)
}

fn point_dd(q: &H) -> V {
    [q[0].div(q[3]), q[1].div(q[3]), q[2].div(q[3])]
}

/// The distance of `q` from the segment from `a` to `b`, as
/// [`super::super::planar::point_segment`] takes it.
fn point_segment_dd(q: [Dd; 2], a: [Dd; 2], b: [Dd; 2]) -> Dd {
    let d = [b[0].sub(a[0]), b[1].sub(a[1])];
    let u = [q[0].sub(a[0]), q[1].sub(a[1])];
    let len2 = d[0].mul(d[0]).add(d[1].mul(d[1]));
    if len2.0 == 0.0 {
        return u[0].mul(u[0]).add(u[1].mul(u[1])).sqrt();
    }
    let mut t = u[0].mul(d[0]).add(u[1].mul(d[1])).div(len2);
    if t.lt(Dd::ZERO) {
        t = Dd::ZERO;
    } else if Dd::ONE.lt(t) {
        t = Dd::ONE;
    }
    let e = [u[0].sub(d[0].mul(t)), u[1].sub(d[1].mul(t))];
    e[0].mul(e[0]).add(e[1].mul(e[1])).sqrt()
}

/// A Bézier piece of a NURBS curve, extracted in double-double.
struct Piece {
    a: f64,
    b: f64,
    points: Vec<H>,
}

/// `NurbsCurve::bezier_segments` in double-double: the blossoms of the
/// curve's control points at each span's ends (Ramshaw 1987).
fn pieces_dd(n: &NurbsCurve) -> Vec<Piece> {
    let p = n.degree();
    let u = n.knots().as_slice();
    let pts = n.homogeneous_points();
    let (lo, hi) = n.domain();
    let mut out = Vec::new();
    for k in p..pts.len() {
        let (a, b) = (u[k], u[k + 1]);
        if b <= a || b <= lo || a >= hi {
            continue;
        }
        let points = (0..=p)
            .map(|j| {
                let mut args = vec![a; p - j];
                args.extend(std::iter::repeat_n(b, j));
                let mut d: Vec<H> = pts[k - p..=k].iter().map(|&q| h(q)).collect();
                for (r, &x) in args.iter().enumerate() {
                    let r = r + 1;
                    for jj in (r..=p).rev() {
                        let i = k - p + jj;
                        let alpha = Dd::of(x)
                            .sub(Dd::of(u[i]))
                            .div(Dd::of(u[i + p + 1 - r]).sub(Dd::of(u[i])));
                        d[jj] = comb(&d[jj - 1], &d[jj], Dd::ONE.sub(alpha), alpha);
                    }
                }
                d[p]
            })
            .collect();
        out.push(Piece { a, b, points });
    }
    out
}

/// The curve's point at `t`, from its exact pieces.
fn curve_dd(pieces: &[Piece], t: f64) -> V {
    let piece = pieces
        .iter()
        .find(|s| s.a <= t && t <= s.b)
        .expect("t in the domain");
    let s = Dd::of(t)
        .sub(Dd::of(piece.a))
        .div(Dd::of(piece.b).sub(Dd::of(piece.a)));
    let u = Dd::ONE.sub(s);
    let mut b = piece.points.clone();
    let n = b.len();
    for r in 1..n {
        for i in 0..n - r {
            b[i] = comb(&b[i], &b[i + 1], u, s);
        }
    }
    point_dd(&b[0])
}

/// The exact control points of the curve's pieces on `[t0, t1]`.
fn clipped_dd(pieces: &[Piece], t0: f64, t1: f64) -> Vec<H> {
    pieces
        .iter()
        .filter(|s| s.b > t0 && s.a < t1)
        .flat_map(|s| {
            let len = Dd::of(s.b).sub(Dd::of(s.a));
            let sa = Dd::of(s.a.max(t0)).sub(Dd::of(s.a)).div(len);
            let sb = Dd::of(s.b.min(t1)).sub(Dd::of(s.a)).div(len);
            let left = if sb.lt(Dd::ONE) {
                split_dd(&s.points, sb).0
            } else {
                s.points.clone()
            };
            if Dd::ZERO.lt(sa) {
                split_dd(&left, sa.div(sb)).1
            } else {
                left
            }
        })
        .collect()
}

/// The rational Bézier `pw` raised by one degree, rounded to `f64`.
fn elevate(pw: &[DVec4]) -> Vec<DVec4> {
    let p = pw.len() - 1;
    let hp: Vec<H> = pw.iter().map(|&q| h(q)).collect();
    (0..=p + 1)
        .map(|i| {
            let a = Dd::of(i as f64).div(Dd::of((p + 1) as f64));
            let q = if i == 0 {
                hp[0]
            } else if i == p + 1 {
                hp[p]
            } else {
                comb(&hp[i - 1], &hp[i], a, Dd::ONE.sub(a))
            };
            DVec4::new(q[0].f(), q[1].f(), q[2].f(), q[3].f())
        })
        .collect()
}

/// `n` weights log-uniform between 1 and `rho`, with 1 and `rho` taken.
fn weights(rng: &mut Rng, n: usize, rho: f64) -> Vec<f64> {
    let mut w: Vec<f64> = (0..n).map(|_| log(rng, 0.0, rho.log10())).collect();
    let i = rng.below(n);
    let j = (i + 1 + rng.below(n - 1)) % n;
    w[i] = 1.0;
    w[j] = rho;
    w
}

/// A random rational curve in the plane through `centre` with unit
/// normal `z`: a conic arc (kind 0), one raised to degree 3 to 5 (1), a
/// Bézier piece of degree 2 to 5 (2), or a curve of two to four pieces
/// (3), its weights' ratio up to 10⁴ and its control points within
/// `size` of `centre`; with a description.
fn rational(
    rng: &mut Rng,
    kind: usize,
    centre: DVec3,
    z: DVec3,
    size: f64,
) -> (NurbsCurve, String) {
    let p = match kind {
        0 => 2,
        1 => 3 + rng.below(3),
        _ => 2 + rng.below(4),
    };
    let fx = across(rng, z);
    let fy = z.cross(fx);
    let in_plane =
        |rng: &mut Rng| centre + fx * rng.range(-size, size) + fy * rng.range(-size, size);
    let lift = |x: DVec3, w: f64| (x * w).extend(w);
    // One span's knots: at 0 or anywhere to ±1000, of length 1 or 10⁻³
    // to 10³.
    let span = |rng: &mut Rng, p: usize| {
        let k0 = if rng.below(2) == 0 {
            0.0
        } else {
            rng.range(-1e3, 1e3)
        };
        let len = if rng.below(2) == 0 {
            1.0
        } else {
            log(rng, -3.0, 3.0)
        };
        let mut k = vec![k0; p + 1];
        k.extend(vec![k0 + len; p + 1]);
        k
    };
    let (knots, pw, label) = match kind {
        0 | 1 => {
            let w = log(rng, -4.0, 4.0);
            let mut pw = vec![
                lift(in_plane(rng), 1.0),
                lift(in_plane(rng), w),
                lift(in_plane(rng), 1.0),
            ];
            while pw.len() < p + 1 {
                pw = elevate(&pw);
            }
            (span(rng, p), pw, format!("conic w {w:.4e}"))
        }
        2 => {
            let rho = log(rng, 0.0, 4.0);
            let ws = weights(rng, p + 1, rho);
            let pw = ws.iter().map(|&w| lift(in_plane(rng), w)).collect();
            (span(rng, p), pw, format!("Bézier w {ws:?}"))
        }
        _ => {
            let spans = 2 + rng.below(3);
            let rho = log(rng, 0.0, 4.0);
            let ws = weights(rng, p + spans, rho);
            let pw = ws.iter().map(|&w| lift(in_plane(rng), w)).collect();
            let k0 = if rng.below(2) == 0 {
                0.0
            } else {
                log(rng, 0.0, 3.0)
            };
            let scale = log(rng, -2.0, 2.0);
            let mut knots = vec![k0; p + 1];
            let mut t = k0;
            for _ in 1..spans {
                t += scale * log(rng, -3.0, 0.0);
                knots.push(t);
            }
            t += scale * log(rng, -3.0, 0.0);
            knots.extend(vec![t; p + 1]);
            let label = format!("{spans} spans w {ws:?} knots {knots:?}");
            (knots, pw, label)
        }
    };
    let n = NurbsCurve::from_homogeneous(p, knots, pw).expect("a valid curve");
    (n, format!("{label}, degree {p}"))
}

/// A stretch of the piece on `[a, b]`: from its whole down to 10⁻¹² of
/// it, at its ends, near them or anywhere inside.
fn stretch(rng: &mut Rng, a: f64, b: f64) -> (f64, f64) {
    let len = b - a;
    let f = if rng.below(8) == 0 {
        1.0
    } else {
        log(rng, -12.0, 0.0)
    };
    let (t0, t1) = match rng.below(5) {
        0 => (a, a + f * len),
        1 => (b - f * len, b),
        2 => {
            let t0 = a + log(rng, -12.0, -1.0) * len;
            (t0, t0 + f * len)
        }
        3 => {
            let t1 = b - log(rng, -12.0, -1.0) * len;
            (t1 - f * len, t1)
        }
        _ => {
            let t0 = a + rng.range(0.0, 1.0 - f) * len;
            (t0, t0 + f * len)
        }
    };
    (t0.max(a), t1.min(b))
}

/// The distances a NURBS curve's clipped control points are computed to
/// have from a chord in a plane's chart ([`Trace::bend`]) differ from the
/// exact ones (the curve's pieces, the clips, the chord's ends and the
/// map, all exact, from the stored control points) by at most the
/// trace's allowance for rounding ([`super::hull_rounding`]), and the
/// exact course lies within its band: rational curves of degree 2 to 5,
/// conic arcs, raised conics, Bézier pieces and curves of several pieces
/// with weights' ratios up to 10⁴, in planes 1 m to 10¹⁰ m from the
/// origin, the chart's origin near the curve, at the foot of the plane's
/// normal through the origin, or as far from the curve as it is from the
/// origin, stretches down to 10⁻¹² of a piece.
///
/// Where the weights are uneven and the curve is large, the rounding of
/// the clip's parameters moves the clipped points by up to `ε (1 + ρ) D`
/// (`D` the control points' diameter, `ρ` the weights' ratio), more than
/// an allowance without the weights could hold (checked below).
#[test]
fn hull_distances_round_within_their_allowance() {
    let mut rng = Rng(0x0B5E_7AB1_E5EE_D004);
    let curves = cases("HULL_ROUNDING_CURVES", 60_000);
    let mut most = Most::default();
    let mut most_band = Most::default();
    let mut most_eps = Most::default();
    let mut most_unweighted = Most::default();
    let mut stretches = 0usize;
    for case in 0..curves {
        let kind = case % 4;
        let far = log(&mut rng, 0.0, 10.0);
        let centre = rng.unit() * far;
        let size = if rng.below(4) == 0 {
            far * log(&mut rng, -1.0, 0.3)
        } else {
            log(&mut rng, -3.0, 1.0)
        };
        let z = rng.unit();
        let (n, label) = rational(&mut rng, kind, centre, z, size);
        let fx = across(&mut rng, z);
        let origin = match rng.below(3) {
            0 => centre + (fx * rng.range(-1.0, 1.0) + z.cross(fx) * rng.range(-1.0, 1.0)) * size,
            1 => z * z.dot(centre),
            _ => centre + (fx * rng.range(-1.0, 1.0) + z.cross(fx) * rng.range(-1.0, 1.0)) * far,
        };
        let chart = Chart::Plane(Frame {
            origin,
            x: fx,
            y: z.cross(fx),
            z,
        });
        let curve = Curve3::Nurbs(n.clone());
        let tr: Trace = trace(&chart, TraceOf::Curve(&curve))
            .expect("traced")
            .expect("curved");
        let m = 1.0 + curve.magnitude();
        // The allowance without its weights' term (`hull_rounding`'s other
        // terms, the box's diagonal for the diameter).
        let points = n.control_points();
        let (lo, hi) = points.iter().fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(lo, hi), q| (lo.min(*q), hi.max(*q)),
        );
        let (r, d) = points.iter().fold((0.0f64, 0.0f64), |(r, d), q| {
            (r.max(q.length()), d.max((*q - origin).length()))
        });
        let p = n.degree() as f64;
        let unweighted = EPS * ((7.0 * p + 2.0) * r + 6.0 * p * (hi - lo).length() + 6.4 * d);
        let segments = n.bezier_segments();
        let exact = pieces_dd(&n);
        for _ in 0..8 {
            let j = rng.below(segments.len());
            let (t0, t1) = stretch(&mut rng, segments[j].a, segments[j].b);
            if t1 <= t0 {
                continue;
            }
            stretches += 1;
            let sag = tr.sag(t0, t1);
            let hull = sag.bound - sag.floor;
            let (a, b) = (tr.at(t0), tr.at(t1));
            let map = |q: V| exact_map(&chart, q);
            let (ea, eb) = (map(curve_dd(&exact, t0)), map(curve_dd(&exact, t1)));
            let truth = clipped_dd(&exact, t0, t1)
                .iter()
                .map(|q| point_segment_dd(map(point_dd(q)), ea, eb))
                .fold(Dd::ZERO, Dd::max);
            let error = Dd::of(hull).sub(truth).f().abs();
            let describe = || {
                format!(
                    "case {case}: {label}, far {far:.3e} size {size:.3e} origin {origin} \
                     t [{t0:?}, {t1:?}] in [{:?}, {:?}]",
                    segments[j].a, segments[j].b
                )
            };
            assert!(
                error <= sag.floor,
                "{}: the hull distance {hull:e} is {error:e} off, over its allowance {:e}",
                describe(),
                sag.floor
            );
            most.see(error / sag.floor, describe);
            most_eps.see(error / (EPS * m), describe);
            most_unweighted.see(error / unweighted, describe);
            if a != b {
                let band = sag.bound + tr.slack();
                for k in 0..=8 {
                    // At most `t1`: rounding could carry the last past the
                    // domain's end.
                    let t = (t0 + (t1 - t0) * k as f64 / 8.0).min(t1);
                    let d = from_line(a, b, map(curve_dd(&exact, t)));
                    assert!(d <= band, "{}: {d:e} beyond the band {band:e}", describe());
                    most_band.see(d / band, describe);
                }
            }
        }
    }
    eprintln!(
        "{stretches} stretches: the hull's rounding at most {:.3} of its allowance ({}); \
         {:.1} ε (1 + magnitude) ({}); {:.1} times the allowance without the weights ({}); \
         the course at most {:.3} of its band ({})",
        most.value,
        most.case,
        most_eps.value,
        most_eps.case,
        most_unweighted.value,
        most_unweighted.case,
        most_band.value,
        most_band.case
    );
    // The weights' term of the allowance is needed: without it, the
    // rounding of some hull distances would exceed the allowance.
    assert!(
        most_unweighted.value > 1.0,
        "{} ({})",
        most_unweighted.value,
        most_unweighted.case
    );
}

// NURBS courses in cylinders' and cones' charts.

/// A rational quadratic arc of `delta` about `center`, its middle toward
/// `map(1, 0)`: the image of a circle's arc by `map`, reweighted (1,
/// `w k`, `k²`), the same curve followed unevenly.
fn arc(center: DVec3, map: impl Fn(f64, f64) -> DVec3, delta: f64, k: f64) -> NurbsCurve {
    let w = (0.5 * delta).cos();
    let half = 0.5 * delta;
    let points = vec![
        center + map(half.cos(), -half.sin()),
        center + map(1.0 / w, 0.0),
        center + map(half.cos(), half.sin()),
    ];
    NurbsCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        points,
        Some(vec![1.0, w * k, k * k]),
    )
    .expect("a valid arc")
}

/// NURBS courses in cylinders' and cones' charts lie within their bands
/// (the bound from the course's derivatives, [`Trace::bend`], and the
/// slack): ellipses cut from cylinders by planes leaning up to 2 along the
/// axis per radius, and circles about cones' axes (half-angles 0.008 to
/// 1.56), as rational quadratic arcs of 0.2 to 2.5 radians, reweighted;
/// the curves near the world's origin or up to 10⁴ m from it, the chart's
/// point on the axis near them or up to 10⁶ m from them along it;
/// stretches down to 10⁻¹² of the parameter. The exact course at nine
/// points of each stretch is checked against the band about the line
/// through the stretch's computed ends.
///
/// The chart's points are computed from their offsets from the chart's
/// point, so their rounding grows with its distance: without the chart's
/// size in the slack, courses leave their bands (checked below).
#[test]
fn developed_courses_lie_within_their_bands() {
    let mut rng = Rng(0xDE7E_10ED_C0DE_0005);
    let curves = cases("DEVELOPED_ROUNDING_CURVES", 20_000);
    let mut most = Most::default();
    let mut most_without = Most::default();
    let mut stretches = 0usize;
    for case in 0..curves {
        let cone = case % 2 == 1;
        let axis = rng.unit();
        let x = across(&mut rng, axis);
        let centre = if rng.below(2) == 0 {
            rng.unit() * log(&mut rng, -3.0, 0.0)
        } else {
            rng.unit() * log(&mut rng, 0.0, 4.0)
        };
        let delta = rng.range(0.2, 2.5);
        let k = log(&mut rng, -1.0, 1.0);
        let near = rng.below(4) == 0;
        let (chart, n, label) = if cone {
            let y = x.cross(axis);
            let alpha = [0.01, 0.2, 0.7, 1.3][rng.below(4)] * rng.range(0.8, 1.2);
            let (sin, cos) = alpha.sin_cos();
            let slope = sin / cos;
            let r0 = log(&mut rng, -2.0, 1.0);
            // A circle about the axis at height h from the centre's.
            let h = rng.range(-0.5, 0.5) * r0 / slope.max(0.1);
            let r = r0 + slope * h;
            let n = arc(
                centre + axis * h,
                |u, v| x * (r * u) + y * (r * v),
                delta,
                k,
            );
            // The chart's point away from the apex, its radius growing.
            let along = if near {
                rng.range(0.0, 1.0) * r0
            } else {
                log(&mut rng, 0.0, 6.0)
            };
            let chart = Chart::Cone {
                point: centre + axis * along,
                axis,
                x,
                y,
                radius: r0 + slope * along,
                sin,
                cos,
            };
            (
                chart,
                n,
                format!("cone α {alpha:.3} radius {r0:.3e} h {h:.3e}"),
            )
        } else {
            let y = axis.cross(x);
            let radius = log(&mut rng, -2.0, 1.0);
            let lean = rng.range(-2.0, 2.0);
            let n = arc(
                centre,
                |u, v| x * (radius * u) + y * (radius * v) + axis * (lean * radius * u),
                delta,
                k,
            );
            let along = if near {
                rng.range(-1.0, 1.0)
            } else {
                log(&mut rng, 0.0, 6.0) * sign(&mut rng)
            };
            let chart = Chart::Cylinder {
                point: centre + axis * along,
                axis,
                x,
                y,
                radius,
            };
            (
                chart,
                n,
                format!("cylinder radius {radius:.3e} lean {lean:.3}"),
            )
        };
        let curve = Curve3::Nurbs(n.clone());
        let tr: Trace = trace(&chart, TraceOf::Curve(&curve))
            .expect("traced")
            .expect("curved");
        let exact = pieces_dd(&n);
        let m = 1.0 + curve.magnitude();
        for _ in 0..8 {
            let (t0, t1) = stretch(&mut rng, 0.0, 1.0);
            if t1 <= t0 {
                continue;
            }
            let sag = tr.sag(t0, t1);
            let (a, b) = (tr.at(t0), tr.at(t1));
            if !sag.bound.is_finite() || a == b {
                continue;
            }
            stretches += 1;
            let band = sag.bound + tr.slack();
            let without = sag.bound + 1e-13 * m;
            for j in 0..=8 {
                // At most `t1`: rounding could carry the last past it.
                let t = (t0 + (t1 - t0) * j as f64 / 8.0).min(t1);
                let d = from_line(a, b, exact_map(&chart, curve_dd(&exact, t)));
                let describe = || {
                    format!(
                        "case {case}: {label}, centre {centre}, chart {chart:?}, \
                         t [{t0:?}, {t1:?}]"
                    )
                };
                assert!(d <= band, "{}: {d:e} beyond the band {band:e}", describe());
                most.see(d / band, describe);
                most_without.see(d / without, describe);
            }
        }
    }
    eprintln!(
        "{stretches} stretches: the course at most {:.3} of its band ({}); {:.1} of the band \
         without the chart's part of the slack ({})",
        most.value, most.case, most_without.value, most_without.case
    );
    // The chart's part of the slack is needed: without it, some course
    // leaves its band.
    assert!(
        most_without.value > 1.0,
        "{} ({})",
        most_without.value,
        most_without.case
    );
}
