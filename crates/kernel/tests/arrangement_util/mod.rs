//! Helpers shared by the arrangement tests: curve builders, an independent
//! curve evaluator, structural checks of arrangements, and end-to-end
//! extrusion of regions.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

pub mod oracle;
pub mod raster;

use opendatum_kernel::arrangement::region_containing;
use opendatum_kernel::{
    Body, BrepKernel, DVec2, ExtrudeSpec, Frame, Kernel, KernelError, Profile, ProfileCurve,
    SketchArrangement, Tolerance,
};

pub fn k() -> BrepKernel {
    BrepKernel::default()
}

pub fn tol() -> f64 {
    Tolerance::default().linear
}

pub fn v(x: f64, y: f64) -> DVec2 {
    DVec2::new(x, y)
}

pub fn line(a: DVec2, b: DVec2) -> ProfileCurve {
    ProfileCurve::Line { start: a, end: b }
}

pub fn circle(c: DVec2, r: f64) -> ProfileCurve {
    ProfileCurve::Circle {
        center: c,
        radius: r,
    }
}

pub fn arc(c: DVec2, r: f64, start: f64, sweep: f64) -> ProfileCurve {
    ProfileCurve::Arc {
        center: c,
        radius: r,
        start_angle: start,
        sweep_angle: sweep,
    }
}

pub fn polygon(p: &[DVec2]) -> Vec<ProfileCurve> {
    (0..p.len())
        .map(|i| line(p[i], p[(i + 1) % p.len()]))
        .collect()
}

pub fn assert_rel(a: f64, b: f64, rel: f64) {
    assert!(
        (a - b).abs() <= rel * b.abs().max(1e-300),
        "{a} vs {b} (relative error {:.3e})",
        ((a - b) / b).abs()
    );
}

pub fn sorted_areas(a: &SketchArrangement) -> Vec<f64> {
    let mut x: Vec<f64> = a.regions.iter().map(|r| r.area).collect();
    x.sort_by(f64::total_cmp);
    x
}

/// Rational de Boor evaluation, written independently of the kernel.
pub fn spline_point(
    degree: usize,
    knots: &[f64],
    points: &[DVec2],
    weights: Option<&[f64]>,
    t: f64,
) -> DVec2 {
    let n = points.len() - 1;
    let p = degree;
    let mut k = p;
    while k < n && knots[k + 1] <= t {
        k += 1;
    }
    let w = |i: usize| weights.map_or(1.0, |w| w[i]);
    let mut d: Vec<[f64; 3]> = (0..=p)
        .map(|j| {
            let i = j + k - p;
            [points[i].x * w(i), points[i].y * w(i), w(i)]
        })
        .collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = j + k - p;
            let den = knots[i + p - r + 1] - knots[i];
            let a = if den > 0.0 { (t - knots[i]) / den } else { 0.0 };
            let prev = d[j - 1];
            for (x, q) in d[j].iter_mut().zip(prev) {
                *x = (1.0 - a) * q + a * *x;
            }
        }
    }
    v(d[p][0] / d[p][2], d[p][1] / d[p][2])
}

/// Point of an input curve at a parameter in the arrangement conventions.
pub fn curve_point(c: &ProfileCurve, t: f64) -> DVec2 {
    match c {
        ProfileCurve::Line { start, end } => *start + (*end - *start) * t,
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let a = start_angle + sweep_angle.signum() * t;
            *center + *radius * v(a.cos(), a.sin())
        }
        ProfileCurve::Circle { center, radius } => *center + *radius * v(t.cos(), t.sin()),
        ProfileCurve::Spline {
            degree,
            knots,
            points,
            weights,
        } => spline_point(*degree, knots, points, weights.as_deref(), t),
    }
}

/// Domain of a curve's parameter.
pub fn domain(c: &ProfileCurve) -> (f64, f64) {
    match c {
        ProfileCurve::Line { .. } => (0.0, 1.0),
        ProfileCurve::Arc { sweep_angle, .. } => (0.0, sweep_angle.abs()),
        ProfileCurve::Circle { .. } => (0.0, std::f64::consts::TAU),
        ProfileCurve::Spline {
            degree,
            knots,
            points,
            ..
        } => (knots[*degree], knots[points.len()]),
    }
}

pub fn spline_samples(c: &ProfileCurve, n: usize) -> Vec<DVec2> {
    let (a, b) = domain(c);
    (0..=n)
        .map(|i| curve_point(c, a + (b - a) * i as f64 / n as f64))
        .collect()
}

/// Samples of an edge piece in the direction of travel.
pub fn piece_samples(c: &ProfileCurve, t0: f64, t1: f64, fwd: bool, n: usize) -> Vec<DVec2> {
    let n = if matches!(c, ProfileCurve::Line { .. }) {
        1
    } else {
        n
    };
    let mut s: Vec<DVec2> = (0..=n)
        .map(|i| curve_point(c, t0 + (t1 - t0) * i as f64 / n as f64))
        .collect();
    if !fwd {
        s.reverse();
    }
    s
}

/// Signed area of a polygon, relative to its first point (which keeps
/// rounding small for small polygons far from the origin).
pub fn shoelace(p: &[DVec2]) -> f64 {
    let Some(&o) = p.first() else {
        return 0.0;
    };
    (0..p.len())
        .map(|i| (p[i] - o).perp_dot(p[(i + 1) % p.len()] - o))
        .sum::<f64>()
        * 0.5
}

/// Polygon approximating a region loop.
pub fn loop_polygon(
    a: &SketchArrangement,
    curves: &[ProfileCurve],
    l: &opendatum_kernel::RegionLoop,
    n: usize,
) -> Vec<DVec2> {
    let mut out = Vec::new();
    for &(e, fwd) in &l.edges {
        let ed = &a.edges[e];
        // Every piece's own end points are kept, so tolerance-sized gaps
        // between pieces are closed by straight connectors.
        out.extend(piece_samples(&curves[ed.curve], ed.t0, ed.t1, fwd, n));
    }
    out
}

/// True if the boundary of the union of `sel` passes through some vertex
/// more than once (a non-manifold profile, which extrude rejects). Bridges
/// (edges with selected regions on both sides) are not boundary.
pub fn pinched_union(a: &SketchArrangement, sel: &[usize]) -> bool {
    let mut seen = HashSet::new();
    for &r in sel {
        let reg = &a.regions[r];
        for l in std::iter::once(&reg.outer).chain(&reg.holes) {
            for &(e, fwd) in &l.edges {
                let ed = &a.edges[e];
                let other = if fwd { ed.right } else { ed.left };
                if other.is_some_and(|o| sel.contains(&o)) {
                    continue;
                }
                if !seen.insert(if fwd { ed.start } else { ed.end }) {
                    return true;
                }
            }
        }
    }
    false
}

pub fn pinched(a: &SketchArrangement, r: usize) -> bool {
    pinched_union(a, &[r])
}

/// Structural and geometric invariants of an arrangement, checked with
/// independent evaluation.
pub fn check_invariants(curves: &[ProfileCurve], a: &SketchArrangement) {
    let t = 4.0 * tol();
    for (i, x) in a.vertices.iter().enumerate() {
        assert!(x.point.is_finite());
        assert_eq!(x.curves.len(), x.params.len(), "vertex {i}");
        assert!(
            x.curves.windows(2).all(|w| w[0] < w[1]),
            "vertex {i} curves unsorted"
        );
        for (&c, &p) in x.curves.iter().zip(&x.params) {
            let (lo, hi) = domain(&curves[c]);
            assert!(
                p >= lo && p <= hi,
                "vertex {i} param {p} outside [{lo}, {hi}]"
            );
            if matches!(curves[c], ProfileCurve::Circle { .. }) {
                assert!(p < hi);
            }
            let d = (curve_point(&curves[c], p) - x.point).length();
            assert!(d <= t, "vertex {i} is {d:e} from curve {c}");
        }
    }
    for (i, e) in a.edges.iter().enumerate() {
        assert!(e.curve < curves.len() && e.t0 < e.t1, "edge {i}");
        let c = &curves[e.curve];
        let (lo, hi) = domain(c);
        if matches!(c, ProfileCurve::Circle { .. }) {
            assert!(
                e.t0 >= 0.0 && e.t0 < hi && e.t1 - e.t0 <= hi + 1e-12,
                "edge {i} params"
            );
        } else {
            assert!(
                e.t0 >= lo && e.t1 <= hi,
                "edge {i} params outside the domain"
            );
        }
        assert!(
            (curve_point(c, e.t0) - a.vertices[e.start].point).length() <= t,
            "edge {i} start"
        );
        assert!(
            (curve_point(c, e.t1) - a.vertices[e.end].point).length() <= t,
            "edge {i} end"
        );
        // An edge with the same face on both sides is a bridge (for
        // example a segment joining two islands); it bounds nothing.
    }
    let mut seen: HashMap<(usize, bool), usize> = HashMap::new();
    for (r, reg) in a.regions.iter().enumerate() {
        assert!(reg.area > 0.0);
        // Polygon areas with 2048 and 4096 samples per piece: their error
        // falls with the square of the spacing, so the second is closer by
        // about four times its difference from the first.
        let mut total = 0.0;
        let mut finer = 0.0;
        for (li, l) in std::iter::once(&reg.outer).chain(&reg.holes).enumerate() {
            for (j, &(e, fwd)) in l.edges.iter().enumerate() {
                let ed = &a.edges[e];
                let (_, end) = if fwd {
                    (ed.start, ed.end)
                } else {
                    (ed.end, ed.start)
                };
                let (n, nf) = l.edges[(j + 1) % l.edges.len()];
                let ns = if nf { a.edges[n].start } else { a.edges[n].end };
                assert_eq!(end, ns, "region {r} loop {li} is not closed");
                let face = if fwd { ed.left } else { ed.right };
                assert_eq!(face, Some(r), "edge {e} side does not name region {r}");
                assert!(seen.insert((e, fwd), r).is_none(), "half-edge used twice");
            }
            let poly = loop_polygon(a, curves, l, 2048);
            let ar = shoelace(&poly);
            finer += shoelace(&loop_polygon(a, curves, l, 4096));
            if li == 0 {
                assert!(
                    ar > 0.0,
                    "outer loop of region {r} is not counter-clockwise"
                );
            } else {
                assert!(ar < 0.0, "hole of region {r} is not clockwise");
            }
            total += ar;
        }
        let extrapolated = finer + (finer - total) / 3.0;
        assert!(
            (extrapolated - reg.area).abs() <= 1e-5 * reg.area + 1e-12 + (finer - total).abs(),
            "region {r}: area {} vs polygon {extrapolated} ({total}, {finer})",
            reg.area
        );
        assert_eq!(
            region_containing(a, curves, reg.interior_point, &Tolerance::default()),
            Some(r),
            "interior point of region {r}"
        );
    }
    for (i, e) in a.edges.iter().enumerate() {
        if let Some(r) = e.left {
            assert_eq!(seen.get(&(i, true)), Some(&r));
        }
        if let Some(r) = e.right {
            assert_eq!(seen.get(&(i, false)), Some(&r));
        }
    }
}

pub fn arrange_checked(curves: &[ProfileCurve]) -> SketchArrangement {
    let a = k()
        .sketch_regions(curves)
        .unwrap_or_else(|e| panic!("arrangement failed: {e}"));
    check_invariants(curves, &a);
    a
}

/// No two edges join the same vertices along the same geometry.
pub fn assert_no_duplicate_edges(curves: &[ProfileCurve], a: &SketchArrangement) {
    for i in 0..a.edges.len() {
        for j in i + 1..a.edges.len() {
            let (p, q) = (&a.edges[i], &a.edges[j]);
            let same = (p.start, p.end) == (q.start, q.end) || (p.start, p.end) == (q.end, q.start);
            if !same {
                continue;
            }
            let mp = curve_point(&curves[p.curve], 0.5 * (p.t0 + p.t1));
            let mq = curve_point(&curves[q.curve], 0.5 * (q.t0 + q.t1));
            assert!(
                (mp - mq).length() > 10.0 * tol(),
                "edges {i} and {j} duplicate each other"
            );
        }
    }
}

pub fn extrude_profile(curves: &[ProfileCurve], h: f64) -> Body {
    try_extrude(curves, h).unwrap_or_else(|e| panic!("extrude failed: {e}"))
}

pub fn try_extrude(curves: &[ProfileCurve], h: f64) -> Result<Body, KernelError> {
    let kern = k();
    let profile = Profile {
        plane: Frame::WORLD,
        curves: curves.to_vec(),
    };
    let body = kern.extrude(
        &profile,
        &ExtrudeSpec {
            distance: h,
            symmetric: false,
        },
    )?;
    let issues = kern.validate(&body);
    assert!(issues.is_empty(), "invalid body: {issues:?}");
    Ok(body)
}

pub fn volume(b: &Body) -> f64 {
    k().mass_properties(b, 1e-10).unwrap().volume
}

/// Extrudes the union of `sel` and checks volume = area × height.
pub fn extrude_region(curves: &[ProfileCurve], a: &SketchArrangement, sel: &[usize]) {
    let prof = k().region_profile(a, curves, sel).unwrap();
    let area: f64 = sel.iter().map(|&r| a.regions[r].area).sum();
    let h = area.sqrt().max(1e-6);
    let body = extrude_profile(&prof, h);
    assert_rel(volume(&body), area * h, 1e-9);
}

/// Distance from `p` to a polyline.
fn to_polyline(p: DVec2, poly: &[DVec2]) -> f64 {
    poly.windows(2)
        .map(|w| {
            let ab = w[1] - w[0];
            let s = ((p - w[0]).dot(ab) / ab.length_squared().max(1e-300)).clamp(0.0, 1.0);
            (w[0] + ab * s - p).length()
        })
        .fold(f64::INFINITY, f64::min)
}

/// Smallest gap between the curves of a profile (sampled): between curves
/// that do not meet, from an end point of one curve to another curve that
/// it does not end on, at a close approach of two curves that meet
/// elsewhere (a local minimum of their distance, farther from the meeting
/// point than ten times its size), and the width of a sliver between two
/// curves joining the same two points.
pub fn profile_min_gap(prof: &[ProfileCurve]) -> f64 {
    let s: Vec<Vec<DVec2>> = prof
        .iter()
        .map(|c| {
            let (a, b) = domain(c);
            piece_samples(c, a, b, true, 256)
        })
        .collect();
    let ends: Vec<[DVec2; 2]> = s.iter().map(|x| [x[0], x[x.len() - 1]]).collect();
    let mut d = f64::INFINITY;
    for i in 0..s.len() {
        for j in 0..s.len() {
            if i == j {
                continue;
            }
            let shared: Vec<DVec2> = ends[i]
                .iter()
                .copied()
                .filter(|p| ends[j].iter().any(|q| (*p - *q).length() < 1e-6))
                .collect();
            let gaps: Vec<f64> = s[i].iter().map(|p| to_polyline(*p, &s[j])).collect();
            if shared.is_empty() {
                d = d.min(gaps.iter().copied().fold(f64::INFINITY, f64::min));
                continue;
            }
            for e in ends[i] {
                if !shared.iter().any(|q| (e - *q).length() < 1e-6) {
                    d = d.min(to_polyline(e, &s[j]));
                }
            }
            for k in 1..gaps.len().saturating_sub(1) {
                let far = shared
                    .iter()
                    .all(|q| (s[i][k] - *q).length() > 10.0 * gaps[k]);
                if far && gaps[k] < gaps[k - 1] && gaps[k] <= gaps[k + 1] {
                    d = d.min(gaps[k]);
                }
            }
            if shared.len() == 2 && (shared[0] - shared[1]).length() > 1e-6 {
                // Two curves joining the same two points: the width of the
                // sliver between them.
                d = d.min(gaps.iter().copied().fold(0.0, f64::max));
            }
        }
    }
    d
}

/// Smallest radius of curvature of `c` near its end at parameter `t`
/// (stepping towards `towards`), from circles through three nearby points.
fn end_radius(c: &ProfileCurve, t: f64, towards: f64) -> f64 {
    [0.01, 0.05, 0.1]
        .iter()
        .map(|f| {
            let h = (towards - t) * f;
            let (p, q, r) = (
                curve_point(c, t),
                curve_point(c, t + h),
                curve_point(c, t + 2.0 * h),
            );
            let cross = (q - p).perp_dot(r - p).abs();
            if cross <= 0.0 {
                f64::INFINITY
            } else {
                (q - p).length() * (r - q).length() * (r - p).length() / (2.0 * cross)
            }
        })
        .fold(f64::INFINITY, f64::min)
}

/// True if a profile comes back within `snap` of itself: two different
/// joints (start points of two different curves) within `snap` of each
/// other, or two curves closer than `snap` away from where they meet (see
/// [`profile_min_gap`]). The kernel's profile builder merges end points
/// closer than ten tolerances, including the ends of the pieces it cuts
/// arcs and splines into, so such a profile reads as a loop touching itself
/// and is rejected as branching ([`KernelError::InvalidInput`]); the
/// arrangement keeps the parts apart because they are more than the
/// tolerance apart.
pub fn profile_pinched_at(prof: &[ProfileCurve], snap: f64) -> bool {
    let starts: Vec<DVec2> = prof.iter().map(|c| curve_point(c, domain(c).0)).collect();
    (0..starts.len())
        .any(|i| (i + 1..starts.len()).any(|j| (starts[i] - starts[j]).length() <= snap))
        || profile_min_gap(prof) <= snap
}

/// Number of extrusions rejected because the region's boundary comes back
/// within the profile builder's resolution of itself (see
/// [`profile_pinched_at`]).
pub static EXTRUDE_PINCHED: AtomicUsize = AtomicUsize::new(0);

/// Number of slivers of perturbed curves whose volume was computed to 1e-6
/// because the integrals could not show 1e-10 (see
/// [`extrude_all_within`]).
pub static SLIVER_MASS: AtomicUsize = AtomicUsize::new(0);

/// Number of regions of perturbed curves whose profile was refused for a
/// loop or edges below the profile builder's resolution (see
/// [`extrude_all_within`]).
pub static PROFILE_BELOW_RESOLUTION: AtomicUsize = AtomicUsize::new(0);

/// Extrudes every region on its own: each must give a valid solid with
/// volume = area × height, unless it is pinched (non-manifold profile) or
/// pinched at the profile builder's resolution (see [`profile_pinched_at`])
/// and refused as invalid input. Returns the number extruded.
pub fn extrude_all(curves: &[ProfileCurve], a: &SketchArrangement) -> usize {
    extrude_all_within(curves, a, 0.0)
}

/// [`extrude_all`] for curves whose ends may miss each other by up to
/// `gap`: the area of a region is then defined only up to about
/// `gap` × perimeter, and the volume check allows for that; with a gap, a
/// region's profile may also be refused for a loop or edges below the
/// profile builder's resolution.
pub fn extrude_all_within(curves: &[ProfileCurve], a: &SketchArrangement, gap: f64) -> usize {
    let mut n = 0;
    for r in 0..a.regions.len() {
        let prof = match k().region_profile(a, curves, &[r]) {
            Ok(prof) => prof,
            // Curves that miss each other can leave loops and edges
            // between the tolerance and the profile builder's resolution
            // (ten tolerances): the tiny triangle three curves leave near a
            // point they would share, or the few microns between two
            // crossings of nearly tangent circles. The builder refuses
            // those by design (`region_profile_with_tolerance`).
            Err(KernelError::InvalidInput(m))
                if gap > 0.0 && m.contains("than the profile resolution") =>
            {
                let n = PROFILE_BELOW_RESOLUTION.fetch_add(1, Ordering::Relaxed) + 1;
                eprintln!("PROFILE_BELOW_RESOLUTION #{n}: {m}");
                continue;
            }
            Err(e) => panic!("region {r} has no profile: {e}"),
        };
        let area = a.regions[r].area;
        let h = area.sqrt().max(1e-6);
        match try_extrude(&prof, h) {
            // A non-manifold profile (the region touches itself at a
            // vertex) may be refused; if it is accepted, it is checked.
            Err(KernelError::InvalidInput(_)) if pinched(a, r) => continue,
            Ok(body) => {
                let reg = &a.regions[r];
                let perimeter: f64 = std::iter::once(&reg.outer)
                    .chain(&reg.holes)
                    .map(|l| {
                        let p = loop_polygon(a, curves, l, 64);
                        (0..p.len())
                            .map(|i| (p[(i + 1) % p.len()] - p[i]).length())
                            .sum::<f64>()
                    })
                    .sum();
                let got = match k().mass_properties(&body, 1e-10) {
                    Ok(m) => m.volume,
                    // A sliver between perturbed curves (a lens 0.3 µm wide
                    // and 3 mm long on an arc of 3 m radius) is evaluated
                    // with rounding of about 1e-16 of its coordinates, some
                    // 1e-9 of its width, while the kernel's noise floor is
                    // scaled to the body's size: its integrals cannot show
                    // 1e-10 and say so. The volume check below allows more
                    // than the sliver's whole volume for the gap anyway, so
                    // it is computed to 1e-6.
                    Err(KernelError::MassProperties(e))
                        if gap > 0.0 && area < 10.0 * gap * perimeter =>
                    {
                        let n = SLIVER_MASS.fetch_add(1, Ordering::Relaxed) + 1;
                        eprintln!("SLIVER_MASS #{n}: region {r}: {e}");
                        k().mass_properties(&body, 1e-6).unwrap().volume
                    }
                    Err(e) => panic!("region {r}: {e}"),
                };
                let want = area * h;
                assert!(
                    (got - want).abs() <= 1e-9 * want + 2.0 * gap * perimeter * h,
                    "region {r}: volume {got} vs area × height {want}"
                );
            }
            Err(e @ KernelError::InvalidInput(_)) if profile_pinched_at(&prof, 10.0 * tol()) => {
                let n = EXTRUDE_PINCHED.fetch_add(1, Ordering::Relaxed) + 1;
                eprintln!("EXTRUDE_PINCHED #{n}: {e}");
                continue;
            }
            Err(e) => panic!("region {r} does not extrude: {e}"),
        }
        n += 1;
    }
    n
}
