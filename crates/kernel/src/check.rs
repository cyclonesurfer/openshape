//! B-rep validation, topology statistics and mass properties.

use glam::{DVec3, DVec4};
use rustc_hash::FxHashMap;

use opendatum_nurbs::NurbsCurve;

use crate::geom::{Curve3, Surface};
use crate::math::Frame;
use crate::tess::{self, TessellationOptions};
use crate::topo::{Brep, ShellId};
use crate::{KernelError, Tolerance, TopologyCounts};

/// A point between an edge's ends where it passes within `tol` of the axis
/// of a surface of revolution (the z axis of `fr`), where the surface is
/// singular: its angle is undefined there, so the edge has no parameters
/// on the face and integrals and tessellations along it break. Points
/// within `tol` of the edge's ends are its vertices, where an edge may
/// meet the axis.
///
/// A line's least distance from the axis is found in closed form
/// ([`line_through_axis`]). A NURBS edge each of whose spans lies, with the
/// convex hull of its control points, farther than `tol` from the axis on
/// one side of it, and a circle that does everywhere, are clear of it
/// ([`hull_clear_of_axis`], [`circle_clear_of_axis`]). Otherwise the
/// distance from the axis is sampled along the edge and each sampled local
/// least refined by golden-section search between its neighbours.
pub(crate) fn through_axis(c: &Curve3, t0: f64, t1: f64, fr: &Frame, tol: f64) -> Option<DVec3> {
    let clear = match c {
        Curve3::Nurbs(n) => hull_clear_of_axis(n, t0.min(t1), t0.max(t1), fr, tol),
        Curve3::Circle {
            center,
            x,
            y,
            radius,
        } => circle_clear_of_axis(*center, *x, *y, *radius, fr, tol),
        Curve3::Line { origin, dir } => return line_through_axis(*origin, *dir, t0, t1, fr, tol),
    };
    if clear {
        return None;
    }
    sampled_through_axis(c, t0, t1, fr, tol)
}

/// [`through_axis`] by samples alone.
fn sampled_through_axis(c: &Curve3, t0: f64, t1: f64, fr: &Frame, tol: f64) -> Option<DVec3> {
    const SAMPLES: usize = 64;
    let dist = |t: f64| (c.point(t) - fr.origin).cross(fr.z).length();
    let at = |k: usize| t0 + (t1 - t0) * k as f64 / SAMPLES as f64;
    let ds: Vec<f64> = (0..=SAMPLES).map(|k| dist(at(k))).collect();
    let (a, b) = (c.point(t0), c.point(t1));
    for k in 0..=SAMPLES {
        let left = k.saturating_sub(1);
        let right = (k + 1).min(SAMPLES);
        if ds[k] > ds[left] || ds[k] > ds[right] {
            continue;
        }
        // Golden-section search for the least distance on the bracket.
        let (mut lo, mut hi) = (at(left), at(right));
        let g = 0.5 * (5f64.sqrt() - 1.0);
        let (mut x1, mut x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
        let (mut f1, mut f2) = (dist(x1), dist(x2));
        for _ in 0..80 {
            if f1 <= f2 {
                hi = x2;
                (x2, f2) = (x1, f1);
                x1 = hi - g * (hi - lo);
                f1 = dist(x1);
            } else {
                lo = x1;
                (x1, f1) = (x2, f2);
                x2 = lo + g * (hi - lo);
                f2 = dist(x2);
            }
        }
        let t = if f1 <= f2 { x1 } else { x2 };
        let p = c.point(t);
        if dist(t) <= tol && (p - a).length() > tol && (p - b).length() > tol {
            return Some(p);
        }
    }
    None
}

/// [`through_axis`] for the line `origin + t dir`. The parts of its points
/// square to the axis are affine in `t`, so their length, the distance from
/// the axis, is convex: its one least on the edge is where that part is
/// square to the line's (clamped to the edge's ends), the point the sampled
/// search converges to. An edge whose ends both lie within `tol` of the
/// axis lies within it all along (convexity again): it runs along the axis
/// within the tolerance, wherever its distance is least (the sampled
/// search, deciding by rounding there, may find an end), and its middle
/// stands for it, reported when farther than `tol` from both ends.
fn line_through_axis(
    origin: DVec3,
    dir: DVec3,
    t0: f64,
    t1: f64,
    fr: &Frame,
    tol: f64,
) -> Option<DVec3> {
    let at = |t: f64| origin + dir * t;
    let dist = |t: f64| (at(t) - fr.origin).cross(fr.z).length();
    let away = |p: DVec3| (p - at(t0)).length() > tol && (p - at(t1)).length() > tol;
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    if dist(lo).max(dist(hi)) <= tol {
        let p = at(0.5 * (lo + hi));
        return away(p).then_some(p);
    }
    let perp = |v: DVec3| v - fr.z * v.dot(fr.z);
    let (u, d) = (perp(origin - fr.origin), perp(dir));
    let dd = d.length_squared();
    if dd == 0.0 {
        return None;
    }
    let t = (-u.dot(d) / dd).clamp(lo, hi);
    let p = at(t);
    (dist(t) <= tol && away(p)).then_some(p)
}

/// True if each span of `n` meeting `[lo, hi]` stays farther than `tol`
/// from the axis of `fr` (its z axis). A span lies in the convex hull of
/// its `p + 1` control points when their weights are positive (the strong
/// convex hull property; Piegl & Tiller, *The NURBS Book*, §4.2), and the
/// distance from the axis is at least the distance along any direction `d`
/// square to it, `(x − o) · d`, which is linear, so least over the hull at
/// one of the points. `d` points from the axis toward their mean; the
/// computed least loses a rounding allowance of 10⁻¹² of the coordinates'
/// size. A cheap sufficient test: where it fails, [`through_axis`] samples
/// the curve.
fn hull_clear_of_axis(n: &NurbsCurve, lo: f64, hi: f64, fr: &Frame, tol: f64) -> bool {
    let p = n.degree();
    let u = n.knots().as_slice();
    let pw = n.homogeneous_points();
    let off = |q: DVec4| {
        let v = q.truncate() / q.w - fr.origin;
        v - fr.z * v.dot(fr.z)
    };
    (p..pw.len())
        .filter(|&i| u[i] < u[i + 1] && u[i + 1] > lo && u[i] < hi)
        .all(|i| {
            let pts = &pw[i - p..=i];
            if pts.iter().any(|q| q.w.is_nan() || q.w <= 0.0) {
                return false;
            }
            let d = pts
                .iter()
                .fold(DVec3::ZERO, |m, &q| m + off(q))
                .normalize_or_zero();
            if d == DVec3::ZERO {
                return false;
            }
            let size = pts.iter().fold(fr.origin.length(), |m, q| {
                m.max((q.truncate() / q.w).length())
            });
            let least = pts.iter().fold(f64::INFINITY, |m, &q| m.min(off(q).dot(d)));
            least - 1e-12 * (1.0 + size) > tol
        })
}

/// True if the whole circle `center + radius (cos t x + sin t y)` stays
/// farther than `tol` from the axis of `fr`. With `u`, `X` and `Y` the parts
/// of `center − o`, `x` and `y` square to the axis, the distance from it is
/// `|u + r (cos t X + sin t Y)|`, and `|cos t X + sin t Y|` lies between the
/// singular values `σ₋ ≤ σ₊` of `[X Y]`, so the distance is at least
/// `max(r σ₋ − |u|, |u| − r σ₊)`: `r` for a circle about the axis. Less a
/// rounding allowance of 10⁻¹² of the coordinates' size. `σ₊` comes from
/// the larger eigenvalue of the Gram matrix and `σ₋` as `|X × Y| / σ₊`
/// (their product): the smaller eigenvalue's square root loses about
/// `√ε` to cancellation where `σ₋` is near zero, a circle in a plane along
/// the axis, which cleared such circles through it (found by review).
fn circle_clear_of_axis(
    center: DVec3,
    x: DVec3,
    y: DVec3,
    radius: f64,
    fr: &Frame,
    tol: f64,
) -> bool {
    let perp = |v: DVec3| v - fr.z * v.dot(fr.z);
    let (u, xa, ya) = (perp(center - fr.origin), perp(x), perp(y));
    let (a, b, c) = (xa.dot(xa), xa.dot(ya), ya.dot(ya));
    // The eigenvalues of the Gram matrix [[a, b], [b, c]] are σ².
    let mean = 0.5 * (a + c);
    let dev = (0.25 * (a - c) * (a - c) + b * b).sqrt();
    let most = (mean + dev).sqrt();
    let least = if most > 0.0 {
        xa.cross(ya).length() / most
    } else {
        0.0
    };
    let (r, ul) = (radius.abs(), u.length());
    let gap = (r * least - ul).max(ul - r * most);
    let size = 1.0 + center.length() + fr.origin.length() + r;
    gap - 1e-12 * size > tol
}

/// Counts entities and derives the genus from the Euler–Poincaré formula
/// `V − E + F − R = 2 (S − H)` where `R` is the number of inner loops.
pub(crate) fn counts(b: &Brep) -> TopologyCounts {
    let loops: usize = b.faces.iter().map(|f| f.loops.len()).sum();
    let rings = loops as i64 - b.faces.len() as i64;
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + b.faces.len() as i64 - rings;
    let s = b.shells.len() as i64;
    let twice_genus = 2 * s - chi;
    TopologyCounts {
        solids: b.solids.len(),
        shells: b.shells.len(),
        faces: b.faces.len(),
        loops,
        edges: b.edges.len(),
        vertices: b.vertices.len(),
        genus: if twice_genus >= 0 && twice_genus % 2 == 0 {
            Some((twice_genus / 2) as usize)
        } else {
            None
        },
    }
}

/// Geometry tolerance for "lies on": a generous multiple of the model
/// resolution to absorb projection round-off in iterative inversions.
pub(crate) fn on_tolerance(tol: &Tolerance) -> f64 {
    tol.linear * 100.0
}

/// Checks topological and geometric consistency. Returns a list of problems;
/// an empty list means the body is valid.
pub(crate) fn validate(b: &Brep, tol: &Tolerance) -> Vec<String> {
    let mut issues = Vec::new();
    let on_tol = on_tolerance(tol);
    // A plane needs no iteration (its projection is exact), so its edges
    // are held to what constructions leave: a boolean puts a vertex at a
    // cluster of points within twice the tolerance of each other, a few
    // tolerances at most off a face's plane (1.06 at most in the boolean
    // fuzz at twenty times its default size). Ten tolerances keeps a margin
    // and still catches a face built across a crease (design note 0012 §9).
    let plane_tol = tol.linear * 10.0;

    // Index bounds.
    let nv = b.vertices.len();
    let ne = b.edges.len();
    let nf = b.faces.len();
    for (i, e) in b.edges.iter().enumerate() {
        if e.v0.idx() >= nv || e.v1.idx() >= nv {
            issues.push(format!("edge {i}: vertex index out of range"));
        }
    }
    for (i, f) in b.faces.iter().enumerate() {
        for l in &f.loops {
            if l.coedges.iter().any(|c| c.edge.idx() >= ne) {
                issues.push(format!("face {i}: edge index out of range"));
            }
        }
    }
    for (i, s) in b.shells.iter().enumerate() {
        if s.faces.iter().any(|f| f.idx() >= nf) {
            issues.push(format!("shell {i}: face index out of range"));
        }
    }
    for (i, s) in b.solids.iter().enumerate() {
        if s.shells.iter().any(|h| h.idx() >= b.shells.len()) || s.shells.is_empty() {
            issues.push(format!("solid {i}: bad shell list"));
        }
    }
    if !issues.is_empty() {
        return issues;
    }

    // Edges: parameter range and end points.
    for (i, e) in b.edges.iter().enumerate() {
        if !e.t0.is_finite() || !e.t1.is_finite() || e.t1 <= e.t0 {
            issues.push(format!("edge {i}: empty or invalid parameter range"));
            continue;
        }
        if (e.curve.point(e.t0) - b.vertex(e.v0).point).length() > on_tol {
            issues.push(format!("edge {i}: start point is off its start vertex"));
        }
        if (e.curve.point(e.t1) - b.vertex(e.v1).point).length() > on_tol {
            issues.push(format!("edge {i}: end point is off its end vertex"));
        }
    }

    // Faces: loops closed, edges on surface.
    for (fi, f) in b.faces.iter().enumerate() {
        if f.loops.is_empty() {
            issues.push(format!("face {fi}: no loops"));
        }
        for (li, l) in f.loops.iter().enumerate() {
            if l.coedges.is_empty() {
                issues.push(format!("face {fi} loop {li}: empty"));
                continue;
            }
            for k in 0..l.coedges.len() {
                let a = l.coedges[k];
                let c = l.coedges[(k + 1) % l.coedges.len()];
                if b.coedge_end(a) != b.coedge_start(c) {
                    issues.push(format!("face {fi} loop {li}: not closed at coedge {k}"));
                }
            }
            let limit = if matches!(f.surface, Surface::Plane(_)) {
                plane_tol
            } else {
                on_tol
            };
            // A face on a surface of revolution lies within the half turn
            // of its frame, its angles in [0, π] (ADR 0008): on the side of
            // the frame's x–z plane that y points to.
            let half = match &f.surface {
                Surface::Revolution { frame, .. } => Some(frame),
                _ => None,
            };
            for c in &l.coedges {
                let e = b.edge(c.edge);
                for j in 0..=4 {
                    let t = e.t0 + (e.t1 - e.t0) * j as f64 / 4.0;
                    let p = e.curve.point(t);
                    let (u, v) = f.surface.project(p);
                    let d = (f.surface.point(u, v) - p).length();
                    if d > limit {
                        issues.push(format!(
                            "face {fi}: edge {} leaves the surface by {d:.3e}",
                            c.edge.0
                        ));
                        break;
                    }
                    if let Some(fr) = half {
                        let below = -(p - fr.origin).dot(fr.y);
                        if below > on_tol {
                            issues.push(format!(
                                "face {fi}: edge {} leaves the half turn of its surface by \
                                 {below:.3e} (ADR 0008)",
                                c.edge.0
                            ));
                            break;
                        }
                    }
                }
                if let Some(fr) = half {
                    if let Some(p) = through_axis(&e.curve, e.t0, e.t1, fr, on_tol) {
                        issues.push(format!(
                            "face {fi}: edge {} runs through the axis of its surface at {p:?} \
                             between its ends (a singular point of a face's surface is a vertex)",
                            c.edge.0
                        ));
                    }
                }
            }
        }
    }

    // Shells: each edge used exactly twice, in opposite directions.
    for (si, s) in b.shells.iter().enumerate() {
        let mut uses: FxHashMap<u32, (u32, u32)> = FxHashMap::default();
        for f in &s.faces {
            for l in &b.faces[f.idx()].loops {
                for c in &l.coedges {
                    let u = uses.entry(c.edge.0).or_default();
                    if c.forward {
                        u.0 += 1
                    } else {
                        u.1 += 1
                    }
                }
            }
        }
        for (e, (fw, bw)) in uses {
            if fw != 1 || bw != 1 {
                issues.push(format!(
                    "shell {si}: edge {e} used {fw} times forward and {bw} times reversed (not a closed manifold)"
                ));
            }
        }
    }

    if counts(b).genus.is_none() {
        issues.push("Euler-Poincaré check failed".into());
    }

    // Orientation: every outer shell must enclose positive volume.
    if issues.is_empty() {
        match shell_signed_volumes(b, tol) {
            Ok(volumes) => {
                for (si, solid) in b.solids.iter().enumerate() {
                    for (k, sh) in solid.shells.iter().enumerate() {
                        let v = volumes[sh.idx()];
                        if k == 0 && v <= 0.0 {
                            issues.push(format!("solid {si}: outer shell is inside-out"));
                        }
                        if k > 0 && v >= 0.0 {
                            issues.push(format!("solid {si}: void shell {k} is not inverted"));
                        }
                    }
                }
            }
            Err(e) => issues.push(format!("cannot measure the shells' volumes: {e}")),
        }
    }
    issues
}

pub(crate) fn bounding_size(b: &Brep) -> f64 {
    let mut lo = DVec3::splat(f64::MAX);
    let mut hi = DVec3::splat(f64::MIN);
    for v in &b.vertices {
        lo = lo.min(v.point);
        hi = hi.max(v.point);
    }
    for e in &b.edges {
        for j in 0..=8 {
            let p = e.curve.point(e.t0 + (e.t1 - e.t0) * j as f64 / 8.0);
            lo = lo.min(p);
            hi = hi.max(p);
        }
    }
    if lo.x > hi.x {
        0.0
    } else {
        (hi - lo).length()
    }
}

/// Each shell's volume (by index) from a coarse mesh, a hundredth of the
/// body's size from its surfaces, with a bound on the volume's error:
/// twice the area of the mesh of its curved faces times that distance (the
/// mesh of a flat face lies in its plane, so the region between the mesh
/// and the body is along curved faces), plus the rounding of the sum.
pub(crate) fn shell_volumes(b: &Brep, tol: &Tolerance) -> Result<Vec<(f64, f64)>, KernelError> {
    let opts = TessellationOptions {
        chord_tolerance: bounding_size(b).max(1e-9) * 1e-2,
        max_angle: 0.6,
    };
    let t = tess::tessellate(b, &opts, tol)?;
    Ok(b.shells
        .iter()
        .map(|shell| {
            let (mut v, mut curved, mut magnitude) = (0.0, 0.0, 0.0);
            for f in &shell.faces {
                let ft = &t.faces[f.idx()];
                let flat = b.faces[f.idx()].surface.plane().is_some();
                for tri in &ft.triangles {
                    let [a, bb, c] = tri.map(|i| ft.positions[i as usize]);
                    let w = a.dot(bb.cross(c)) / 6.0;
                    v += w;
                    magnitude += w.abs();
                    if !flat {
                        curved += 0.5 * (bb - a).cross(c - a).length();
                    }
                }
            }
            let rounding = 64.0 * f64::EPSILON * magnitude;
            (v, 2.0 * curved * opts.chord_tolerance + rounding)
        })
        .collect())
}

/// Each shell's signed volume (by index), as orientation is decided: from
/// [`shell_volumes`]' coarse mesh where its error bound cannot reach the
/// sign, else from the mass properties of the shell alone (exact on the
/// surfaces, to a millionth). A thin shell between two curved faces close
/// together (the ring a ball keeps outside a cone that crosses its sphere
/// twice close together) can mesh to nothing at the coarse mesh's
/// distance from the surfaces (found 2026-10-13: such a ring was refused
/// as inside out).
pub(crate) fn shell_signed_volumes(b: &Brep, tol: &Tolerance) -> Result<Vec<f64>, KernelError> {
    let size = bounding_size(b).max(tol.linear);
    shell_volumes(b, tol)?
        .into_iter()
        .enumerate()
        .map(|(i, (v, bound))| {
            if v.abs() > bound {
                return Ok(v);
            }
            let shell = shell_brep(b, ShellId(i as u32));
            Ok(crate::mass::mass_properties(&shell, size, 1e-6)?.volume)
        })
        .collect()
}

/// The shell `shell` of `b` as a body of its own (no solid), sharing all of
/// `b`'s vertices and edges.
pub(crate) fn shell_brep(b: &Brep, shell: ShellId) -> Brep {
    let mut out = Brep {
        vertices: b.vertices.clone(),
        edges: b.edges.clone(),
        ..Brep::default()
    };
    let ids = b.shells[shell.idx()]
        .faces
        .iter()
        .map(|f| {
            let face = &b.faces[f.idx()];
            out.add_face(face.surface.clone(), face.forward, face.loops.clone())
        })
        .collect();
    out.add_shell(ids);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Frame;
    use crate::{BrepKernel, Kernel};

    /// A face's plane moved off its edges: within ten tolerances it is
    /// still valid, beyond them not (the hundred tolerances curved
    /// surfaces get would pass it).
    #[test]
    fn planes_are_held_to_ten_tolerances() {
        let tol = Tolerance::default();
        let body = BrepKernel::default()
            .make_box(&Frame::WORLD, DVec3::ONE)
            .expect("a box");
        for (lift, valid) in [(5.0, true), (20.0, false)] {
            let mut b = (*body.brep).clone();
            let mut moved = 0;
            for f in &mut b.faces {
                if let Surface::Plane(fr) = &mut f.surface {
                    if fr.z.z.abs() > 0.5 && (fr.origin.z - 1.0).abs() < 1e-12 {
                        fr.origin += fr.z * lift * tol.linear;
                        moved += 1;
                    }
                }
            }
            assert_eq!(moved, 1);
            let issues = validate(&b, &tol);
            assert_eq!(issues.is_empty(), valid, "{lift}: {issues:?}");
        }
    }

    /// Edges meeting the axis of a surface of revolution between their
    /// ends are found, wherever along them and however they cross it;
    /// edges meeting it at an end, or passing it at a distance, are not.
    #[test]
    fn edges_through_the_axis_of_a_surface_of_revolution_are_found() {
        let tol = 1e-7;
        let fr = Frame::WORLD;
        let line = |a: DVec3, b: DVec3| Curve3::Line {
            origin: a,
            dir: b - a,
        };
        let circle = |center: DVec3, x: DVec3, y: DVec3, radius: f64| Curve3::Circle {
            center,
            x,
            y,
            radius,
        };
        let through = [
            // Diameters of a flat cap, the axis at either end's side.
            (
                line(DVec3::new(1.0, 0.0, 2.0), DVec3::new(-0.04, 0.0, 2.0)),
                0.0,
                1.0,
            ),
            (
                line(DVec3::new(0.03, 0.0, 0.0), DVec3::new(-1.0, 0.0, 0.0)),
                0.0,
                1.0,
            ),
            (
                line(DVec3::new(0.3, 0.4, 0.0), DVec3::new(-0.3, -0.4, 0.0)),
                0.0,
                1.0,
            ),
            // A meridian of a sphere over a pole.
            (circle(DVec3::ZERO, DVec3::X, DVec3::Z, 1.0), 1.0, 2.0),
            // A circle in a plane across the axis.
            (
                circle(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 0.5),
                2.0,
                4.0,
            ),
        ];
        for (c, t0, t1) in &through {
            assert!(through_axis(c, *t0, *t1, &fr, tol).is_some(), "{c:?}");
        }
        let clear = [
            // A radius ending at the axis, and one leaving it.
            (
                line(DVec3::new(1.0, 0.0, 2.0), DVec3::new(0.0, 0.0, 2.0)),
                0.0,
                1.0,
            ),
            (
                line(DVec3::new(0.0, 0.0, 2.0), DVec3::new(0.0, 1.0, 2.0)),
                0.0,
                1.0,
            ),
            // A chord passing the axis a micrometre away.
            (
                line(DVec3::new(1.0, 1e-6, 0.0), DVec3::new(-1.0, 1e-6, 0.0)),
                0.0,
                1.0,
            ),
            // A circle about the axis, and a meridian ending at a pole.
            (circle(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0), 0.0, 3.0),
            (
                circle(DVec3::ZERO, DVec3::X, DVec3::Z, 1.0),
                0.0,
                0.5 * std::f64::consts::PI,
            ),
        ];
        for (c, t0, t1) in &clear {
            assert!(through_axis(c, *t0, *t1, &fr, tol).is_none(), "{c:?}");
        }
    }

    /// The tests of NURBS edges' spans and of circles against the axis
    /// clear only edges that sampling clears too: arcs of circles and
    /// ellipses written as rational quadratics, and circles, in planes at
    /// any angle to the axis, through it, near it and away from it, whole
    /// and in part; and they clear most of those away from it, and every
    /// circle about the axis.
    #[test]
    fn the_quick_tests_of_the_axis_clear_only_clear_edges() {
        use crate::boolean::testing::Rng;
        use std::f64::consts::{FRAC_1_SQRT_2, TAU};
        let tol = 1e-7;
        let fr = Frame::WORLD;
        let mut rng = Rng(0xA815_0001);
        let mut counts = [[0usize; 3]; 2];
        let dir = |rng: &mut Rng| {
            DVec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            )
            .normalize_or(DVec3::X)
        };
        for case in 0..6000 {
            let center = DVec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ) * 10f64.powf(rng.range(-6.0, 0.0));
            let x = dir(&mut rng);
            let y = x
                .cross(dir(&mut rng))
                .normalize_or(x.any_orthonormal_vector());
            let (ra, rb) = (rng.range(0.1, 2.0), rng.range(0.1, 2.0));
            // A quarter of an ellipse (from centre + a to centre + b, the
            // middle control point at centre + a + b with weight cos 45°),
            // or a circle.
            let nurbs = case % 2 == 0;
            let make = |center: DVec3| {
                if nurbs {
                    Curve3::Nurbs(
                        NurbsCurve::new(
                            2,
                            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                            vec![center + x * ra, center + x * ra + y * rb, center + y * rb],
                            Some(vec![1.0, FRAC_1_SQRT_2, 1.0]),
                        )
                        .expect("an arc"),
                    )
                } else {
                    Curve3::Circle {
                        center,
                        x,
                        y,
                        radius: ra,
                    }
                }
            };
            let whole = if nurbs { 1.0 } else { TAU };
            let (t0, t1) = if case % 4 < 2 {
                (0.0, whole)
            } else {
                let t0 = rng.range(0.0, 0.8);
                (t0 * whole, rng.range(t0 + 0.1, 1.0) * whole)
            };
            // A third through the axis between the ends, a third passing
            // it 10 to 1,000 tolerances away.
            let curve = match case % 3 {
                0 => make(center),
                k => {
                    let at = make(center).point(t0 + (t1 - t0) * rng.range(0.1, 0.9));
                    let gap = if k == 1 {
                        DVec3::ZERO
                    } else {
                        dir(&mut rng).cross(DVec3::Z).normalize_or(DVec3::X)
                            * (tol * 10f64.powf(rng.range(1.0, 3.0)))
                    };
                    make(center - DVec3::new(at.x, at.y, 0.0) + gap)
                }
            };
            let quick = match &curve {
                Curve3::Nurbs(n) => hull_clear_of_axis(n, t0, t1, &fr, tol),
                &Curve3::Circle {
                    center,
                    x,
                    y,
                    radius,
                } => circle_clear_of_axis(center, x, y, radius, &fr, tol),
                Curve3::Line { .. } => unreachable!(),
            };
            let found = sampled_through_axis(&curve, t0, t1, &fr, tol);
            let k = usize::from(!nurbs);
            counts[k][2] += usize::from(found.is_some());
            if quick {
                counts[k][0] += 1;
                assert!(found.is_none(), "case {case}: {curve:?} {t0} {t1}");
                assert_eq!(through_axis(&curve, t0, t1, &fr, tol), None);
            } else {
                counts[k][1] += 1;
                assert_eq!(through_axis(&curve, t0, t1, &fr, tol), found);
            }
        }
        // Per kind: cleared quickly, sampled, found through the axis.
        for [quick, sampled, through] in counts {
            assert!(quick > 300 && sampled > 700 && through > 400, "{counts:?}");
        }
        // Circles about the axis, square to it.
        for k in 0..100 {
            let r = 10f64.powf(-6.0 + 0.07 * k as f64);
            let tilt = DVec3::new(1e-13, -2e-13, 1.0).normalize();
            let x = tilt.any_orthonormal_vector();
            assert!(circle_clear_of_axis(
                DVec3::new(0.0, 0.0, r),
                x,
                tilt.cross(x),
                r,
                &fr,
                tol
            ));
        }
    }

    /// Circles in planes along a leaning axis, about a point of it (a
    /// ball's meridians), and in planes tilted from those so that they pass
    /// within ten times the tolerance either side of it: the quick test
    /// clears none that sampling finds through the axis, and none of the
    /// meridians (it cleared some from a kilometre down, whose least
    /// singular value lost `√ε` to cancellation; found by review).
    #[test]
    fn circles_along_a_leaning_axis_are_not_cleared_through_it() {
        use crate::boolean::testing::Rng;
        use std::f64::consts::TAU;
        let tol = 1e-5;
        // The case review found: a meridian of a ball of radius 700.
        let z = DVec3::new(0.6132052342058546, -0.05385148810082932, -0.788085882357925);
        let x = DVec3::new(0.3917140890063436, 0.7540871621540453, -0.5271741878624077);
        let fr = Frame {
            origin: DVec3::ZERO,
            x,
            y: z.cross(x),
            z,
        };
        let (x, y) = (
            x,
            DVec3::new(-0.47214594291149464, 0.6565146549758518, 0.5882743546969849),
        );
        assert!(!circle_clear_of_axis(DVec3::ZERO, x, y, 700.0, &fr, tol));
        let meridian = Curve3::Circle {
            center: DVec3::ZERO,
            x,
            y,
            radius: 700.0,
        };
        assert!(through_axis(&meridian, 0.0, TAU, &fr, tol).is_some());
        let mut rng = Rng(0xA815_0003);
        let unit = |rng: &mut Rng| {
            DVec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            )
            .normalize_or(DVec3::X)
        };
        let (mut cleared, mut through) = (0, 0);
        for case in 0..3000 {
            let z = unit(&mut rng);
            let fx = z.any_orthonormal_vector();
            let fr = Frame {
                origin: unit(&mut rng) * rng.range(0.0, 100.0),
                x: fx,
                y: z.cross(fx),
                z,
            };
            let center = fr.origin + z * rng.range(-100.0, 100.0);
            let r = if case % 2 == 0 {
                rng.range(500.0, 1000.0)
            } else {
                10f64.powf(rng.range(0.0, 3.0))
            };
            // The plane's normal, square to the axis for a meridian, else
            // tilted toward it so that the circle passes the axis at
            // r sin θ.
            let side = z.cross(unit(&mut rng)).normalize_or(fx);
            let theta = if case % 3 == 0 {
                0.0
            } else {
                (tol * 10f64.powf(rng.range(-1.0, 1.0)) / r).asin()
            };
            let n = side * theta.cos() + z * theta.sin();
            // The start a third to a sixth of a turn from the poles.
            let phi = rng.range(0.3, 1.2);
            let w = n.cross(z).normalize();
            let x = (z * phi.cos() + w * phi.sin() - n * n.dot(w * phi.sin())).normalize();
            let y = n.cross(x).normalize();
            let quick = circle_clear_of_axis(center, x, y, r, &fr, tol);
            let circle = Curve3::Circle {
                center,
                x,
                y,
                radius: r,
            };
            let sampled = sampled_through_axis(&circle, 0.0, TAU, &fr, tol);
            if case % 3 == 0 {
                assert!(!quick, "case {case}: a meridian of radius {r} cleared");
                assert!(sampled.is_some(), "case {case}: {circle:?}");
            }
            if quick {
                cleared += 1;
                assert!(sampled.is_none(), "case {case}: {circle:?} cleared");
            }
            through += usize::from(sampled.is_some());
        }
        assert!(cleared > 300 && through > 1500, "{cleared} {through}");
    }

    /// The closed form for lines finds what sampling finds: lines through
    /// the axis between their ends and lines passing it 10 to 1,000
    /// tolerances away, lines from a point on it (a cone's ruling from its
    /// apex), lines that would meet it beyond their ends, and lines
    /// parallel to it, off it; about axes at any angle and
    /// place, either way along the edge. Along the axis, where sampling
    /// decides by its rounding, the closed form finds the edge's middle.
    #[test]
    fn lines_meet_the_axis_where_sampling_finds_them() {
        use crate::boolean::testing::Rng;
        let tol = 1e-7;
        let mut rng = Rng(0xA815_0002);
        let unit = |rng: &mut Rng| {
            DVec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            )
            .normalize_or(DVec3::X)
        };
        let mut through = 0;
        for case in 0..6000 {
            let z = unit(&mut rng);
            let x = z.any_orthonormal_vector();
            let fr = Frame {
                origin: unit(&mut rng) * 10f64.powf(rng.range(-3.0, 1.0)),
                x,
                y: z.cross(x),
                z,
            };
            let side = z.cross(unit(&mut rng)).normalize_or(x);
            let len = 10f64.powf(rng.range(-1.0, 0.5));
            let on_axis = fr.origin + z * rng.range(-2.0, 2.0);
            let (origin, dir) = match case % 6 {
                // The gap along the common perpendicular of the line and
                // the axis, so that it is the line's least distance.
                0 | 1 => {
                    let dir = unit(&mut rng) * len;
                    let gap = if case % 6 == 0 {
                        DVec3::ZERO
                    } else {
                        z.cross(dir).normalize_or(side) * (tol * 10f64.powf(rng.range(1.0, 3.0)))
                    };
                    (on_axis + gap - dir * rng.range(0.1, 0.9), dir)
                }
                2 => (on_axis, unit(&mut rng) * len),
                3 => (
                    on_axis + side * (tol * 10f64.powf(rng.range(1.0, 5.0))),
                    z * len,
                ),
                4 => (on_axis, z * len),
                _ => {
                    let dir = unit(&mut rng) * len;
                    let at = rng.range(0.2, 2.0);
                    (
                        on_axis - dir * if case % 4 < 2 { -at } else { 1.0 + at },
                        dir,
                    )
                }
            };
            let (t0, t1) = if case % 2 == 0 {
                (0.0, 1.0)
            } else {
                (1.0, 0.0)
            };
            let line = Curve3::Line { origin, dir };
            let found = line_through_axis(origin, dir, t0, t1, &fr, tol);
            assert_eq!(through_axis(&line, t0, t1, &fr, tol), found);
            if case % 6 == 4 {
                let middle = origin + dir * 0.5;
                assert!(found.is_some_and(|p| (p - middle).length() <= 1e-12 * len));
                continue;
            }
            let sampled = sampled_through_axis(&line, t0, t1, &fr, tol);
            assert_eq!(found.is_some(), case % 6 == 0, "case {case}: {line:?}");
            assert_eq!(found.is_some(), sampled.is_some(), "case {case}: {line:?}");
            if let (Some(p), Some(q)) = (found, sampled) {
                assert!((p - q).length() <= 2.0 * tol, "case {case}: {p} {q}");
                through += 1;
            }
        }
        assert_eq!(through, 1000);
        // Rulings from a point of the axis, leaning half a tolerance off it
        // at their far ends: along it within the tolerance, the middle is
        // reported when farther than the tolerance from both ends.
        let fr = Frame::WORLD;
        for (len, along) in [(5.0, true), (1.5, false)] {
            let dir = DVec3::new(0.5 * tol, 0.0, len * tol);
            let found = line_through_axis(DVec3::ZERO, dir, 0.0, 1.0, &fr, tol);
            assert_eq!(found, along.then_some(dir * 0.5), "{len}");
        }
    }

    /// A face of a revolve turned half a turn on its surface lies on that
    /// surface still, outside its frame's half turn: refused (ADR 0008).
    #[test]
    fn faces_of_revolution_keep_to_their_half_turn() {
        use crate::profile::Profile;
        use crate::sweep::{self, RevolveSpec};
        use glam::DVec2;
        let tol = Tolerance::default();
        let pipe = Profile::polygon(
            Frame::WORLD,
            &[
                DVec2::new(1.0, 0.0),
                DVec2::new(2.0, 0.0),
                DVec2::new(2.0, 1.0),
                DVec2::new(1.0, 1.0),
            ],
        );
        let spec = RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: std::f64::consts::TAU,
        };
        let mut b = sweep::revolve(&pipe, &spec, &tol).expect("a pipe").0;
        assert_eq!(validate(&b, &tol), Vec::<String>::new());
        let f = b
            .faces
            .iter_mut()
            .find(|f| matches!(f.surface, Surface::Revolution { .. }))
            .expect("a face of revolution");
        if let Surface::Revolution { frame, .. } = &mut f.surface {
            frame.x = -frame.x;
            frame.y = -frame.y;
        }
        let issues = validate(&b, &tol);
        assert!(issues.iter().any(|m| m.contains("half turn")), "{issues:?}");
    }
}
