//! B-rep validation, topology statistics and mass properties.

use glam::DVec3;
use rustc_hash::FxHashMap;

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
/// The distance from the axis is sampled along the edge and each sampled
/// local least refined by golden-section search between its neighbours.
fn through_axis(c: &Curve3, t0: f64, t1: f64, fr: &Frame, tol: f64) -> Option<DVec3> {
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

/// Checks topological and geometric consistency. Returns a list of problems;
/// an empty list means the body is valid.
pub(crate) fn validate(b: &Brep, tol: &Tolerance) -> Vec<String> {
    let mut issues = Vec::new();
    // Geometry tolerance for "lies on": a generous multiple of the model
    // resolution to absorb projection round-off in iterative inversions.
    let on_tol = tol.linear * 100.0;
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
