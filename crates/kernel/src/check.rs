//! B-rep validation, topology statistics and mass properties.

use glam::DVec3;
use rustc_hash::FxHashMap;

use crate::geom::Surface;
use crate::tess::{self, TessellationOptions};
use crate::topo::Brep;
use crate::{KernelError, Tolerance, TopologyCounts};

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
        match shell_volumes(b, tol) {
            Ok(volumes) => {
                for (si, solid) in b.solids.iter().enumerate() {
                    for (k, sh) in solid.shells.iter().enumerate() {
                        let v = volumes[sh.idx()].0;
                        if k == 0 && v <= 0.0 {
                            issues.push(format!("solid {si}: outer shell is inside-out"));
                        }
                        if k > 0 && v >= 0.0 {
                            issues.push(format!("solid {si}: void shell {k} is not inverted"));
                        }
                    }
                }
            }
            Err(e) => issues.push(format!("cannot tessellate: {e}")),
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
/// body's size from its surfaces, with a bound on the volume's error: twice
/// the mesh's area times that distance.
pub(crate) fn shell_volumes(b: &Brep, tol: &Tolerance) -> Result<Vec<(f64, f64)>, KernelError> {
    let opts = TessellationOptions {
        chord_tolerance: bounding_size(b).max(1e-9) * 1e-2,
        max_angle: 0.6,
    };
    let t = tess::tessellate(b, &opts, tol)?;
    Ok(b.shells
        .iter()
        .map(|shell| {
            let (mut v, mut area) = (0.0, 0.0);
            for f in &shell.faces {
                let ft = &t.faces[f.idx()];
                for tri in &ft.triangles {
                    let [a, bb, c] = tri.map(|i| ft.positions[i as usize]);
                    v += a.dot(bb.cross(c)) / 6.0;
                    area += 0.5 * (bb - a).cross(c - a).length();
                }
            }
            (v, 2.0 * area * opts.chord_tolerance)
        })
        .collect())
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
}
