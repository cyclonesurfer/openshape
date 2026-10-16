//! An oracle for results that touch themselves (ADR 0024, design note
//! 0012 §18.5), shared by the booleans' tests: a point of an edge of a
//! result within half the tolerance of a face that does not bound the
//! result there. It reads only the result's mesh and adjacency, so it is
//! independent of the booleans' own decision of contacts.
//!
//! A pinch the booleans miss leaves an edge or a vertex of one operand on
//! a face of the other in the result: along an edge, the edge lies on the
//! face; at a point, the vertex does. Every such edge or vertex is an
//! entity of the result, so sampling the result's edges finds it. Faces
//! touching faces with no edge between them (surfaces that touch, slice 2)
//! are not found by this oracle.

#![allow(dead_code)]

use opendatum_kernel::{
    Body, BodyMesh, BrepKernel, CurveKind, DVec3, FaceIx, Kernel, TessellationOptions,
};

/// The distance from `p` to the triangle `a b c`, by the closest point's
/// Voronoi region of the triangle (Ericson, *Real-Time Collision
/// Detection*, 2005, §5.1.5).
pub fn triangle_distance(p: DVec3, a: DVec3, b: DVec3, c: DVec3) -> f64 {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.length();
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return bp.length();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return (p - (a + ab * (d1 / (d1 - d3)))).length();
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return cp.length();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return (p - (a + ac * (d2 / (d2 - d6)))).length();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).length();
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    (p - (a + ab * v + ac * w)).length()
}

/// Where `r` touches itself, if it does: the first point (by edge, then
/// along it) of an edge's mesh polyline, or the middle of one of its
/// segments, within half the tolerance of a face that is not along the
/// edge there (a face the edge does not bound and, for a point within the
/// tolerance of one of the edge's ends, that is not at that end: every
/// point that near a vertex is that near each face at it), meshed by
/// `opts`.
pub fn touches_itself(
    k: &BrepKernel,
    r: &Body,
    opts: &TessellationOptions,
) -> Result<Option<String>, String> {
    let mesh = k.tessellate(r, opts).map_err(|e| e.to_string())?;
    Ok(touches_itself_in(k, r, &mesh, opts.chord_tolerance))
}

/// [`touches_itself`] on `r`'s mesh, whose chords are within `chord` of
/// the curves and surfaces. A flat face's mesh is first moved onto its
/// plane: its corners are the result's vertices, which may be off the
/// plane by up to the tolerance, so the mesh as made can lean towards a
/// point the face itself is a tolerance from. Moved, the mesh of a flat
/// face bounded by straight edges is the face, and decides alone. For any
/// other face, a
/// point within `chord` and half the tolerance of its mesh is a candidate,
/// confirmed by its distance to the face's surface; so a point on the
/// surface just outside such a face, within `chord` of its boundary, would
/// be reported too (none has been seen).
pub fn touches_itself_in(k: &BrepKernel, r: &Body, mesh: &BodyMesh, chord: f64) -> Option<String> {
    let near = 0.5 * k.tolerance().linear;
    let adj = k.adjacency(r);
    let kinds = k.edge_kinds(r);
    let mut ends: Vec<Vec<usize>> = vec![Vec::new(); adj.edge_faces.len()];
    for (v, edges) in adj.vertex_edges.iter().enumerate() {
        for e in edges {
            ends[e.0 as usize].push(v);
        }
    }
    let at_vertex: Vec<Vec<u32>> = (adj.vertex_edges.iter())
        .map(|edges| {
            let mut f: Vec<u32> = (edges.iter())
                .flat_map(|e| adj.edge_faces[e.0 as usize].iter().map(|f| f.0))
                .collect();
            f.sort_unstable();
            f.dedup();
            f
        })
        .collect();
    // Per face: whether its mesh is the face, and how far from it a point
    // of the face can be.
    let planes: Vec<_> = (0..mesh.faces.len())
        .map(|f| k.face_plane(r, FaceIx(f as u32)))
        .collect();
    let positions: Vec<Vec<DVec3>> = (mesh.faces.iter().zip(&planes))
        .map(|(face, plane)| {
            (face.positions.iter())
                .map(|&p| match plane {
                    Some(pl) => p - pl.normal * (p - pl.point).dot(pl.normal),
                    None => p,
                })
                .collect()
        })
        .collect();
    let mut exact: Vec<bool> = planes.iter().map(Option::is_some).collect();
    for (e, faces) in adj.edge_faces.iter().enumerate() {
        if kinds[e] != CurveKind::Line {
            for f in faces {
                exact[f.0 as usize] = false;
            }
        }
    }
    let pad = |f: usize| if exact[f] { near } else { chord + near };
    let bounds: Vec<(DVec3, DVec3)> = (positions.iter().enumerate())
        .map(|(i, ps)| {
            let lo = ps.iter().fold(DVec3::INFINITY, |m, p| m.min(*p));
            let hi = ps.iter().fold(DVec3::NEG_INFINITY, |m, p| m.max(*p));
            (lo - DVec3::splat(pad(i)), hi + DVec3::splat(pad(i)))
        })
        .collect();
    let report = |e: usize, f: usize, d: f64, p: DVec3| {
        Some(format!(
            "edge {e} comes within {d:.3e} m of face {f} at ({:.6}, {:.6}, {:.6}) mm",
            p.x * 1e3,
            p.y * 1e3,
            p.z * 1e3
        ))
    };
    let tol = k.tolerance().linear;
    for (e, line) in mesh.edges.iter().enumerate() {
        let along: Vec<u32> = adj.edge_faces[e].iter().map(|f| f.0).collect();
        // The faces a point of the edge may be near: those along it, and
        // those at an end within the tolerance of the point, which every
        // point that near the end is near.
        let allowed = |p: DVec3| -> Vec<u32> {
            let mut a = along.clone();
            for &v in &ends[e] {
                if (p - mesh.vertices[v]).length() <= tol {
                    a.extend(&at_vertex[v]);
                }
            }
            a
        };
        let mut samples: Vec<DVec3> = Vec::new();
        for (i, p) in line.iter().enumerate() {
            samples.push(*p);
            if i + 1 < line.len() {
                samples.push(0.5 * (*p + line[i + 1]));
            }
        }
        for p in samples {
            let allowed = allowed(p);
            for (f, face) in mesh.faces.iter().enumerate() {
                let (lo, hi) = bounds[f];
                if allowed.contains(&(f as u32)) || p.cmplt(lo).any() || p.cmpgt(hi).any() {
                    continue;
                }
                let Some(d) = (face.triangles.iter())
                    .map(|t| {
                        let [a, b, c] = t.map(|i| positions[f][i as usize]);
                        triangle_distance(p, a, b, c)
                    })
                    .find(|&d| d <= pad(f))
                else {
                    continue;
                };
                if exact[f] {
                    return report(e, f, d, p);
                }
                if let Some(t) = k.tangent_plane(r, FaceIx(f as u32), p) {
                    let d = (p - t.point).length();
                    if d <= near {
                        return report(e, f, d, p);
                    }
                }
            }
        }
    }
    None
}
