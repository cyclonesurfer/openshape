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

/// The closest points of the segments `p1 q1` and `p2 q2`: the point on the
/// first and its distance from the second (Ericson, *Real-Time Collision
/// Detection*, 2005, §5.1.9).
pub fn segment_distance(p1: DVec3, q1: DVec3, p2: DVec3, q2: DVec3) -> (DVec3, f64) {
    let (d1, d2, r) = (q1 - p1, q2 - p2, p1 - p2);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let tiny = f64::MIN_POSITIVE;
    let (s, t) = if a <= tiny && e <= tiny {
        (0.0, 0.0)
    } else if a <= tiny {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(r);
        if e <= tiny {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(d2);
            let denom = a * e - b * b;
            let s = if denom > 0.0 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + f) / e;
            if t < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t > 1.0 {
                (((b - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s, t)
            }
        }
    };
    let x = p1 + d1 * s;
    (x, (x - (p2 + d2 * t)).length())
}

/// The point of the segment `p q` nearest the triangle `a b c`, and its
/// distance. Where they do not meet, the nearest points are an end of the
/// segment and the triangle, or the segment and a side of the triangle;
/// where the segment crosses the triangle's plane, the crossing is tried
/// too, and is at no distance if it is inside the triangle.
pub fn segment_triangle_distance((p, q): (DVec3, DVec3), [a, b, c]: [DVec3; 3]) -> (DVec3, f64) {
    let mut best = (p, triangle_distance(p, a, b, c));
    let mut take = |x: (DVec3, f64)| {
        if x.1 < best.1 {
            best = x;
        }
    };
    take((q, triangle_distance(q, a, b, c)));
    for (u, v) in [(a, b), (b, c), (c, a)] {
        take(segment_distance(p, q, u, v));
    }
    let n = (b - a).cross(c - a);
    let (dp, dq) = ((p - a).dot(n), (q - a).dot(n));
    if (dp <= 0.0) != (dq <= 0.0) || dp == 0.0 {
        let x = if dp == dq {
            p
        } else {
            p + (q - p) * (dp / (dp - dq))
        };
        take((x, triangle_distance(x, a, b, c)));
    }
    best
}

/// Where `r` touches itself, if it does: the first point (by edge, then
/// along it) of an edge's mesh polyline, or the middle of one of its
/// segments, within half the tolerance of a face that is not along the
/// edge there (a face the edge does not bound and, for a point within the
/// tolerance of one of the edge's ends, that is not at that end: every
/// point that near a vertex is that near each face at it; nor a flat face
/// that a straight edge leaves at a corner, [`touches_itself_in`]), meshed
/// by `opts`.
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
/// face bounded by straight edges is the face, and decides alone.
///
/// A straight edge is measured against such a face whole, not at sample
/// points, where the face has no corner at the edge's ends: its nearest
/// point to the face's triangles, wherever along it that is. (The booleans'
/// own check of their results, `check_apart`, tests the ends and the middle
/// of each edge; this keeps the oracle from sharing that blind spot.) Where
/// the face has a corner at an end, the edge leaves that corner and comes
/// within half the tolerance of the face for half the tolerance over the
/// sine of the angle between them, which no fixed distance from the corner
/// excuses. Where the edge's other end is more than the tolerance from the
/// face's plane, the edge leaves the face at that corner (crossing its
/// plane at a slant, if at all, where the tolerance puts the crossing at
/// the corner) and the face is let off: a straight line is near a plane
/// only over one stretch, and this one is the corner's. The booleans' own
/// check lets these corners off for the same reason, so the oracle cannot
/// find a fault there; it reported three results of the contact fuzz at a
/// hundred times its default run that are of this kind (contacts 1335,
/// 2714 and 5358, `sharp_corners_the_check_of_places_lets_off`). Where the
/// other end is within the tolerance of the plane, the edge lies along the
/// face's plane, and the face is tested at the sample points, as curved
/// edges and faces are.
///
/// For any other face, a point within `chord` and half the tolerance of
/// its mesh is a candidate, confirmed by its distance to the face's
/// surface; so a point on the surface just outside such a face, within
/// `chord` of its boundary, would be reported too (none has been seen).
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
        // Flat faces at one end of a straight edge whose other end is more
        // than the tolerance from their planes: the edge leaves them at
        // that corner, and is near them only about it.
        let mut leaving: Vec<u32> = Vec::new();
        if let (CurveKind::Line, &[v, w]) = (kinds[e], ends[e].as_slice()) {
            for (v, w) in [(v, w), (w, v)] {
                for &f in &at_vertex[v] {
                    let off = (planes[f as usize].as_ref())
                        .is_some_and(|pl| (mesh.vertices[w] - pl.point).dot(pl.normal).abs() > tol);
                    if exact[f as usize] && off {
                        leaving.push(f);
                    }
                }
            }
        }
        // The faces a point of the edge may be near: those along it, those
        // it leaves at a corner, and those at an end within the tolerance
        // of the point, which every point that near the end is near.
        let allowed = |p: DVec3| -> Vec<u32> {
            let mut a = along.clone();
            a.extend(&leaving);
            for &v in &ends[e] {
                if (p - mesh.vertices[v]).length() <= tol {
                    a.extend(&at_vertex[v]);
                }
            }
            a
        };
        // A straight edge against the flat faces with no corner at its
        // ends, whole.
        if kinds[e] == CurveKind::Line {
            let (p, q) = (line[0], line[line.len() - 1]);
            let (lo, hi) = (p.min(q), p.max(q));
            for (f, face) in mesh.faces.iter().enumerate() {
                let (flo, fhi) = bounds[f];
                if !exact[f]
                    || along.contains(&(f as u32))
                    || ends[e].iter().any(|&v| at_vertex[v].contains(&(f as u32)))
                    || hi.cmplt(flo).any()
                    || lo.cmpgt(fhi).any()
                {
                    continue;
                }
                let nearest = (face.triangles.iter())
                    .map(|t| segment_triangle_distance((p, q), t.map(|i| positions[f][i as usize])))
                    .find(|x| x.1 <= near);
                if let Some((x, d)) = nearest {
                    return report(e, f, d, x);
                }
            }
        }
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
