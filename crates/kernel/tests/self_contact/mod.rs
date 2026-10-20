//! An oracle for results that touch themselves (ADR 0024, design note
//! 0012 §18.5), shared by the booleans' tests: a point of an edge of a
//! result within half the tolerance of a face that does not bound the
//! result there. It reads only the result (its mesh, adjacency, edge
//! kinds, and faces' planes and tangent planes), not the operands or the
//! booleans' records, so it is independent of the booleans' own decision
//! of contacts.
//!
//! A pinch the booleans miss leaves an edge or a vertex of one operand on
//! a face of the other in the result: along an edge, the edge lies on the
//! face; at a point, the vertex does. Every such edge or vertex is an
//! entity of the result, so sampling the result's edges finds it. Faces
//! touching faces with no edge between them (surfaces that touch, slice 2)
//! are not found by this oracle.
//!
//! Nor are two faces that cross each other away from their edges: a flap
//! or a fold of the result under a tolerance thick, whose faces pass
//! through each other along a line no edge of the result follows. A
//! second check, [`faces_cross`], finds those between flat faces: the
//! common line of two faces' planes passing through both faces more than
//! the tolerance from each one's boundary.

#![allow(dead_code)]

use std::fmt;

use opendatum_kernel::{
    Body, BodyMesh, BrepKernel, CurveKind, DVec2, DVec3, FaceIx, Kernel, KernelError,
    TessellationOptions,
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

/// The distance from `p` to the segment `a b`.
pub fn point_segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let d = b - a;
    let l = d.length_squared();
    let t = if l > 0.0 {
        ((p - a).dot(d) / l).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + d * t)).length()
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

/// Where `r` touches itself, if it does: a point of an edge within half
/// the tolerance of a face that is not along the edge there (a face the
/// edge does not bound and, for a point within the tolerance of one of the
/// edge's ends, that is not at that end: every point that near a vertex is
/// that near each face at it; nor a flat face that a straight edge leaves
/// at a corner, or whose edge at a corner it leaves, [`touches_itself_in`]),
/// meshed by `opts`. Edge by edge: for a straight edge, against each flat
/// face bounded by straight edges with no corner at the edge's ends, the
/// edge's point nearest the first of the face's triangles within that
/// distance (the edge measured whole); then the first point of the edge's
/// mesh polyline, or the middle of one of its segments, within that
/// distance of a face.
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
/// point to the first of the face's triangles within half the tolerance,
/// wherever along it that is. (The booleans' own check of their results,
/// `check_apart`, tests the ends and the middle of each edge; this keeps
/// the oracle from sharing that blind spot.) Where the face has a corner at
/// an end, the edge leaves that corner and comes within half the tolerance
/// of the face for half the tolerance over the sine of the angle between
/// them, which no fixed distance from the corner excuses. Where the edge's
/// other end is more than the tolerance from the face's plane, the edge
/// leaves the face at that corner (crossing its plane at a slant, if at
/// all, where the tolerance puts the crossing at the corner) and the face
/// is let off: a straight line is near a plane only over one stretch, and
/// this one is the corner's. The booleans' own check lets these corners off
/// for the same reason, so the oracle cannot find a fault there; it
/// reported three results of the contact fuzz at a hundred times its
/// default run that are of this kind (contacts 1335, 2714 and 5358,
/// `sharp_corners_the_check_of_places_lets_off`). Where the other end is
/// within the tolerance of the plane, the edge lies along the face's plane,
/// and the face is tested at the sample points, as curved edges and faces
/// are.
///
/// There a face is let off at a point of the edge where the edge leaves
/// one of the face's edges at their shared corner: where the point's
/// nearest point on the face is on a straight edge of the face that ends
/// at a vertex W the two edges share (its distance from that edge within
/// the tolerance over 1,024 of its distance from the face, for rounding),
/// and the edge's other end is more than the tolerance from that edge's
/// line. Two straight edges leaving W at an angle θ are within half the
/// tolerance of each other for half the tolerance over sin θ from W and no
/// farther, and a straight line is near another only over one stretch,
/// here the corner's: the closeness is the corner's, not a contact, as for
/// a face the edge leaves at a corner above. It is found where a meeting
/// splits a face at the acute tip of a notch: the face's pieces either side
/// of the notch then share only the tip, and each notch edge is within half
/// the tolerance of the piece across the notch for half the tolerance over
/// the sine of the notch's angle from the tip, as near as it is in the
/// operand to the notch's other wall, which the rule above lets off there.
/// The 5,733 reports of the notch-tip probes (design note 0012 §21, sweeps
/// 6, 8, 10 and 11) that this lets off are all of this kind: the shared vertex
/// is the notch's tip, the edge bounds a face in the face's plane, and the
/// point is nearer the tip than half the tolerance over the sine of the
/// angle. A point whose nearest point is inside the face (an edge lying
/// over a face it does not bound), or on an edge of the face at no vertex
/// the two share, or on an edge the other runs along (its far end within
/// the tolerance of that edge's line: a slit), is still reported.
///
/// For any other face, a point within `chord` and half the tolerance of its
/// mesh is a candidate, confirmed by its distance to the face's surface; so
/// a point on the surface just outside such a face, within `chord` of its
/// boundary, would be reported too. Its reports on the curved fuzzes are
/// not yet sorted into contacts and such points, so those fuzzes do not run
/// it yet (design note 0012 §19).
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
        // Whether the edge leaves an edge of exact face `f` at a corner they
        // share, and the point `p` of it is nearest that edge of the face.
        let corner = |e: usize, f: usize, p: DVec3| -> bool {
            let (CurveKind::Line, &[v0, v1], Some(pl)) = (kinds[e], ends[e].as_slice(), &planes[f])
            else {
                return false;
            };
            let onto = |q: DVec3| q - pl.normal * (q - pl.point).dot(pl.normal);
            let from_face = (mesh.faces[f].triangles.iter())
                .map(|t| {
                    let [a, b, c] = t.map(|i| positions[f][i as usize]);
                    triangle_distance(p, a, b, c)
                })
                .fold(f64::INFINITY, f64::min);
            for (w, far) in [(v0, v1), (v1, v0)] {
                if !at_vertex[w].contains(&(f as u32)) {
                    continue;
                }
                for g in &adj.vertex_edges[w] {
                    let g = g.0 as usize;
                    let Some(&u) = ends[g].iter().find(|&&u| u != w) else {
                        continue;
                    };
                    if g == e
                        || kinds[g] != CurveKind::Line
                        || !adj.edge_faces[g].iter().any(|x| x.0 as usize == f)
                    {
                        continue;
                    }
                    let (a, b) = (mesh.vertices[w], mesh.vertices[u]);
                    let line = (b - a).normalize_or_zero();
                    if line == DVec3::ZERO
                        || (mesh.vertices[far] - a)
                            .reject_from_normalized(line)
                            .length()
                            <= tol
                    {
                        continue;
                    }
                    if point_segment_distance(p, onto(a), onto(b)) <= from_face + tol / 1024.0 {
                        return true;
                    }
                }
            }
            false
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
                    if corner(e, f, p) {
                        continue;
                    }
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

/// The sine of the angle under which [`faces_cross`] takes two planes for
/// parallel (the model's angular tolerance): where faces at a smaller
/// angle reach their planes' common line, they are within the tolerance
/// of each other over a band hundreds of metres wide about it, a contact
/// rather than a crossing.
pub const PARALLEL: f64 = 1e-9;

/// A flat face as [`pair_crossing`] reads it: its plane and the segments
/// of its boundary.
#[derive(Clone, Debug)]
pub struct FlatFace {
    /// A point of the face's plane.
    pub point: DVec3,
    /// The plane's unit normal.
    pub normal: DVec3,
    /// The face's boundary: each straight edge as the segment between its
    /// vertices, each curved edge as the segments of its polyline, as many
    /// times as the edge bounds the face. Each loop must close exactly
    /// (the same point ending one segment and starting the next).
    pub segments: Vec<(DVec3, DVec3)>,
}

/// Where the common line of two flat faces' planes passes through both
/// faces deepest ([`pair_crossing`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PairCrossing {
    /// The point of the common part farthest from the nearer of the two
    /// faces' boundaries.
    pub point: DVec3,
    /// Its distance from each face's boundary, in the face's plane.
    pub depths: [f64; 2],
    /// The smaller of the two: how deep the faces cross.
    pub depth: f64,
    /// The angle between the two planes, in radians (0 to π/2).
    pub angle: f64,
    /// Whether the faces' outward normals point the same way (a positive
    /// dot product).
    pub same_way: bool,
}

/// The deepest crossing of a body's flat faces ([`faces_cross`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Crossing {
    /// The two faces, in face order.
    pub faces: [usize; 2],
    /// Where and how deep they cross.
    pub at: PairCrossing,
    /// The distance from the point to the body's nearest vertex.
    pub vertex: f64,
    /// Whether the two faces share an edge, a vertex only, or nothing.
    pub shared: Shared,
    /// How many pairs of faces cross.
    pub pairs: usize,
    /// The tolerance the depths are measured against.
    pub tolerance: f64,
}

/// What two faces share ([`Crossing::shared`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shared {
    /// An edge bounds both.
    Edge,
    /// A vertex but no edge.
    Vertex,
    /// Neither.
    Nothing,
}

impl fmt::Display for Crossing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (e, p, [a, b]) = (self.tolerance, self.at.point * 1e3, self.faces);
        write!(
            f,
            "faces {a} and {b} cross: their planes' common line passes ({:.6}, {:.6}, {:.6}) mm, \
             {:.3} tolerances inside face {a} and {:.3} inside face {b}; planes at {:.4}° \
             (normals {}), {:.3} tolerances from the nearest vertex, the faces sharing {}; \
             {} pairs cross",
            p.x,
            p.y,
            p.z,
            self.at.depths[0] / e,
            self.at.depths[1] / e,
            self.at.angle.to_degrees(),
            if self.at.same_way { "alike" } else { "opposed" },
            self.vertex / e,
            match self.shared {
                Shared::Edge => "an edge",
                Shared::Vertex => "a vertex",
                Shared::Nothing => "nothing",
            },
            self.pairs
        )
    }
}

/// Where two flat faces of `r` cross each other, if any do: a point of
/// the common line of their planes that lies in both faces more than the
/// tolerance ε from each face's boundary, measured in the face's plane.
/// The deepest such pair is reported, with how many pairs cross. The
/// boundaries are read from `r`'s mesh made with `chord` (a curved edge
/// is its polyline; a straight edge is the segment between its vertices,
/// whatever the chord). Curved faces are not checked: a crossing with a
/// curved face is not found. Planes within the angular tolerance of
/// parallel are not paired ([`PARALLEL`]).
///
/// The rule, pair by pair, is [`pair_crossing`]'s: exact for straight
/// edges up to rounding (about 1e-16 of the coordinates); a curved edge's
/// polyline is within `chord` of the curve, so where one bounds a face a
/// crossing is found or missed within `chord` of its true depth.
///
/// Why a valid solid does not trip it: its faces meet only along the
/// edges and at the vertices they share. A point of the common line of
/// two faces' planes that lies in both faces is a point the faces share,
/// so on an edge or a vertex of both, at no distance from either face's
/// boundary; two faces sharing an edge along their common line are
/// reported only if they also cross elsewhere.
///
/// Where the tolerance still lets it trip: a valid solid's vertices may
/// stand off its faces' planes by about ε, which moves the planes' common
/// line off the edge or vertex two faces share. Let the faces meet there
/// at an angle α (the dihedral angle on the side under 180°). A vertex
/// there a off one plane and b off the other lies
/// √(a² + b² + 2ab cos α) / sin α from the common line, and the line runs
/// inside both faces by up to δ cot(α/2) for a = b = δ (δ cot α for one
/// offset δ alone). For α ≥ 90° that is at most δ: no report for δ ≤ ε. So a
/// valid solid trips the check only at a fin or a slit sharper than a
/// right angle (normals nearly opposed) with a vertex off a face's plane,
/// and only about as deep as that vertex lies from the common line; a
/// report deeper than that is a defect. Run over the kernel's test suites
/// at commit 77a27b7, every report (1441; planes at 0.5° to 51°, vertices
/// up to 2.1ε off) was such a fold, no deeper than its off-plane vertex
/// lies from the common line (to 0.1 %).
pub fn faces_cross(k: &BrepKernel, r: &Body, chord: f64) -> Result<Option<String>, KernelError> {
    Ok(faces_crossing(k, r, chord)?.map(|c| c.to_string()))
}

/// [`faces_cross`]'s finding as data.
pub fn faces_crossing(
    k: &BrepKernel,
    r: &Body,
    chord: f64,
) -> Result<Option<Crossing>, KernelError> {
    let opts = TessellationOptions {
        chord_tolerance: chord,
        max_angle: 0.5,
    };
    let mesh = k.tessellate(r, &opts)?;
    Ok(faces_cross_in(k, r, &mesh))
}

/// [`faces_crossing`] on `r`'s mesh `mesh`. A face whose boundary does not
/// close (its crossings with a common line are odd in number) is reported
/// as crossing at no depth, with `pairs` 0, as no valid body has one.
pub fn faces_cross_in(k: &BrepKernel, r: &Body, mesh: &BodyMesh) -> Option<Crossing> {
    let eps = k.tolerance().linear;
    let adj = k.adjacency(r);
    let kinds = k.edge_kinds(r);
    let mut ends: Vec<Vec<usize>> = vec![Vec::new(); adj.edge_faces.len()];
    for (v, edges) in adj.vertex_edges.iter().enumerate() {
        for e in edges {
            ends[e.0 as usize].push(v);
        }
    }
    let mut flat: Vec<Option<FlatFace>> = (0..mesh.faces.len())
        .map(|f| {
            k.face_plane(r, FaceIx(f as u32)).map(|pl| FlatFace {
                point: pl.point,
                normal: pl.normal,
                segments: Vec::new(),
            })
        })
        .collect();
    // Each face's vertices and edges, for what two faces share.
    let mut face_edges: Vec<Vec<usize>> = vec![Vec::new(); mesh.faces.len()];
    let mut face_vertices: Vec<Vec<usize>> = vec![Vec::new(); mesh.faces.len()];
    for (e, faces) in adj.edge_faces.iter().enumerate() {
        let line = &mesh.edges[e];
        let segments: Vec<(DVec3, DVec3)> = match (kinds[e], ends[e].as_slice()) {
            (CurveKind::Line, &[v, w]) => vec![(mesh.vertices[v], mesh.vertices[w])],
            _ => {
                // The polyline, its ends on the edge's vertices (the mesh
                // pins them there; this keeps the loops closed exactly).
                let mut pts = line.clone();
                let snap = |p: DVec3| {
                    (ends[e].iter())
                        .map(|&v| mesh.vertices[v])
                        .min_by(|a, b| (*a - p).length().total_cmp(&(*b - p).length()))
                        .unwrap_or(p)
                };
                if let (Some(first), Some(last)) = (pts.first().copied(), pts.last().copied()) {
                    let n = pts.len();
                    pts[0] = snap(first);
                    pts[n - 1] = snap(last);
                }
                pts.windows(2).map(|w| (w[0], w[1])).collect()
            }
        };
        for f in faces {
            let f = f.0 as usize;
            face_edges[f].push(e);
            face_vertices[f].extend(&ends[e]);
            if let Some(face) = flat[f].as_mut() {
                face.segments.extend(&segments);
            }
        }
    }
    let bounds: Vec<Option<(DVec3, DVec3)>> = (flat.iter())
        .map(|face| {
            let face = face.as_ref()?;
            let pts = face.segments.iter().flat_map(|&(a, b)| [a, b]);
            let lo = pts.clone().fold(DVec3::INFINITY, |m, p| m.min(p));
            let hi = pts.fold(DVec3::NEG_INFINITY, |m, p| m.max(p));
            Some((lo - DVec3::splat(eps), hi + DVec3::splat(eps)))
        })
        .collect();
    let mut best: Option<Crossing> = None;
    let mut pairs = 0;
    for f in 0..flat.len() {
        let (Some(ff), Some((flo, fhi))) = (&flat[f], bounds[f]) else {
            continue;
        };
        for g in f + 1..flat.len() {
            let (Some(gg), Some((glo, ghi))) = (&flat[g], bounds[g]) else {
                continue;
            };
            if fhi.cmplt(glo).any() || ghi.cmplt(flo).any() {
                continue;
            }
            let at = match pair_crossing(ff, gg, eps) {
                Ok(Some(at)) => at,
                Ok(None) => continue,
                Err(_) => {
                    return Some(Crossing {
                        faces: [f, g],
                        at: PairCrossing {
                            point: DVec3::ZERO,
                            depths: [0.0; 2],
                            depth: 0.0,
                            angle: 0.0,
                            same_way: false,
                        },
                        vertex: 0.0,
                        shared: Shared::Nothing,
                        pairs: 0,
                        tolerance: eps,
                    })
                }
            };
            pairs += 1;
            if best.as_ref().is_some_and(|b| b.at.depth >= at.depth) {
                continue;
            }
            let shared = if face_edges[f].iter().any(|e| face_edges[g].contains(e)) {
                Shared::Edge
            } else if face_vertices[f]
                .iter()
                .any(|v| face_vertices[g].contains(v))
            {
                Shared::Vertex
            } else {
                Shared::Nothing
            };
            let vertex = (mesh.vertices.iter())
                .map(|v| (*v - at.point).length())
                .fold(f64::INFINITY, f64::min);
            best = Some(Crossing {
                faces: [f, g],
                at,
                vertex,
                shared,
                pairs: 0,
                tolerance: eps,
            });
        }
    }
    best.map(|b| Crossing { pairs, ..b })
}

/// The interval of the line `y = 0` within `r` of the segment `a b`, in
/// `x`: the line meets the segment's `r`-neighbourhood (a stadium, which
/// is convex) in one interval, the hull of where it meets the discs about
/// the ends and the rectangle along the segment.
fn near_segment(a: DVec2, b: DVec2, r: f64) -> Option<(f64, f64)> {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in [a, b] {
        if p.y.abs() <= r {
            let h = (r * r - p.y * p.y).sqrt();
            lo = lo.min(p.x - h);
            hi = hi.max(p.x + h);
        }
    }
    let u = b - a;
    let len2 = u.length_squared();
    if len2 > 0.0 {
        // α x + β within [lo, hi], solved for x.
        let solve = |alpha: f64, beta: f64, l: f64, h: f64| -> Option<(f64, f64)> {
            if alpha == 0.0 {
                return (l <= beta && beta <= h).then_some((f64::NEG_INFINITY, f64::INFINITY));
            }
            let (x, y) = ((l - beta) / alpha, (h - beta) / alpha);
            Some((x.min(y), x.max(y)))
        };
        // The foot of (x, 0) on the segment's line within the segment, and
        // the point within r of that line.
        let along = solve(u.x, -a.x * u.x - a.y * u.y, 0.0, len2);
        let rl = r * len2.sqrt();
        let across = solve(-u.y, a.x * u.y - a.y * u.x, -rl, rl);
        if let (Some((l1, h1)), Some((l2, h2))) = (along, across) {
            let (l, h) = (l1.max(l2), h1.min(h2));
            if l <= h {
                lo = lo.min(l);
                hi = hi.max(h);
            }
        }
    }
    (lo <= hi).then_some((lo, hi))
}

/// The distance from `(x, 0)` to the nearest of `segments`.
fn boundary_distance(x: f64, segments: &[(DVec2, DVec2)]) -> f64 {
    let p = DVec2::new(x, 0.0);
    (segments.iter())
        .map(|&(a, b)| {
            let u = b - a;
            let len2 = u.length_squared();
            let s = if len2 > 0.0 {
                ((p - a).dot(u) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            (p - (a + u * s)).length()
        })
        .fold(f64::INFINITY, f64::min)
}

/// Where the common line of the planes of the flat faces `f` and `g`
/// passes through both more than `eps` from each one's boundary, if it
/// does, at its deepest; `Err` where a face's boundary crosses the line an
/// odd number of times (a loop that does not close). Planes within
/// [`PARALLEL`] of parallel give `None`.
///
/// The common line L is solved from the two planes (refined once by its
/// residuals, so it lies on both up to rounding). Each face's boundary is
/// taken into its plane's coordinates along and across L. The part of L
/// inside a face is found from the boundary's crossings with L by the
/// even–odd rule, a boundary point on L counted on the side `across ≥ 0`
/// (as if L were moved off it by an infinitesimal): with each loop closed
/// exactly, it crosses L an even number of times, and the parity is right
/// at every point of L off the boundary. The common part is the part
/// inside both faces.
///
/// A point of L is within r of a segment of a boundary over one interval
/// of L ([`near_segment`]), so the points of the common part more than r
/// from both boundaries are what is left of it when one interval per
/// segment of either face is taken away: an exact test, by intervals, with
/// no sampling, for any r. The faces cross if anything is left at
/// r = `eps`. The depth reported, the largest r that leaves something, is
/// found by bisecting r with the same test (to a millionth of `eps`), and
/// the point is the middle of the longest piece left at that r; its two
/// depths are measured directly.
pub fn pair_crossing(f: &FlatFace, g: &FlatFace, eps: f64) -> Result<Option<PairCrossing>, String> {
    let d = f.normal.cross(g.normal);
    let s = d.length();
    if s <= PARALLEL {
        return Ok(None);
    }
    let dir = d / s;
    let solve = |rf: f64, rg: f64| (g.normal.cross(d) * rf + d.cross(f.normal) * rg) / (s * s);
    let mut p0 = solve(f.normal.dot(f.point), g.normal.dot(g.point));
    p0 += solve(f.normal.dot(f.point - p0), g.normal.dot(g.point - p0));
    // Each face's boundary in its plane, along (x) and across (y) L.
    let local = |face: &FlatFace| -> Vec<(DVec2, DVec2)> {
        let across = face.normal.cross(dir);
        let q = |p: DVec3| DVec2::new((p - p0).dot(dir), (p - p0).dot(across));
        face.segments.iter().map(|&(a, b)| (q(a), q(b))).collect()
    };
    let (sf, sg) = (local(f), local(g));
    // The part of L inside a face: between the 1st and 2nd crossing, the
    // 3rd and 4th, and so on.
    let inside = |segments: &[(DVec2, DVec2)]| -> Result<Vec<(f64, f64)>, String> {
        let mut xs: Vec<f64> = (segments.iter())
            .filter(|(a, b)| (a.y >= 0.0) != (b.y >= 0.0))
            .map(|(a, b)| a.x + (b.x - a.x) * (a.y / (a.y - b.y)))
            .collect();
        if xs.len() % 2 == 1 {
            return Err(format!(
                "a face's boundary crosses the common line {} times",
                xs.len()
            ));
        }
        xs.sort_by(f64::total_cmp);
        Ok(xs.chunks(2).map(|c| (c[0], c[1])).collect())
    };
    let (inf, ing) = (inside(&sf)?, inside(&sg)?);
    let mut common: Vec<(f64, f64)> = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < inf.len() && j < ing.len() {
        let (lo, hi) = (inf[i].0.max(ing[j].0), inf[i].1.min(ing[j].1));
        if lo < hi {
            common.push((lo, hi));
        }
        if inf[i].1 < ing[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    if common.is_empty() {
        return Ok(None);
    }
    // The middle of the longest piece of the common part more than r from
    // both boundaries, if any is left.
    let left = |r: f64| -> Option<f64> {
        let mut near: Vec<(f64, f64)> = (sf.iter().chain(&sg))
            .filter_map(|&(a, b)| near_segment(a, b, r))
            .collect();
        near.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut best: Option<(f64, f64)> = None;
        let mut take = |lo: f64, hi: f64| {
            if lo < hi && best.is_none_or(|(l, h)| hi - lo > h - l) {
                best = Some((lo, hi));
            }
        };
        for &(c0, c1) in &common {
            let mut at = c0;
            for &(lo, hi) in &near {
                if lo >= c1 {
                    break;
                }
                if lo > at {
                    take(at, lo);
                }
                at = at.max(hi);
                if at >= c1 {
                    break;
                }
            }
            if at < c1 {
                take(at, c1);
            }
        }
        best.map(|(lo, hi)| 0.5 * (lo + hi))
    };
    let Some(mut x) = left(eps) else {
        return Ok(None);
    };
    // The deepest: bisect r between eps and half the longest common piece
    // (each piece ends on a boundary).
    let (mut lo, mut hi) = (
        eps,
        common.iter().map(|c| 0.5 * (c.1 - c.0)).fold(eps, f64::max),
    );
    for _ in 0..100 {
        if hi - lo <= 1e-6 * eps {
            break;
        }
        let mid = 0.5 * (lo + hi);
        match left(mid) {
            Some(t) => {
                lo = mid;
                x = t;
            }
            None => hi = mid,
        }
    }
    let depths = [boundary_distance(x, &sf), boundary_distance(x, &sg)];
    let depth = depths[0].min(depths[1]);
    if depth <= eps {
        return Ok(None);
    }
    let cos = f.normal.dot(g.normal);
    Ok(Some(PairCrossing {
        point: p0 + dir * x,
        depths,
        depth,
        angle: s.atan2(cos.abs()),
        same_way: cos > 0.0,
    }))
}
