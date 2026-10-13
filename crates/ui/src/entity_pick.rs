//! Picking the faces, edges and vertices of parts in the 3D view (ADR 0019):
//! of the entities within a few points of the pointer, a vertex wins over
//! an edge, and an edge over the face under the pointer. An edge or vertex
//! hidden behind a face nearer the eye is not picked. Seams, and points
//! inside an edge or on no edge (a sphere's pole), are no entities the user
//! sees (design note 0013 §11), so they are never picked.

use egui::{Pos2, Rect};
use opendatum_core_model::{DVec3, FeatureId};
use opendatum_kernel::{BodyMesh, EdgeIx, FaceIx, Frame, VertexIx};
use opendatum_naming::{BodyNames, Pick};
use opendatum_render::Camera;

use crate::picking::ray_triangle;
use crate::sketch_mode::ScreenMap;

/// How near the pointer a vertex is picked, in screen points.
pub const VERTEX_PICK: f32 = 7.0;
/// How near the pointer an edge is picked, in screen points.
pub const EDGE_PICK: f32 = 5.0;

/// Which kinds of entity a pick may return.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PickKinds {
    /// Faces, edges and vertices.
    All,
    /// Faces only (a face to extrude up to).
    Faces,
}

/// A part's entity under the pointer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EntityHit {
    /// The part, by the feature that made it.
    pub feature: FeatureId,
    /// The kernel entity: one piece of the entity the user sees.
    pub pick: Pick,
    /// The picked point on it, in model coordinates.
    pub point: DVec3,
    /// The face's outward normal or the edge's direction there (unit);
    /// zero for a vertex.
    pub direction: DVec3,
}

impl EntityHit {
    /// The part and entity, without where it was picked.
    pub fn key(&self) -> (FeatureId, Pick) {
        (self.feature, self.pick)
    }
}

/// A part as the picker sees it.
#[derive(Clone, Copy, Debug)]
pub struct PickPart<'a> {
    /// The feature that made it.
    pub feature: FeatureId,
    /// Its kernel mesh.
    pub mesh: &'a BodyMesh,
    /// Its names, when it has them: edges and vertices without names (and
    /// every edge and vertex of a part without names) are not picked.
    pub names: Option<&'a BodyNames>,
}

/// The nearest face along the ray: part index, face, distance, and the
/// hit triangle's outward normal.
fn first_face(parts: &[PickPart], o: DVec3, d: DVec3) -> Option<(usize, FaceIx, f64, DVec3)> {
    let mut best: Option<(usize, FaceIx, f64, DVec3)> = None;
    for (pi, p) in parts.iter().enumerate() {
        for (fi, fm) in p.mesh.faces.iter().enumerate() {
            for t in &fm.triangles {
                let [a, b, c] = t.map(|i| fm.positions[i as usize]);
                if let Some(s) = ray_triangle(o, d, a, b, c) {
                    if best.is_none_or(|(_, _, bs, _)| s < bs) {
                        let n = (b - a).cross(c - a).normalize_or_zero();
                        best = Some((pi, FaceIx(fi as u32), s, n));
                    }
                }
            }
        }
    }
    best
}

/// The parts' size: the largest diagonal of their bounds.
fn size(parts: &[PickPart]) -> f64 {
    parts
        .iter()
        .filter_map(|p| p.mesh.bounds())
        .map(|(lo, hi)| (hi - lo).length())
        .fold(0.0, f64::max)
}

/// Distance from `p` to the segment `ab` on screen, and where along it.
fn screen_segment(p: Pos2, a: Pos2, b: Pos2) -> (f32, f32) {
    let ab = b - a;
    let len2 = ab.length_sq();
    let s = if len2 > 0.0 {
        ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    ((a + ab * s - p).length(), s)
}

/// The entity under `pos` (a point of `rect`, where `camera` draws).
pub fn pick_entity(
    parts: &[PickPart],
    camera: &Camera,
    rect: Rect,
    pos: Pos2,
    kinds: PickKinds,
) -> Option<EntityHit> {
    let ray = |s: Pos2| {
        let (o, d) = camera.pixel_ray(
            s.x - rect.min.x,
            s.y - rect.min.y,
            rect.width(),
            rect.height(),
        );
        (o, d.normalize())
    };
    let (o, d) = ray(pos);
    let face = first_face(parts, o, d).map(|(pi, f, t, n)| EntityHit {
        feature: parts[pi].feature,
        pick: Pick::Face(f),
        point: o + d * t,
        direction: if n.dot(d) > 0.0 { -n } else { n },
    });
    if kinds == PickKinds::Faces {
        return face;
    }
    let map = ScreenMap::new(camera, Frame::WORLD, rect);
    // A point is hidden when a face is nearer the eye along the ray to it,
    // by more than a small fraction of the model (the faces it lies on
    // meet that ray where it is).
    let slack = 1e-3 * size(parts);
    let visible = |q: DVec3| {
        let Some(s) = map.world_to_screen(q) else {
            return false;
        };
        let (qo, qd) = ray(s);
        let depth = (q - qo).dot(qd);
        first_face(parts, qo, qd).is_none_or(|(_, _, t, _)| t >= depth - slack)
    };
    // Vertices, nearest first.
    let mut vertices: Vec<(f32, EntityHit)> = Vec::new();
    for p in parts {
        let Some(names) = p.names else { continue };
        for (vi, q) in p.mesh.vertices.iter().enumerate() {
            let v = VertexIx(vi as u32);
            if names.vertex(v).is_none() {
                continue;
            }
            if let Some(s) = map.world_to_screen(*q) {
                let dist = (s - pos).length();
                if dist <= VERTEX_PICK {
                    vertices.push((
                        dist,
                        EntityHit {
                            feature: p.feature,
                            pick: Pick::Vertex(v),
                            point: *q,
                            direction: DVec3::ZERO,
                        },
                    ));
                }
            }
        }
    }
    vertices.sort_by(|a, b| a.0.total_cmp(&b.0));
    if let Some((_, hit)) = vertices.into_iter().find(|(_, h)| visible(h.point)) {
        return Some(hit);
    }
    // Edges, nearest first.
    let mut edges: Vec<(f32, EntityHit)> = Vec::new();
    for p in parts {
        let Some(names) = p.names else { continue };
        for (ei, line) in p.mesh.edges.iter().enumerate() {
            let e = EdgeIx(ei as u32);
            if p.mesh.seams.get(ei).copied().unwrap_or(false) || names.edge(e).is_none() {
                continue;
            }
            let mut best: Option<(f32, DVec3, DVec3)> = None;
            for w in line.windows(2) {
                let (Some(a), Some(b)) = (map.world_to_screen(w[0]), map.world_to_screen(w[1]))
                else {
                    continue;
                };
                let (dist, s) = screen_segment(pos, a, b);
                if dist <= EDGE_PICK && best.is_none_or(|(bd, _, _)| dist < bd) {
                    let q = w[0] + (w[1] - w[0]) * f64::from(s);
                    best = Some((dist, q, (w[1] - w[0]).normalize_or_zero()));
                }
            }
            if let Some((dist, point, direction)) = best {
                edges.push((
                    dist,
                    EntityHit {
                        feature: p.feature,
                        pick: Pick::Edge(e),
                        point,
                        direction,
                    },
                ));
            }
        }
    }
    edges.sort_by(|a, b| a.0.total_cmp(&b.0));
    if let Some((_, hit)) = edges.into_iter().find(|(_, h)| visible(h.point)) {
        return Some(hit);
    }
    face
}

/// The kernel entities that make up the entity the user sees as `pick`
/// (all pieces of a face split at seams, of an edge, or the vertex), or
/// just `pick` when the part has no names.
pub fn logical_pieces(names: Option<&BodyNames>, pick: Pick) -> Vec<Pick> {
    let Some(names) = names else {
        return vec![pick];
    };
    match pick {
        Pick::Face(f) => names
            .face(f)
            .pieces
            .iter()
            .map(|&p| Pick::Face(p))
            .collect(),
        Pick::Edge(e) => names.edge(e).map_or(vec![pick], |l| {
            l.pieces.iter().map(|&p| Pick::Edge(p)).collect()
        }),
        Pick::Vertex(v) => names.vertex(v).map_or(vec![pick], |l| {
            l.pieces.iter().map(|&p| Pick::Vertex(p)).collect()
        }),
    }
}
