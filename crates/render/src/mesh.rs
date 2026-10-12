//! Neutral triangle-mesh data consumed by the renderer.
//!
//! The renderer never sees kernel types. Whoever tessellates a body (the CAD
//! kernel, an STL importer, a test) produces a [`MeshData`] and hands it to
//! [`crate::SceneRenderer::set_meshes`].

use glam::{Mat3, Mat4, Vec3};
use serde::{Deserialize, Serialize};

/// A displayable body: a triangle mesh plus its B-rep edges as line segments.
///
/// Conventions:
/// - Coordinates are in model space (OpenDatum uses metres, Z up), but the
///   renderer is unit-agnostic: the camera derives its clip planes from the
///   scene bounds, so 1 mm parts and 10 m assemblies both render correctly.
/// - Triangles are counter-clockwise when seen from outside. Back faces are
///   not culled (open shells and sheet bodies are common in CAD), and
///   lighting is two-sided.
/// - [`MeshData::color`] is **sRGB-encoded** RGB (as in colour pickers) plus
///   linear alpha. Alpha below 1 renders the body translucent.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct MeshData {
    /// Vertex positions.
    pub positions: Vec<[f32; 3]>,
    /// Per-vertex unit normals, one per entry of [`MeshData::positions`].
    ///
    /// For smooth shading share vertices between triangles; for faceted
    /// (flat) shading duplicate vertices per face. If this does not have the
    /// same length as `positions`, the renderer computes area-weighted smooth
    /// normals at upload time (see [`MeshData::compute_smooth_normals`]).
    pub normals: Vec<[f32; 3]>,
    /// Triangle-list indices into [`MeshData::positions`], three per triangle.
    pub indices: Vec<u32>,
    /// B-rep edges as a line list: consecutive vertex pairs `[a0, b0, a1, b1, ...]`.
    ///
    /// These are drawn as thin dark lines on top of the shaded surface.
    pub edge_lines: Vec<[f32; 3]>,
    /// Body colour: sRGB-encoded RGB and linear alpha, each in `0.0..=1.0`.
    pub color: [f32; 4],
    /// Colour of this mesh's edge lines (sRGB-encoded RGB and linear alpha),
    /// or `None` for the scene's edge colour
    /// ([`crate::SceneRenderer::set_edge_color`]).
    #[serde(default)]
    pub edge_color: Option<[f32; 4]>,
}

impl MeshData {
    /// A neutral light blue-grey, the default body colour (sRGB, opaque).
    pub const DEFAULT_COLOR: [f32; 4] = [0.70, 0.75, 0.82, 1.0];

    /// Axis-aligned bounding box `(min, max)` over all finite triangle
    /// vertices and edge-line vertices, or `None` if there are none.
    pub fn bounds(&self) -> Option<([f32; 3], [f32; 3])> {
        let mut min = Vec3::splat(f32::INFINITY);
        let mut max = Vec3::splat(f32::NEG_INFINITY);
        let mut any = false;
        for p in self.positions.iter().chain(self.edge_lines.iter()) {
            let v = Vec3::from_array(*p);
            if v.is_finite() {
                min = min.min(v);
                max = max.max(v);
                any = true;
            }
        }
        any.then(|| (min.to_array(), max.to_array()))
    }

    /// Number of complete triangles in [`MeshData::indices`].
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Number of complete edge segments in [`MeshData::edge_lines`].
    pub fn edge_segment_count(&self) -> usize {
        self.edge_lines.len() / 2
    }

    /// True if the mesh has neither triangles nor edge segments.
    pub fn is_empty(&self) -> bool {
        self.triangle_count() == 0 && self.edge_segment_count() == 0
    }

    /// True if [`MeshData::normals`] has one entry per position.
    pub fn has_normals(&self) -> bool {
        self.normals.len() == self.positions.len()
    }

    /// Replaces [`MeshData::normals`] with area-weighted vertex normals
    /// computed from the triangles. Vertices not used by any triangle get
    /// `+Z`. Out-of-range indices are ignored.
    pub fn compute_smooth_normals(&mut self) {
        let mut acc = vec![Vec3::ZERO; self.positions.len()];
        for tri in self.indices.chunks_exact(3) {
            let [a, b, c] = [tri[0] as usize, tri[1] as usize, tri[2] as usize];
            if a >= acc.len() || b >= acc.len() || c >= acc.len() {
                continue;
            }
            let pa = Vec3::from_array(self.positions[a]);
            let pb = Vec3::from_array(self.positions[b]);
            let pc = Vec3::from_array(self.positions[c]);
            // Cross product length is twice the area: area weighting for free.
            let n = (pb - pa).cross(pc - pa);
            if n.is_finite() {
                acc[a] += n;
                acc[b] += n;
                acc[c] += n;
            }
        }
        self.normals = acc
            .into_iter()
            .map(|n| n.try_normalize().unwrap_or(Vec3::Z).to_array())
            .collect();
    }

    /// Returns a copy with positions and edges transformed by `m` and normals
    /// by its inverse transpose (renormalised).
    pub fn transformed(&self, m: &Mat4) -> MeshData {
        let normal_matrix = Mat3::from_mat4(*m).inverse().transpose();
        let xf = |p: &[f32; 3]| m.transform_point3(Vec3::from_array(*p)).to_array();
        MeshData {
            positions: self.positions.iter().map(xf).collect(),
            normals: self
                .normals
                .iter()
                .map(|n| {
                    (normal_matrix * Vec3::from_array(*n))
                        .try_normalize()
                        .unwrap_or(Vec3::Z)
                        .to_array()
                })
                .collect(),
            indices: self.indices.clone(),
            edge_lines: self.edge_lines.iter().map(xf).collect(),
            color: self.color,
            edge_color: self.edge_color,
        }
    }

    /// An axis-aligned box from `min` to `max` with flat-shaded faces:
    /// 24 vertices (4 per face, each with the face normal), 12 triangles and
    /// its 12 edges as line segments.
    pub fn cuboid(min: [f32; 3], max: [f32; 3], color: [f32; 4]) -> MeshData {
        let (lo, hi) = (Vec3::from_array(min), Vec3::from_array(max));
        let corner = |x: bool, y: bool, z: bool| -> [f32; 3] {
            [
                if x { hi.x } else { lo.x },
                if y { hi.y } else { lo.y },
                if z { hi.z } else { lo.z },
            ]
        };
        // Each face: outward normal and its 4 corners counter-clockwise when
        // seen from outside.
        let faces: [([f32; 3], [[bool; 3]; 4]); 6] = [
            (
                [1.0, 0.0, 0.0],
                [
                    [true, false, false],
                    [true, true, false],
                    [true, true, true],
                    [true, false, true],
                ],
            ),
            (
                [-1.0, 0.0, 0.0],
                [
                    [false, true, false],
                    [false, false, false],
                    [false, false, true],
                    [false, true, true],
                ],
            ),
            (
                [0.0, 1.0, 0.0],
                [
                    [true, true, false],
                    [false, true, false],
                    [false, true, true],
                    [true, true, true],
                ],
            ),
            (
                [0.0, -1.0, 0.0],
                [
                    [false, false, false],
                    [true, false, false],
                    [true, false, true],
                    [false, false, true],
                ],
            ),
            (
                [0.0, 0.0, 1.0],
                [
                    [false, false, true],
                    [true, false, true],
                    [true, true, true],
                    [false, true, true],
                ],
            ),
            (
                [0.0, 0.0, -1.0],
                [
                    [false, true, false],
                    [true, true, false],
                    [true, false, false],
                    [false, false, false],
                ],
            ),
        ];
        let mut mesh = MeshData {
            color,
            ..MeshData::default()
        };
        for (normal, quad) in faces {
            let base = mesh.positions.len() as u32;
            for [x, y, z] in quad {
                mesh.positions.push(corner(x, y, z));
                mesh.normals.push(normal);
            }
            mesh.indices
                .extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        }
        // 12 edges: for each axis, the 4 edges parallel to it.
        for axis in 0..3 {
            for (a, b) in [(false, false), (true, false), (false, true), (true, true)] {
                let mut start = [false; 3];
                let mut end = [true; 3];
                let (o1, o2) = ((axis + 1) % 3, (axis + 2) % 3);
                start[o1] = a;
                start[o2] = b;
                end[o1] = a;
                end[o2] = b;
                mesh.edge_lines.push(corner(start[0], start[1], start[2]));
                mesh.edge_lines.push(corner(end[0], end[1], end[2]));
            }
        }
        mesh
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cuboid_counts_and_bounds() {
        let m = MeshData::cuboid([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0], MeshData::DEFAULT_COLOR);
        assert_eq!(m.positions.len(), 24);
        assert_eq!(m.normals.len(), 24);
        assert_eq!(m.triangle_count(), 12);
        assert_eq!(m.edge_segment_count(), 12);
        assert_eq!(m.bounds(), Some(([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0])));
        // Every triangle's winding agrees with its vertex normals (outward).
        for tri in m.indices.chunks_exact(3) {
            let p: Vec<Vec3> = tri
                .iter()
                .map(|&i| Vec3::from_array(m.positions[i as usize]))
                .collect();
            let n = (p[1] - p[0]).cross(p[2] - p[0]).normalize();
            let vn = Vec3::from_array(m.normals[tri[0] as usize]);
            assert!(n.dot(vn) > 0.99, "winding disagrees with normal");
        }
        // Every edge is axis-parallel with full length.
        for seg in m.edge_lines.chunks_exact(2) {
            let d = Vec3::from_array(seg[1]) - Vec3::from_array(seg[0]);
            let nonzero = d.to_array().iter().filter(|c| c.abs() > 0.0).count();
            assert_eq!(nonzero, 1);
        }
    }

    #[test]
    fn empty_mesh_has_no_bounds() {
        assert_eq!(MeshData::default().bounds(), None);
        assert!(MeshData::default().is_empty());
    }

    #[test]
    fn smooth_normals_of_cube_point_outward() {
        let mut m = MeshData::cuboid([0.0; 3], [1.0; 3], MeshData::DEFAULT_COLOR);
        m.normals.clear();
        m.compute_smooth_normals();
        assert!(m.has_normals());
        for (p, n) in m.positions.iter().zip(&m.normals) {
            let outward = Vec3::from_array(*p) - Vec3::splat(0.5);
            assert!(Vec3::from_array(*n).dot(outward) > 0.0);
        }
    }

    #[test]
    fn transformed_moves_bounds_and_keeps_unit_normals() {
        let m = MeshData::cuboid([0.0; 3], [1.0; 3], MeshData::DEFAULT_COLOR);
        let t = m.transformed(&Mat4::from_scale_rotation_translation(
            Vec3::new(2.0, 1.0, 1.0),
            glam::Quat::IDENTITY,
            Vec3::new(10.0, 0.0, 0.0),
        ));
        assert_eq!(t.bounds(), Some(([10.0, 0.0, 0.0], [12.0, 1.0, 1.0])));
        for n in &t.normals {
            assert!((Vec3::from_array(*n).length() - 1.0).abs() < 1e-5);
        }
    }
}
