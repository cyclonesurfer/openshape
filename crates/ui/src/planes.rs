//! The default geometry of a Part Studio (the three standard planes through
//! the origin): their frames, display meshes and picking.

use openshape_core_model::{DVec3, PlaneRef};
use openshape_kernel::Frame;
use openshape_render::MeshData;

/// One of the three default planes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StdPlane {
    /// World XY (normal +Z).
    Top,
    /// World XZ (normal −Y).
    Front,
    /// World YZ (normal +X).
    Right,
}

impl StdPlane {
    /// All three, in feature-list order.
    pub const ALL: [StdPlane; 3] = [StdPlane::Top, StdPlane::Front, StdPlane::Right];

    /// Display name.
    pub fn label(self) -> &'static str {
        match self {
            StdPlane::Top => "Top",
            StdPlane::Front => "Front",
            StdPlane::Right => "Right",
        }
    }

    /// The model's reference to this plane.
    pub fn plane_ref(self) -> PlaneRef {
        match self {
            StdPlane::Top => PlaneRef::Top,
            StdPlane::Front => PlaneRef::Front,
            StdPlane::Right => PlaneRef::Right,
        }
    }

    /// The plane's frame.
    pub fn frame(self) -> Frame {
        match self {
            StdPlane::Top => Frame::top(),
            StdPlane::Front => Frame::front(),
            StdPlane::Right => Frame::right(),
        }
    }

    /// Index in [`StdPlane::ALL`].
    pub fn index(self) -> usize {
        self as usize
    }
}

/// Half the side of the displayed planes: large enough to frame the
/// sketches (which lie on planes) with a margin, and 0.6 of the parts'
/// farthest coordinate from the origin, so the planes mark the origin's
/// neighbourhood without covering the model; at least 25 mm.
pub fn half_size(
    parts: Option<([f32; 3], [f32; 3])>,
    sketches: Option<([f32; 3], [f32; 3])>,
) -> f64 {
    let reach = |b: Option<([f32; 3], [f32; 3])>| {
        b.map_or(0.0, |(lo, hi)| {
            lo.iter()
                .chain(hi.iter())
                .map(|v| f64::from(v.abs()))
                .fold(0.0, f64::max)
        })
    };
    (reach(parts) * 0.6).max(reach(sketches) * 1.15).max(0.025)
}

/// The four corners of a plane's square of half side `h`, counter-clockwise.
pub fn corners(plane: StdPlane, h: f64) -> [DVec3; 4] {
    let f = plane.frame();
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .map(|(u, v)| f.origin + f.x * (u * h) + f.y * (v * h))
}

/// A translucent square with an outline for the renderer (sRGB colours,
/// linear alpha).
pub fn mesh(plane: StdPlane, h: f64, fill: [f32; 4], outline: [f32; 4]) -> MeshData {
    let n = plane.frame().z.as_vec3().to_array();
    let c = corners(plane, h).map(|p| p.as_vec3().to_array());
    MeshData {
        positions: c.to_vec(),
        normals: vec![n; 4],
        indices: vec![0, 1, 2, 0, 2, 3],
        edge_lines: (0..4).flat_map(|i| [c[i], c[(i + 1) % 4]]).collect(),
        color: fill,
        edge_color: Some(outline),
    }
}

/// True if the plane is seen nearly edge-on when looking along `forward`.
pub fn edge_on(plane: StdPlane, forward: DVec3) -> bool {
    plane.frame().z.dot(forward.normalize_or_zero()).abs() < 0.12
}

/// Where the ray `origin + t·dir` (t > 0) meets the plane's square, as `t`.
pub fn hit(plane: StdPlane, h: f64, origin: DVec3, dir: DVec3) -> Option<f64> {
    let f = plane.frame();
    let denom = dir.dot(f.z);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (f.origin - origin).dot(f.z) / denom;
    if t <= 0.0 {
        return None;
    }
    let p = origin + dir * t - f.origin;
    (p.dot(f.x).abs() <= h && p.dot(f.y).abs() <= h).then_some(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planes_have_the_documented_normals() {
        assert_eq!(StdPlane::Top.frame().z, DVec3::Z);
        assert_eq!(StdPlane::Front.frame().z, -DVec3::Y);
        assert_eq!(StdPlane::Right.frame().z, DVec3::X);
    }

    #[test]
    fn size_follows_the_parts() {
        assert_eq!(half_size(None, None), 0.025);
        let parts = Some(([0.0, -0.1, 0.0], [0.05, 0.02, 0.03]));
        assert!((half_size(parts, None) - 0.06).abs() < 1e-6);
        assert_eq!(half_size(Some(([0.0; 3], [0.01; 3])), None), 0.025);
        let sketch = Some(([-0.04, -0.025, 0.0], [0.04, 0.025, 0.0]));
        assert!((half_size(parts, sketch) - 0.06).abs() < 1e-6);
        assert!((half_size(None, sketch) - 0.046).abs() < 1e-6);
    }

    #[test]
    fn rays_hit_inside_the_square_only() {
        let o = DVec3::new(0.0, 0.0, 1.0);
        assert!((hit(StdPlane::Top, 0.1, o, -DVec3::Z).unwrap() - 1.0).abs() < 1e-12);
        assert!(hit(StdPlane::Top, 0.1, o + DVec3::X * 0.2, -DVec3::Z).is_none());
        assert!(
            hit(StdPlane::Top, 0.1, o, DVec3::Z).is_none(),
            "behind the ray"
        );
        assert!(hit(StdPlane::Top, 0.1, o, DVec3::X).is_none(), "parallel");
        let m = mesh(StdPlane::Right, 0.1, [1.0; 4], [0.5; 4]);
        assert_eq!(m.triangle_count(), 2);
        assert_eq!(m.edge_segment_count(), 4);
        assert!(m.positions.iter().all(|p| p[0] == 0.0));
        assert!(edge_on(StdPlane::Right, DVec3::Y));
        assert!(!edge_on(StdPlane::Right, DVec3::new(1.0, 1.0, 0.0)));
    }
}
