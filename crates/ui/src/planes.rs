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

/// Half the side of the displayed planes for parts within `bounds`: a little
/// beyond the parts' farthest coordinate from the origin, at least 25 mm.
pub fn half_size(bounds: Option<([f32; 3], [f32; 3])>) -> f64 {
    let reach = bounds.map_or(0.0, |(lo, hi)| {
        lo.iter()
            .chain(hi.iter())
            .map(|v| f64::from(v.abs()))
            .fold(0.0, f64::max)
    });
    (reach * 1.1).max(0.025)
}

/// The four corners of a plane's square of half side `h`, counter-clockwise.
pub fn corners(plane: StdPlane, h: f64) -> [DVec3; 4] {
    let f = plane.frame();
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .map(|(u, v)| f.origin + f.x * (u * h) + f.y * (v * h))
}

/// A translucent square for the renderer (sRGB colour, linear alpha).
pub fn mesh(plane: StdPlane, h: f64, color: [f32; 4]) -> MeshData {
    let n = plane.frame().z.as_vec3().to_array();
    MeshData {
        positions: corners(plane, h).map(|p| p.as_vec3().to_array()).to_vec(),
        normals: vec![n; 4],
        indices: vec![0, 1, 2, 0, 2, 3],
        edge_lines: Vec::new(),
        color,
    }
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
        assert_eq!(half_size(None), 0.025);
        let h = half_size(Some(([0.0, -0.1, 0.0], [0.05, 0.02, 0.03])));
        assert!((h - 0.11).abs() < 1e-6);
    }

    #[test]
    fn rays_hit_inside_the_square_only() {
        let o = DVec3::new(0.0, 0.0, 1.0);
        assert!((hit(StdPlane::Top, 0.1, o, -DVec3::Z).unwrap() - 1.0).abs() < 1e-12);
        assert!(hit(StdPlane::Top, 0.1, o + DVec3::X * 0.2, -DVec3::Z).is_none());
        assert!(hit(StdPlane::Top, 0.1, o, DVec3::Z).is_none(), "behind the ray");
        assert!(hit(StdPlane::Top, 0.1, o, DVec3::X).is_none(), "parallel");
        let m = mesh(StdPlane::Right, 0.1, [1.0; 4]);
        assert_eq!(m.triangle_count(), 2);
        assert!(m.positions.iter().all(|p| p[0] == 0.0));
    }
}
