//! Picking parts in the 3D view: the nearest display-mesh triangle along a
//! view ray (Möller and Trumbore, "Fast, minimum storage ray/triangle
//! intersection", Journal of Graphics Tools 2(1), 1997).

use openshape_core_model::DVec3;
use openshape_render::MeshData;

/// Distance along `dir` from `origin` to the triangle `abc`, if the ray
/// (t > 0) crosses it. Both sides count.
pub fn ray_triangle(origin: DVec3, dir: DVec3, a: DVec3, b: DVec3, c: DVec3) -> Option<f64> {
    let e1 = b - a;
    let e2 = c - a;
    let p = dir.cross(e2);
    let det = e1.dot(p);
    let scale = e1.length() * e2.length() * dir.length();
    if det.abs() <= 1e-12 * scale {
        return None;
    }
    let inv = 1.0 / det;
    let s = origin - a;
    let u = s.dot(p) * inv;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = dir.dot(q) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(q) * inv;
    (t > 0.0).then_some(t)
}

/// The mesh hit first by the ray and the distance to it.
pub fn pick_mesh<'a>(
    meshes: impl IntoIterator<Item = &'a MeshData>,
    origin: DVec3,
    dir: DVec3,
) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64)> = None;
    for (i, m) in meshes.into_iter().enumerate() {
        if let Some((lo, hi)) = m.bounds() {
            if !ray_box(origin, dir, lo, hi) {
                continue;
            }
        }
        let vertex = |k: u32| m.positions.get(k as usize).map(|p| DVec3::from(p.map(f64::from)));
        for t in m.indices.chunks_exact(3) {
            let (Some(a), Some(b), Some(c)) = (vertex(t[0]), vertex(t[1]), vertex(t[2])) else {
                continue;
            };
            if let Some(d) = ray_triangle(origin, dir, a, b, c) {
                if best.is_none_or(|(_, bd)| d < bd) {
                    best = Some((i, d));
                }
            }
        }
    }
    best
}

/// True if the ray meets the box (slab test, slightly enlarged).
fn ray_box(origin: DVec3, dir: DVec3, lo: [f32; 3], hi: [f32; 3]) -> bool {
    let (mut t0, mut t1) = (0.0f64, f64::INFINITY);
    for k in 0..3 {
        let pad = 1e-6 * (1.0 + f64::from(hi[k] - lo[k]).abs());
        let (l, h) = (f64::from(lo[k]) - pad, f64::from(hi[k]) + pad);
        if dir[k].abs() < 1e-300 {
            if origin[k] < l || origin[k] > h {
                return false;
            }
            continue;
        }
        let (a, b) = ((l - origin[k]) / dir[k], (h - origin[k]) / dir[k]);
        t0 = t0.max(a.min(b));
        t1 = t1.min(a.max(b));
        if t0 > t1 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_hits_a_triangle_from_both_sides() {
        let (a, b, c) = (DVec3::ZERO, DVec3::X, DVec3::Y);
        let down = ray_triangle(DVec3::new(0.2, 0.2, 1.0), -DVec3::Z, a, b, c).unwrap();
        assert!((down - 1.0).abs() < 1e-12);
        let up = ray_triangle(DVec3::new(0.2, 0.2, -2.0), DVec3::Z, a, b, c).unwrap();
        assert!((up - 2.0).abs() < 1e-12);
        assert!(ray_triangle(DVec3::new(0.8, 0.8, 1.0), -DVec3::Z, a, b, c).is_none());
        assert!(ray_triangle(DVec3::new(0.2, 0.2, 1.0), DVec3::Z, a, b, c).is_none());
        assert!(ray_triangle(DVec3::new(0.2, 0.2, 1.0), DVec3::X, a, b, c).is_none());
    }

    #[test]
    fn the_nearest_mesh_wins() {
        let near = MeshData::cuboid([0.0, 0.0, 1.0], [1.0, 1.0, 2.0], [1.0; 4]);
        let far = MeshData::cuboid([0.0, 0.0, -2.0], [1.0, 1.0, -1.0], [1.0; 4]);
        let o = DVec3::new(0.5, 0.5, 5.0);
        let (i, t) = pick_mesh([&far, &near], o, -DVec3::Z).unwrap();
        assert_eq!(i, 1);
        assert!((t - 3.0).abs() < 1e-9);
        assert!(pick_mesh([&far, &near], DVec3::new(3.0, 3.0, 5.0), -DVec3::Z).is_none());
    }
}
