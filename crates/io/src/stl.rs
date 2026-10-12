//! STL export (ASCII and binary).
//!
//! STL has no units; callers choose a scale (OpenDatum stores metres, most
//! slicers expect millimetres, so `MILLIMETRES` is the usual choice).

use std::io::{self, Write};

use opendatum_kernel::{BodyMesh, DVec3};

/// Scale factor from metres to millimetres.
pub const MILLIMETRES: f64 = 1000.0;

fn triangles(meshes: &[BodyMesh]) -> impl Iterator<Item = [DVec3; 3]> + '_ {
    meshes.iter().flat_map(|m| {
        m.faces.iter().flat_map(|f| {
            f.triangles
                .iter()
                .map(move |t| t.map(|i| f.positions[i as usize]))
        })
    })
}

fn facet_normal(t: &[DVec3; 3]) -> DVec3 {
    (t[1] - t[0]).cross(t[2] - t[0]).normalize_or_zero()
}

/// Writes ASCII STL.
pub fn write_ascii<W: Write>(
    out: &mut W,
    name: &str,
    meshes: &[BodyMesh],
    scale: f64,
) -> io::Result<()> {
    let name: String = name.chars().filter(|c| !c.is_whitespace()).collect();
    writeln!(out, "solid {name}")?;
    for t in triangles(meshes) {
        let n = facet_normal(&t);
        writeln!(out, "  facet normal {:e} {:e} {:e}", n.x, n.y, n.z)?;
        writeln!(out, "    outer loop")?;
        for p in t {
            let p = p * scale;
            writeln!(out, "      vertex {:e} {:e} {:e}", p.x, p.y, p.z)?;
        }
        writeln!(out, "    endloop")?;
        writeln!(out, "  endfacet")?;
    }
    writeln!(out, "endsolid {name}")
}

/// Writes binary STL (80-byte header, little-endian f32 records).
pub fn write_binary<W: Write>(out: &mut W, meshes: &[BodyMesh], scale: f64) -> io::Result<()> {
    let count: usize = meshes.iter().map(|m| m.triangle_count()).sum();
    let count = u32::try_from(count)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "too many triangles for STL"))?;
    let mut header = [0u8; 80];
    let tag = b"OpenDatum binary STL";
    header[..tag.len()].copy_from_slice(tag);
    out.write_all(&header)?;
    out.write_all(&count.to_le_bytes())?;
    for t in triangles(meshes) {
        let n = facet_normal(&t);
        for v in [n, t[0] * scale, t[1] * scale, t[2] * scale] {
            for c in [v.x, v.y, v.z] {
                out.write_all(&(c as f32).to_le_bytes())?;
            }
        }
        out.write_all(&0u16.to_le_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use opendatum_kernel::{BrepKernel, Frame, Kernel, TessellationOptions};

    fn cube_mesh() -> BodyMesh {
        let k = BrepKernel::default();
        let b = k.make_box(&Frame::WORLD, DVec3::splat(0.01)).unwrap();
        k.tessellate(&b, &TessellationOptions::default()).unwrap()
    }

    #[test]
    fn binary_layout_and_volume() {
        let m = cube_mesh();
        let mut buf = Vec::new();
        write_binary(&mut buf, std::slice::from_ref(&m), MILLIMETRES).unwrap();
        let n = u32::from_le_bytes(buf[80..84].try_into().unwrap()) as usize;
        assert_eq!(n, 12);
        assert_eq!(buf.len(), 84 + 50 * n);
        // Re-read the facets and integrate the volume: a 10 mm cube.
        let mut vol = 0.0f64;
        for i in 0..n {
            let rec = &buf[84 + 50 * i..84 + 50 * (i + 1)];
            let f = |k: usize| f32::from_le_bytes(rec[4 * k..4 * k + 4].try_into().unwrap()) as f64;
            let p = |j: usize| DVec3::new(f(3 + 3 * j), f(4 + 3 * j), f(5 + 3 * j));
            vol += p(0).dot(p(1).cross(p(2))) / 6.0;
        }
        assert!((vol - 1000.0).abs() < 1e-3);
    }

    #[test]
    fn ascii_has_matching_facets() {
        let m = cube_mesh();
        let mut buf = Vec::new();
        write_ascii(&mut buf, "my cube", std::slice::from_ref(&m), 1.0).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("solid mycube\n"));
        assert!(s.trim_end().ends_with("endsolid mycube"));
        assert_eq!(s.matches("facet normal").count(), 12);
        assert_eq!(s.matches("vertex").count(), 36);
    }
}
