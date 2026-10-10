//! Linear (extrude) and rotational (revolve) sweeps of planar profiles.
//!
//! Both build B-reps directly from the oriented profile regions, so the
//! topology is known exactly. Periodic surfaces are split so that no face
//! spans more than π in a periodic direction; this keeps every face's
//! parameter domain simple (see ADR 0001 and `docs/design/0001-brep.md`).

use std::f64::consts::{PI, TAU};

use glam::{DMat4, DVec3};

use crate::geom::{Curve3, Surface};
use crate::math::{rotate_vec, rotation_about, Frame};
use crate::profile::{regions, Profile, Region, Seg};
use crate::topo::{Brep, Coedge, EdgeId, Loop, VertexId};
use crate::{KernelError, Tolerance};

/// Extrude parameters.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExtrudeSpec {
    /// Signed distance along the profile plane normal (metres).
    pub distance: f64,
    /// Extrude half the distance to each side of the profile plane.
    pub symmetric: bool,
}

/// Revolve parameters.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RevolveSpec {
    /// A point on the axis, in profile-plane coordinates.
    pub axis_origin: glam::DVec2,
    /// Axis direction, in profile-plane coordinates.
    pub axis_dir: glam::DVec2,
    /// Signed sweep angle in radians, `0 < |angle| ≤ 2π` (right-hand rule
    /// about `axis_dir`).
    pub angle: f64,
}

fn loop_from(coedges: Vec<Coedge>, reverse: bool) -> Loop {
    if reverse {
        Loop {
            coedges: coedges
                .into_iter()
                .rev()
                .map(|c| Coedge {
                    edge: c.edge,
                    forward: !c.forward,
                })
                .collect(),
        }
    } else {
        Loop { coedges }
    }
}

fn fw(edge: EdgeId) -> Coedge {
    Coedge {
        edge,
        forward: true,
    }
}

fn bw(edge: EdgeId) -> Coedge {
    Coedge {
        edge,
        forward: false,
    }
}

/// Builds a solid by extruding a profile.
pub(crate) fn extrude(
    profile: &Profile,
    spec: &ExtrudeSpec,
    tol: &Tolerance,
) -> Result<Brep, KernelError> {
    if !spec.distance.is_finite() || spec.distance.abs() <= tol.linear {
        return Err(KernelError::InvalidInput(
            "extrude distance must be non-zero".into(),
        ));
    }
    let regs = regions(profile, tol)?;
    let plane = profile.plane;
    let d = plane.z * spec.distance;
    let base = if spec.symmetric {
        -d * 0.5
    } else {
        DVec3::ZERO
    };
    let m0 = DMat4::from_translation(base) * plane.matrix();
    let m1 = DMat4::from_translation(base + d) * plane.matrix();
    let mut brep = Brep::default();
    for r in &regs {
        extrude_region(&mut brep, r, &plane, d, base, &m0, &m1);
    }
    Ok(brep)
}

fn extrude_region(
    brep: &mut Brep,
    region: &Region,
    plane: &Frame,
    d: DVec3,
    base: DVec3,
    m0: &DMat4,
    m1: &DMat4,
) {
    let mut faces = Vec::new();
    let mut bottom_loops = Vec::new();
    let mut top_loops = Vec::new();
    let up = d.dot(plane.z) > 0.0;
    for lp in &region.loops {
        let n = lp.len();
        let vb: Vec<VertexId> = lp
            .iter()
            .map(|s| brep.add_vertex(m0.transform_point3(s.start().extend(0.0))))
            .collect();
        let vt: Vec<VertexId> = lp
            .iter()
            .map(|s| brep.add_vertex(m1.transform_point3(s.start().extend(0.0))))
            .collect();
        let eb: Vec<EdgeId> = (0..n)
            .map(|i| {
                brep.add_edge(
                    lp[i].curve.transformed(m0),
                    lp[i].t0,
                    lp[i].t1,
                    vb[i],
                    vb[(i + 1) % n],
                )
            })
            .collect();
        let et: Vec<EdgeId> = (0..n)
            .map(|i| {
                brep.add_edge(
                    lp[i].curve.transformed(m1),
                    lp[i].t0,
                    lp[i].t1,
                    vt[i],
                    vt[(i + 1) % n],
                )
            })
            .collect();
        let ev: Vec<EdgeId> = (0..n)
            .map(|i| {
                let a = brep.vertex(vb[i]).point;
                brep.add_edge(Curve3::Line { origin: a, dir: d }, 0.0, 1.0, vb[i], vt[i])
            })
            .collect();
        for i in 0..n {
            let s = &lp[i];
            let bottom_curve = brep.edge(eb[i]).curve.clone();
            let surface = match &bottom_curve {
                Curve3::Line { origin, dir } => {
                    let o = *origin + *dir * s.t0;
                    let x = dir.normalize();
                    // Parameterize so that u runs along the edge and v along d.
                    let frame = Frame {
                        origin: o - x * (s.t0 * dir.length()),
                        x,
                        y: d.normalize(),
                        z: x.cross(d.normalize()),
                    };
                    let _ = o;
                    Surface::Plane(frame)
                }
                _ => Surface::Extrusion {
                    curve: bottom_curve.clone(),
                    dir: d,
                    base_origin: plane.origin + base,
                    base_normal: plane.z,
                },
            };
            // Outward direction: material is on the left of travel (seen from
            // +plane.z), so outward is tangent × z.
            let tm = 0.5 * (s.t0 + s.t1);
            let tangent = m0.transform_vector3(s.curve.deriv(tm));
            let outward = tangent.cross(plane.z);
            let (u, v) = surface.project(bottom_curve.point(tm) + d * 0.5);
            let normal = surface.normal(u, v).unwrap_or(DVec3::ZERO);
            let forward = normal.dot(outward) > 0.0;
            let lp_side = vec![fw(eb[i]), fw(ev[(i + 1) % n]), bw(et[i]), bw(ev[i])];
            faces.push(brep.add_face(surface, forward, vec![loop_from(lp_side, !forward)]));
        }
        bottom_loops.push(eb.iter().map(|&e| fw(e)).collect::<Vec<_>>());
        top_loops.push(et.iter().map(|&e| fw(e)).collect::<Vec<_>>());
    }
    // Caps. Profile loops are counter-clockwise about plane.z.
    let bottom_frame = plane.translated(base);
    let top_frame = plane.translated(base + d);
    // Bottom outward is -d, top outward is +d.
    let bottom_out_is_z = !up;
    let top_out_is_z = up;
    faces.push(
        brep.add_face(
            Surface::Plane(bottom_frame),
            bottom_out_is_z,
            bottom_loops
                .into_iter()
                .map(|l| loop_from(l, !bottom_out_is_z))
                .collect(),
        ),
    );
    faces.push(
        brep.add_face(
            Surface::Plane(top_frame),
            top_out_is_z,
            top_loops
                .into_iter()
                .map(|l| loop_from(l, !top_out_is_z))
                .collect(),
        ),
    );
    let shell = brep.add_shell(faces);
    brep.add_solid(vec![shell]);
}

/// Builds a solid by revolving a profile about an axis in its plane.
pub(crate) fn revolve(
    profile: &Profile,
    spec: &RevolveSpec,
    tol: &Tolerance,
) -> Result<Brep, KernelError> {
    let theta = spec.angle.abs();
    if !theta.is_finite() || theta <= tol.angular || theta > TAU + 1e-9 {
        return Err(KernelError::InvalidInput(
            "revolve angle must be in (0, 2π]".into(),
        ));
    }
    if spec.axis_dir.length() <= 1e-12
        || !spec.axis_dir.is_finite()
        || !spec.axis_origin.is_finite()
    {
        return Err(KernelError::InvalidInput(
            "revolve axis is degenerate".into(),
        ));
    }
    let full = theta >= TAU - 1e-9;
    let regs = regions(profile, tol)?;
    let plane = profile.plane;
    let o = plane.point2(spec.axis_origin);
    let z = plane.vector2(spec.axis_dir).normalize() * spec.angle.signum();
    let w = z.cross(plane.z).normalize(); // in-plane, perpendicular to the axis
                                          // Which side of the axis is the profile on?
    let on_tol = tol.linear * 10.0;
    let (mut min_s, mut max_s) = (f64::MAX, f64::MIN);
    for r in &regs {
        for l in &r.loops {
            for s in l {
                for p in s.sample(32) {
                    let sd = (plane.point2(p) - o).dot(w);
                    min_s = min_s.min(sd);
                    max_s = max_s.max(sd);
                }
            }
        }
    }
    let x = if min_s >= -on_tol {
        w
    } else if max_s <= on_tol {
        -w
    } else {
        return Err(KernelError::InvalidInput(
            "profile crosses the revolve axis".into(),
        ));
    };
    if max_s.abs().max(min_s.abs()) <= on_tol {
        return Err(KernelError::InvalidInput(
            "profile lies on the revolve axis".into(),
        ));
    }
    let frame0 = Frame {
        origin: o,
        x,
        y: z.cross(x),
        z,
    };
    let k = if full {
        2usize.max((theta / PI - 1e-9).ceil() as usize)
    } else {
        ((theta / PI) - 1e-9).ceil().max(1.0) as usize
    };
    let angles: Vec<f64> = (0..=k).map(|j| theta * j as f64 / k as f64).collect();
    let mut brep = Brep::default();
    for r in &regs {
        revolve_region(&mut brep, r, &plane, &frame0, &angles, full, on_tol)?;
    }
    Ok(brep)
}

fn rho(frame: &Frame, p: DVec3) -> f64 {
    let d = p - frame.origin;
    (d - frame.z * d.dot(frame.z)).length()
}

#[allow(clippy::too_many_arguments)]
fn revolve_region(
    brep: &mut Brep,
    region: &Region,
    plane: &Frame,
    f0: &Frame,
    angles: &[f64],
    full: bool,
    on_tol: f64,
) -> Result<(), KernelError> {
    let k = angles.len() - 1;
    let nj = if full { k } else { k + 1 }; // number of distinct copies
    let copy = |j: usize| if full { j % k } else { j };
    let rot = |j: usize| rotation_about(f0.origin, f0.z, angles[j]);
    let pm = plane.matrix();

    let mut faces = Vec::new();
    let mut start_cap: Vec<Vec<Coedge>> = Vec::new();
    let mut end_cap: Vec<Vec<Coedge>> = Vec::new();
    for lp in &region.loops {
        let n = lp.len();
        let starts: Vec<DVec3> = lp.iter().map(|s| plane.point2(s.start())).collect();
        let on_axis: Vec<bool> = starts.iter().map(|p| rho(f0, *p) <= on_tol).collect();
        // Segments lying on the axis (straight lines with both ends on it).
        let axis_seg: Vec<bool> = (0..n)
            .map(|i| {
                on_axis[i] && on_axis[(i + 1) % n] && matches!(lp[i].curve, Curve3::Line { .. })
            })
            .collect();
        // Reject curves touching the axis between their ends.
        for (i, s) in lp.iter().enumerate() {
            if axis_seg[i] {
                continue;
            }
            for (q, p) in s.sample(64).iter().enumerate() {
                if q == 0 || q == 64 {
                    continue;
                }
                if rho(f0, plane.point2(*p)) <= on_tol {
                    return Err(KernelError::InvalidInput(
                        "profile touches the revolve axis away from a vertex".into(),
                    ));
                }
            }
        }
        // Vertices: per copy for off-axis points, shared for on-axis points.
        let verts: Vec<Vec<VertexId>> = (0..n)
            .map(|i| {
                if on_axis[i] {
                    let v = brep.add_vertex(starts[i]);
                    vec![v; nj]
                } else {
                    (0..nj)
                        .map(|j| brep.add_vertex(rot(j).transform_point3(starts[i])))
                        .collect()
                }
            })
            .collect();
        let vtx = |i: usize, j: usize| verts[i % n][copy(j)];
        // Profile copies (and single axis edges for partial revolves).
        let mut copies: Vec<Vec<Option<EdgeId>>> = vec![vec![None; nj]; n];
        for i in 0..n {
            let base = lp[i].curve.transformed(&pm);
            if axis_seg[i] {
                if !full {
                    let e = brep.add_edge(base, lp[i].t0, lp[i].t1, vtx(i, 0), vtx(i + 1, 0));
                    copies[i] = vec![Some(e); nj];
                }
                continue;
            }
            for (j, slot) in copies[i].iter_mut().enumerate() {
                let c = base.transformed(&rot(j));
                *slot = Some(brep.add_edge(c, lp[i].t0, lp[i].t1, vtx(i, j), vtx(i + 1, j)));
            }
        }
        // Circular edges swept by off-axis vertices, one per piece.
        let mut arcs: Vec<Vec<Option<EdgeId>>> = vec![vec![None; k]; n];
        for i in 0..n {
            if on_axis[i] {
                continue;
            }
            for j in 0..k {
                let p = rot(j).transform_point3(starts[i]);
                let center = f0.origin + f0.z * (p - f0.origin).dot(f0.z);
                let r = (p - center).length();
                let cx = (p - center) / r;
                let curve = Curve3::Circle {
                    center,
                    x: cx,
                    y: f0.z.cross(cx),
                    radius: r,
                };
                arcs[i][j] = Some(brep.add_edge(
                    curve,
                    0.0,
                    angles[j + 1] - angles[j],
                    vtx(i, j),
                    vtx(i, j + 1),
                ));
            }
        }
        // Side faces.
        for i in 0..n {
            if axis_seg[i] {
                continue;
            }
            let s: &Seg = &lp[i];
            let tm = 0.5 * (s.t0 + s.t1);
            let tangent2 = s.curve.deriv(tm).truncate();
            let tangent = plane.vector2(tangent2);
            let outward0 = tangent.cross(plane.z);
            for j in 0..k {
                let rj = rot(j);
                let frame = Frame {
                    origin: f0.origin,
                    x: rj.transform_vector3(f0.x),
                    y: rj.transform_vector3(f0.y),
                    z: f0.z,
                };
                let profile_curve = brep
                    .edge(copies[i][copy(j)].expect("copy exists"))
                    .curve
                    .clone();
                let surface = Surface::Revolution {
                    frame,
                    profile: profile_curve,
                };
                let du = angles[j + 1] - angles[j];
                let normal = surface.normal(0.5 * du, tm).unwrap_or(DVec3::ZERO);
                let outward = rotate_vec(outward0, f0.z, angles[j] + 0.5 * du);
                let forward = normal.dot(outward) > 0.0;
                let mut ce = Vec::new();
                if let Some(a) = arcs[i][j] {
                    ce.push(fw(a));
                }
                ce.push(fw(copies[i][copy(j + 1)].expect("copy exists")));
                if let Some(a) = arcs[(i + 1) % n][j] {
                    ce.push(bw(a));
                }
                ce.push(bw(copies[i][copy(j)].expect("copy exists")));
                faces.push(brep.add_face(surface, forward, vec![loop_from(ce, !forward)]));
            }
        }
        if !full {
            start_cap.push(
                (0..n)
                    .map(|i| fw(copies[i][0].expect("copy exists")))
                    .collect(),
            );
            end_cap.push(
                (0..n)
                    .map(|i| fw(copies[i][k].expect("copy exists")))
                    .collect(),
            );
        }
    }
    if !full {
        // Start cap: the profile plane; outward is against the rotation (-y0).
        let start_out_is_z = (-f0.y).dot(plane.z) > 0.0;
        faces.push(
            brep.add_face(
                Surface::Plane(*plane),
                start_out_is_z,
                start_cap
                    .into_iter()
                    .map(|l| loop_from(l, !start_out_is_z))
                    .collect(),
            ),
        );
        let rk = rot(k);
        let end_plane = Frame {
            origin: rk.transform_point3(plane.origin),
            x: rk.transform_vector3(plane.x),
            y: rk.transform_vector3(plane.y),
            z: rk.transform_vector3(plane.z),
        };
        let yk = rk.transform_vector3(f0.y);
        let end_out_is_z = yk.dot(end_plane.z) > 0.0;
        faces.push(
            brep.add_face(
                Surface::Plane(end_plane),
                end_out_is_z,
                end_cap
                    .into_iter()
                    .map(|l| loop_from(l, !end_out_is_z))
                    .collect(),
            ),
        );
    }
    let shell = brep.add_shell(faces);
    brep.add_solid(vec![shell]);
    Ok(())
}
