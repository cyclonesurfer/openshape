//! The OpenDatum geometry kernel: a pure-Rust boundary-representation
//! modeller.
//!
//! The rest of the application talks to the kernel only through the
//! [`Kernel`] trait and the neutral data types in this module ([`Profile`],
//! [`BodyMesh`], [`MassProperties`], …). Bodies are opaque handles; their
//! topology and geometry are private so the kernel can evolve (or be
//! replaced) without touching features, rendering or I/O. See
//! `docs/adr/0001-geometry-kernel.md`.
//!
//! Units: lengths are metres, angles radians, all arithmetic `f64`.

pub mod arrangement;
mod boolean;
mod bvh;
mod check;
mod geom;
mod mass;
mod math;
mod origin;
mod profile;
mod sweep;
mod tess;
mod topo;

use std::sync::Arc;

pub use arrangement::{
    ArrangementEdge, ArrangementVertex, CurvatureSample, CurveContact, ProfileEdge, Region,
    RegionLoop, SketchArrangement,
};
pub use geom::SurfaceKind;
pub use glam::{DMat4, DVec2, DVec3};
pub use math::{Frame, Tolerance};
pub use origin::{
    Adjacency, EdgeIx, EdgeOrigin, FaceIx, FaceOrigin, Merged, Operand, Origins, PieceEnd,
    ProfilePiece, ProfileVertex, VertexIx, VertexOrigin,
};
pub use profile::{Profile, ProfileCurve};
pub use sweep::{ExtrudeSpec, RevolveSpec};
pub use tess::TessellationOptions;

/// Errors reported by kernel operations. An operation either returns a valid
/// body or one of these; it never returns a corrupt body.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum KernelError {
    /// The inputs are invalid (bad dimensions, open or self-intersecting
    /// profile, …).
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// The operation produced a body that failed validation. This is a kernel
    /// bug or an unsupported configuration; the result is discarded.
    #[error("operation produced an invalid body: {0}")]
    InvalidResult(String),
    /// Tessellation failed.
    #[error("tessellation failed: {0}")]
    Tessellation(String),
    /// Mass properties could not be computed to the requested accuracy.
    #[error("mass properties: {0}")]
    MassProperties(String),
    /// The inputs are valid but of a kind the operation does not handle
    /// yet (ADR 0020).
    #[error("not supported yet: {0}")]
    Unsupported(String),
    /// The operation would pass one of the kernel's work limits
    /// (docs/PERFORMANCE.md), which keep any input from running for minutes
    /// or exhausting memory. Nothing is produced.
    #[error("too large: {0}")]
    TooLarge(String),
}

/// A boolean operation (design note 0012).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum BooleanOp {
    /// Everything in either body.
    Union,
    /// What is in both bodies.
    Intersect,
    /// What is in the first body and not in the second.
    Subtract,
}

/// Entity counts of a body, plus the genus derived from Euler–Poincaré.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TopologyCounts {
    /// Solids (lumps).
    pub solids: usize,
    /// Shells.
    pub shells: usize,
    /// Faces.
    pub faces: usize,
    /// Loops (outer and inner).
    pub loops: usize,
    /// Edges.
    pub edges: usize,
    /// Vertices.
    pub vertices: usize,
    /// Total number of through-holes, if the Euler characteristic is consistent.
    pub genus: Option<usize>,
}

/// Mass properties for unit density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    /// Volume (m³).
    pub volume: f64,
    /// Surface area (m²).
    pub area: f64,
    /// Centre of mass.
    pub centroid: DVec3,
    /// Inertia tensor about the centroid (column-major rows), unit density.
    pub inertia: [[f64; 3]; 3],
}

/// Triangle mesh of one face.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FaceMesh {
    /// Vertex positions.
    pub positions: Vec<DVec3>,
    /// Unit outward normals per vertex.
    pub normals: Vec<DVec3>,
    /// Counter-clockwise (outward) triangles.
    pub triangles: Vec<[u32; 3]>,
}

/// The plane a planar face lies in, see [`Kernel::face_plane`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FacePlane {
    /// A point of the plane (not necessarily of the face).
    pub point: DVec3,
    /// The face's outward unit normal.
    pub normal: DVec3,
}

/// Display/export mesh of a body: one mesh per face plus edge polylines.
/// Face and edge order is stable for a given body and is the basis of
/// picking.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BodyMesh {
    /// Per-face meshes.
    pub faces: Vec<FaceMesh>,
    /// Per-edge polylines.
    pub edges: Vec<Vec<DVec3>>,
    /// Per edge: true for a seam, where two pieces of one surface meet (the
    /// split of a periodic face, ADR 0008, or the faces a sweep makes of
    /// two pieces of one curve). A seam bounds faces but is no crease in
    /// the shape; displays hide it.
    pub seams: Vec<bool>,
    /// Per vertex, its position (in the vertex order of
    /// [`Kernel::origins`]), for picking vertices.
    pub vertices: Vec<DVec3>,
}

impl BodyMesh {
    /// Total triangle count.
    pub fn triangle_count(&self) -> usize {
        self.faces.iter().map(|f| f.triangles.len()).sum()
    }

    /// Axis-aligned bounds of all vertices.
    pub fn bounds(&self) -> Option<(DVec3, DVec3)> {
        let mut it = self.faces.iter().flat_map(|f| f.positions.iter());
        let first = *it.next()?;
        Some(it.fold((first, first), |(lo, hi), p| (lo.min(*p), hi.max(*p))))
    }

    /// Merges the per-face meshes into one indexed triangle soup, welding
    /// vertices that are bit-identical (shared boundary samples). Returns
    /// positions and triangles; useful for export and watertightness checks.
    pub fn welded(&self) -> (Vec<DVec3>, Vec<[u32; 3]>) {
        let mut map = rustc_hash::FxHashMap::default();
        let mut pos = Vec::new();
        let mut tris = Vec::new();
        for f in &self.faces {
            let idx: Vec<u32> = f
                .positions
                .iter()
                .map(|p| {
                    let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
                    *map.entry(key).or_insert_with(|| {
                        pos.push(*p);
                        pos.len() as u32 - 1
                    })
                })
                .collect();
            for t in &f.triangles {
                tris.push(t.map(|i| idx[i as usize]));
            }
        }
        (pos, tris)
    }
}

/// An opaque solid body produced by a kernel, with the origins of its
/// entities (ADR 0018).
#[derive(Clone, Debug)]
pub struct Body {
    brep: Arc<topo::Brep>,
    origins: Arc<Origins>,
}

impl PartialEq for Body {
    fn eq(&self, other: &Self) -> bool {
        (Arc::ptr_eq(&self.brep, &other.brep) || self.brep == other.brep)
            && (Arc::ptr_eq(&self.origins, &other.origins) || self.origins == other.origins)
    }
}

/// The interface between the application and a geometry kernel.
///
/// Every operation validates its result; an invalid result is reported as
/// [`KernelError::InvalidResult`] rather than returned.
pub trait Kernel: Send + Sync {
    /// The kernel's opaque body handle.
    type Body: Clone + std::fmt::Debug + Send + Sync + PartialEq;

    /// The model tolerances in effect.
    fn tolerance(&self) -> Tolerance;

    /// An axis-aligned box with one corner at `frame.origin`, extending
    /// `size` along the frame axes. Negative sizes extend the other way.
    fn make_box(&self, frame: &Frame, size: DVec3) -> Result<Self::Body, KernelError>;

    /// Extrudes a planar profile along its plane normal.
    fn extrude(&self, profile: &Profile, spec: &ExtrudeSpec) -> Result<Self::Body, KernelError>;

    /// Revolves a planar profile about an axis lying in its plane.
    fn revolve(&self, profile: &Profile, spec: &RevolveSpec) -> Result<Self::Body, KernelError>;

    /// The union, intersection or difference (`a − b`) of two bodies
    /// (design note 0012, ADR 0020). The result's origins refer to the
    /// operands' faces, edges and vertices. Fails with
    /// [`KernelError::InvalidInput`] when the result would be empty or
    /// non-manifold (the bodies touching only along an edge or at a
    /// point), and with [`KernelError::Unsupported`] for what is not
    /// handled yet: faces other than planes and circular cylinders, faces
    /// that touch without crossing (tangent contacts), and cylinders that
    /// meet at an angle.
    fn boolean(
        &self,
        a: &Self::Body,
        b: &Self::Body,
        op: BooleanOp,
    ) -> Result<Self::Body, KernelError>;

    /// Tessellates a body for display or export.
    fn tessellate(
        &self,
        body: &Self::Body,
        opts: &TessellationOptions,
    ) -> Result<BodyMesh, KernelError>;

    /// Checks a body; returns the list of problems (empty when valid).
    fn validate(&self, body: &Self::Body) -> Vec<String>;

    /// Entity counts.
    fn topology(&self, body: &Self::Body) -> TopologyCounts;

    /// Mass properties by exact boundary integration (design note 0003).
    /// Each integral is computed to within `relative_accuracy` of the
    /// boundary integral of its integrand's magnitude (for the volume of a
    /// solid that is star-shaped about its vertex centroid, the volume
    /// itself), but no finer than rounding allows: about 1e-13 of the
    /// body's size to the integral's dimension, or, where the body's
    /// points are computed from larger coordinates (far from the origin,
    /// or on an arc of large radius), 1e-13 of those coordinates' size
    /// times the body's size to one dimension less. Fails rather than
    /// return a less accurate value.
    fn mass_properties(
        &self,
        body: &Self::Body,
        relative_accuracy: f64,
    ) -> Result<MassProperties, KernelError>;

    /// Analytic kinds of the body's faces, in face order.
    fn face_kinds(&self, body: &Self::Body) -> Vec<SurfaceKind>;

    /// How each face, edge and vertex of the body arose from the inputs of
    /// the operation that made it, in handle order (ADR 0018): the raw
    /// material of persistent names (design note 0013).
    fn origins<'a>(&self, body: &'a Self::Body) -> &'a Origins;

    /// Which faces each edge bounds, which edges end at each vertex, and
    /// which edges are seams.
    fn adjacency(&self, body: &Self::Body) -> Adjacency;

    /// The plane of a face whose surface is a plane (see
    /// [`Kernel::face_kinds`]), with the face's outward normal; `None` for
    /// a face of any other kind or an index past the body's faces.
    fn face_plane(&self, body: &Self::Body, face: FaceIx) -> Option<FacePlane>;

    /// The tangent plane of a face's surface at the point of the surface
    /// nearest `near` (exactly for planes, cylinders, cones and spheres,
    /// closely for other surfaces and points near them), with the face's
    /// outward normal there; `None` for an index past the body's faces or
    /// a singular point of the surface (a cone's apex).
    fn tangent_plane(&self, body: &Self::Body, face: FaceIx, near: DVec3) -> Option<FacePlane>;

    /// The planar arrangement of sketch curves: vertices, edges and bounded
    /// regions (design note 0006).
    fn sketch_regions(&self, curves: &[ProfileCurve]) -> Result<SketchArrangement, KernelError>;

    /// The boundary of the union of the given regions of `arrangement` as
    /// profile curves, ready for [`Kernel::extrude`] or [`Kernel::revolve`].
    fn region_profile(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        regions: &[usize],
    ) -> Result<Vec<ProfileCurve>, KernelError>;

    /// [`Kernel::region_profile`] with, beside each curve, the arrangement
    /// edge it is a piece of and the direction the profile runs along it
    /// (an edge cut in two gives two consecutive curves of one edge): how
    /// features name the faces swept from the profile.
    fn region_profile_edges(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        regions: &[usize],
    ) -> Result<Vec<(ProfileCurve, ProfileEdge)>, KernelError>;

    /// The region of `arrangement` containing `point` (at least the model
    /// tolerance inside it), if any: for picking regions and for re-finding
    /// a named region after a sketch edit.
    fn region_containing(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        point: DVec2,
    ) -> Option<usize>;

    /// Where two profile curves meet: crossings, tangencies, end points on
    /// the other curve and the ends of overlaps, as parameter pairs in the
    /// arrangement's conventions, sorted along `a` (design note 0007 §6).
    /// For sketch tools such as trim and extend.
    fn curve_contacts(
        &self,
        a: &ProfileCurve,
        b: &ProfileCurve,
    ) -> Result<Vec<CurveContact>, KernelError>;

    /// The point of a profile curve at a parameter in the arrangement's
    /// conventions (clamped to the curve's domain).
    fn curve_point(&self, curve: &ProfileCurve, t: f64) -> Result<DVec2, KernelError>;

    /// Points, unit tangents and signed curvatures along a profile curve,
    /// `n + 1` evenly spaced by length plus both sides of any curvature jump
    /// at a spline knot: for curvature combs.
    fn curvature_samples(
        &self,
        curve: &ProfileCurve,
        n: usize,
    ) -> Result<Vec<CurvatureSample>, KernelError>;
}

/// The built-in OpenDatum kernel.
#[derive(Clone, Debug, Default)]
pub struct BrepKernel {
    tol: Tolerance,
}

impl BrepKernel {
    /// A kernel with the given tolerances.
    pub fn with_tolerance(tol: Tolerance) -> Self {
        Self { tol }
    }

    fn finish(&self, (brep, origins): (topo::Brep, Origins)) -> Result<Body, KernelError> {
        let mut issues = check::validate(&brep, &self.tol);
        let counts = (brep.faces.len(), brep.edges.len(), brep.vertices.len());
        let table = (
            origins.faces.len(),
            origins.edges.len(),
            origins.vertices.len(),
        );
        if counts != table {
            issues.push(format!(
                "origin table covers {table:?} faces, edges and vertices of {counts:?}"
            ));
        }
        if issues.is_empty() {
            Ok(Body {
                brep: Arc::new(brep),
                origins: Arc::new(origins),
            })
        } else {
            Err(KernelError::InvalidResult(issues.join("; ")))
        }
    }
}

impl Kernel for BrepKernel {
    type Body = Body;

    fn tolerance(&self) -> Tolerance {
        self.tol
    }

    fn make_box(&self, frame: &Frame, size: DVec3) -> Result<Body, KernelError> {
        if !size.is_finite() || size.abs().min_element() <= self.tol.linear {
            return Err(KernelError::InvalidInput(
                "box dimensions must be non-zero".into(),
            ));
        }
        let profile = Profile::rectangle(
            *frame,
            DVec2::ZERO.min(size.truncate()),
            DVec2::ZERO.max(size.truncate()),
        );
        self.extrude(
            &profile,
            &ExtrudeSpec {
                distance: size.z,
                symmetric: false,
            },
        )
    }

    fn extrude(&self, profile: &Profile, spec: &ExtrudeSpec) -> Result<Body, KernelError> {
        self.finish(sweep::extrude(profile, spec, &self.tol)?)
    }

    fn revolve(&self, profile: &Profile, spec: &RevolveSpec) -> Result<Body, KernelError> {
        self.finish(sweep::revolve(profile, spec, &self.tol)?)
    }

    fn boolean(&self, a: &Body, b: &Body, op: BooleanOp) -> Result<Body, KernelError> {
        self.finish(boolean::boolean(&a.brep, &b.brep, op, &self.tol)?)
    }

    fn tessellate(&self, body: &Body, opts: &TessellationOptions) -> Result<BodyMesh, KernelError> {
        let t = tess::tessellate(&body.brep, opts, &self.tol)?;
        Ok(BodyMesh {
            faces: t
                .faces
                .into_iter()
                .map(|f| FaceMesh {
                    positions: f.positions,
                    normals: f.normals,
                    triangles: f.triangles,
                })
                .collect(),
            edges: t.edges,
            seams: t.seams,
            vertices: body.brep.vertices.iter().map(|v| v.point).collect(),
        })
    }

    fn validate(&self, body: &Body) -> Vec<String> {
        check::validate(&body.brep, &self.tol)
    }

    fn topology(&self, body: &Body) -> TopologyCounts {
        check::counts(&body.brep)
    }

    fn mass_properties(
        &self,
        body: &Body,
        relative_accuracy: f64,
    ) -> Result<MassProperties, KernelError> {
        if !relative_accuracy.is_finite() || relative_accuracy <= 0.0 {
            return Err(KernelError::InvalidInput(
                "accuracy must be positive".into(),
            ));
        }
        let size = check::bounding_size(&body.brep).max(self.tol.linear);
        mass::mass_properties(&body.brep, size, relative_accuracy)
    }

    fn face_kinds(&self, body: &Body) -> Vec<SurfaceKind> {
        body.brep.faces.iter().map(|f| f.surface.kind()).collect()
    }

    fn origins<'a>(&self, body: &'a Body) -> &'a Origins {
        &body.origins
    }

    fn adjacency(&self, body: &Body) -> Adjacency {
        let b = &body.brep;
        let mut edge_faces = vec![Vec::new(); b.edges.len()];
        for (fi, f) in b.faces.iter().enumerate() {
            for c in f.loops.iter().flat_map(|l| &l.coedges) {
                edge_faces[c.edge.idx()].push(FaceIx(fi as u32));
            }
        }
        let mut vertex_edges = vec![Vec::new(); b.vertices.len()];
        for (ei, e) in b.edges.iter().enumerate() {
            vertex_edges[e.v0.idx()].push(EdgeIx(ei as u32));
            if e.v1 != e.v0 {
                vertex_edges[e.v1.idx()].push(EdgeIx(ei as u32));
            }
        }
        Adjacency {
            edge_faces,
            vertex_edges,
            seams: tess::seams(b),
        }
    }

    fn tangent_plane(&self, body: &Body, face: FaceIx, near: DVec3) -> Option<FacePlane> {
        let f = body.brep.faces.get(face.0 as usize)?;
        // A plane is regular everywhere, even where its parameters are not
        // (the centre of a flat cap of a revolve).
        if let Some((o, n)) = f.surface.plane() {
            return Some(FacePlane {
                point: near - n * (near - o).dot(n),
                normal: if f.forward { n } else { -n },
            });
        }
        let (u, v) = f.surface.project(near);
        let n = f.surface.normal_or_limit(u, v)?;
        Some(FacePlane {
            point: f.surface.point(u, v),
            normal: if f.forward { n } else { -n },
        })
    }

    fn face_plane(&self, body: &Body, face: FaceIx) -> Option<FacePlane> {
        let f = body.brep.faces.get(face.0 as usize)?;
        let (point, n) = f.surface.plane()?;
        let normal = if f.forward { n } else { -n };
        Some(FacePlane { point, normal })
    }

    fn sketch_regions(&self, curves: &[ProfileCurve]) -> Result<SketchArrangement, KernelError> {
        arrangement::arrange(curves, &self.tol)
    }

    fn region_profile(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        regions: &[usize],
    ) -> Result<Vec<ProfileCurve>, KernelError> {
        arrangement::region_profile_with_tolerance(arrangement, curves, regions, &self.tol)
    }

    fn region_profile_edges(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        regions: &[usize],
    ) -> Result<Vec<(ProfileCurve, ProfileEdge)>, KernelError> {
        arrangement::region_profile_edges_with_tolerance(arrangement, curves, regions, &self.tol)
    }

    fn region_containing(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        point: DVec2,
    ) -> Option<usize> {
        arrangement::region_containing(arrangement, curves, point, &self.tol)
    }

    fn curve_contacts(
        &self,
        a: &ProfileCurve,
        b: &ProfileCurve,
    ) -> Result<Vec<CurveContact>, KernelError> {
        arrangement::curve_contacts(a, b, &self.tol)
    }

    fn curve_point(&self, curve: &ProfileCurve, t: f64) -> Result<DVec2, KernelError> {
        arrangement::curve_point(curve, t, &self.tol)
    }

    fn curvature_samples(
        &self,
        curve: &ProfileCurve,
        n: usize,
    ) -> Result<Vec<CurvatureSample>, KernelError> {
        arrangement::curvature_samples(curve, n, &self.tol)
    }
}
