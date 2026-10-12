//! The OpenShape geometry kernel: a pure-Rust boundary-representation
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
mod check;
mod geom;
mod mass;
mod math;
mod profile;
mod sweep;
mod tess;
mod topo;

use std::sync::Arc;

pub use arrangement::{ArrangementEdge, ArrangementVertex, Region, RegionLoop, SketchArrangement};
pub use geom::SurfaceKind;
pub use glam::{DMat4, DVec2, DVec3};
pub use math::{Frame, Tolerance};
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
    /// split of a periodic face, ADR 0008). A seam bounds faces but is no
    /// crease in the shape; displays hide it.
    pub seams: Vec<bool>,
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
        let mut map = std::collections::HashMap::new();
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

/// An opaque solid body produced by a kernel.
#[derive(Clone, Debug)]
pub struct Body {
    brep: Arc<topo::Brep>,
}

impl PartialEq for Body {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.brep, &other.brep) || self.brep == other.brep
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

    /// Mass properties by exact boundary integration (design note 0003),
    /// to the given accuracy relative to the body's size.
    fn mass_properties(
        &self,
        body: &Self::Body,
        relative_accuracy: f64,
    ) -> Result<MassProperties, KernelError>;

    /// Analytic kinds of the body's faces, in face order.
    fn face_kinds(&self, body: &Self::Body) -> Vec<SurfaceKind>;

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

    /// The region of `arrangement` containing `point` (at least the model
    /// tolerance inside it), if any: for picking regions and for re-finding
    /// a named region after a sketch edit.
    fn region_containing(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        point: DVec2,
    ) -> Option<usize>;
}

/// The built-in OpenShape kernel.
#[derive(Clone, Debug, Default)]
pub struct OsKernel {
    tol: Tolerance,
}

impl OsKernel {
    /// A kernel with the given tolerances.
    pub fn with_tolerance(tol: Tolerance) -> Self {
        Self { tol }
    }

    fn finish(&self, brep: topo::Brep) -> Result<Body, KernelError> {
        let issues = check::validate(&brep, &self.tol);
        if issues.is_empty() {
            Ok(Body {
                brep: Arc::new(brep),
            })
        } else {
            Err(KernelError::InvalidResult(issues.join("; ")))
        }
    }
}

impl Kernel for OsKernel {
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

    fn tessellate(&self, body: &Body, opts: &TessellationOptions) -> Result<BodyMesh, KernelError> {
        let t = tess::tessellate(&body.brep, opts)?;
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
        Ok(mass::mass_properties(&body.brep, size * relative_accuracy))
    }

    fn face_kinds(&self, body: &Body) -> Vec<SurfaceKind> {
        body.brep.faces.iter().map(|f| f.surface.kind()).collect()
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

    fn region_containing(
        &self,
        arrangement: &SketchArrangement,
        curves: &[ProfileCurve],
        point: DVec2,
    ) -> Option<usize> {
        arrangement::region_containing(arrangement, curves, point, &self.tol)
    }
}
