//! The OpenShape application shell.
//!
//! [`Model`] is the UI-independent controller: an editing [`Session`] (the op
//! log and its derived document), the geometry kernel, and the regenerated
//! parts with their display meshes. The egui front end lives in [`app`];
//! [`demo`] drives the same model offscreen to prove the pipeline end to end.

pub mod app;
pub mod demo;
pub mod dialogs;

use std::path::Path;

use openshape_core_model::{
    ElementId, Feature, FeatureId, InsertAt, ModelOp, PartStudio, ReplicaId, Session,
};
use openshape_features::{regenerate, FeatureStatus, StudioResult};
use openshape_kernel::{Body, BodyMesh, Kernel, OsKernel, TessellationOptions};
use openshape_render::MeshData;

/// Colours cycled through for successive parts (sRGB).
const PART_COLORS: [[f32; 4]; 6] = [
    [0.62, 0.70, 0.80, 1.0],
    [0.80, 0.66, 0.52, 1.0],
    [0.60, 0.78, 0.62, 1.0],
    [0.78, 0.62, 0.76, 1.0],
    [0.82, 0.78, 0.55, 1.0],
    [0.58, 0.76, 0.78, 1.0],
];

/// Converts a kernel mesh into renderer mesh data.
pub fn to_mesh_data(mesh: &BodyMesh, color: [f32; 4]) -> MeshData {
    let mut out = MeshData {
        color,
        ..MeshData::default()
    };
    for f in &mesh.faces {
        let base = out.positions.len() as u32;
        out.positions
            .extend(f.positions.iter().map(|p| p.as_vec3().to_array()));
        out.normals
            .extend(f.normals.iter().map(|n| n.as_vec3().to_array()));
        for t in &f.triangles {
            out.indices.extend(t.iter().map(|i| base + i));
        }
    }
    for e in &mesh.edges {
        for w in e.windows(2) {
            out.edge_lines.push(w[0].as_vec3().to_array());
            out.edge_lines.push(w[1].as_vec3().to_array());
        }
    }
    out
}

/// The application model: document session, kernel, and regenerated parts.
pub struct Model {
    session: Session,
    kernel: OsKernel,
    result: StudioResult<Body>,
    meshes: Vec<MeshData>,
    tessellation: TessellationOptions,
    /// Incremented whenever the display meshes change.
    pub generation: u64,
}

impl Default for Model {
    fn default() -> Self {
        Self::new()
    }
}

impl Model {
    /// A new, empty document.
    pub fn new() -> Self {
        Self::from_session(Session::new_document(ReplicaId::generate()))
    }

    fn from_session(session: Session) -> Self {
        let mut m = Model {
            session,
            kernel: OsKernel::default(),
            result: StudioResult {
                statuses: Vec::new(),
                bodies: Vec::new(),
            },
            meshes: Vec::new(),
            tessellation: TessellationOptions {
                chord_tolerance: 2e-5,
                max_angle: 12f64.to_radians(),
            },
            generation: 0,
        };
        m.regenerate();
        m
    }

    /// Opens a document saved with [`Model::save`].
    pub fn open(path: &Path) -> Result<Self, String> {
        let bytes = std::fs::read(path).map_err(|e| e.to_string())?;
        let session = Session::open(&bytes, ReplicaId::generate()).map_err(|e| e.to_string())?;
        Ok(Self::from_session(session))
    }

    /// Saves the document (the op log) in the native format.
    pub fn save(&self, path: &Path) -> Result<(), String> {
        let bytes = self.session.save().map_err(|e| e.to_string())?;
        std::fs::write(path, bytes).map_err(|e| e.to_string())
    }

    /// Exports all parts as binary STL in millimetres.
    pub fn export_stl(&self, path: &Path) -> Result<(), String> {
        let mut meshes = Vec::new();
        for b in &self.result.bodies {
            meshes.push(
                self.kernel
                    .tessellate(&b.body, &self.tessellation)
                    .map_err(|e| e.to_string())?,
            );
        }
        let mut f =
            std::io::BufWriter::new(std::fs::File::create(path).map_err(|e| e.to_string())?);
        openshape_io::stl::write_binary(&mut f, &meshes, openshape_io::stl::MILLIMETRES)
            .map_err(|e| e.to_string())
    }

    /// The Part Studio being edited.
    pub fn studio_id(&self) -> ElementId {
        self.session
            .first_studio()
            .expect("a document always has its genesis Part Studio")
    }

    /// The Part Studio's feature list.
    pub fn studio(&self) -> &PartStudio {
        &self
            .session
            .document()
            .element(self.studio_id())
            .expect("studio exists")
            .studio
    }

    /// The editing session (op log and document).
    pub fn session(&self) -> &Session {
        &self.session
    }

    /// The kernel.
    pub fn kernel(&self) -> &OsKernel {
        &self.kernel
    }

    /// Regeneration status of a feature.
    pub fn status(&self, id: FeatureId) -> Option<&FeatureStatus> {
        self.result.status(id)
    }

    /// Regenerated bodies.
    pub fn bodies(&self) -> &StudioResult<Body> {
        &self.result
    }

    /// Display meshes, one per part.
    pub fn meshes(&self) -> &[MeshData] {
        &self.meshes
    }

    /// Bounds of all display meshes.
    pub fn bounds(&self) -> Option<([f32; 3], [f32; 3])> {
        self.meshes
            .iter()
            .filter_map(|m| m.bounds())
            .reduce(|(a, b), (c, d)| {
                (
                    [a[0].min(c[0]), a[1].min(c[1]), a[2].min(c[2])],
                    [b[0].max(d[0]), b[1].max(d[1]), b[2].max(d[2])],
                )
            })
    }

    /// Applies an edit through the op log and regenerates.
    pub fn apply(&mut self, op: ModelOp) {
        self.session.apply(op);
        self.regenerate();
    }

    /// Appends a feature at the end of the Part Studio.
    pub fn add_feature(&mut self, name: &str, feature: Feature) -> FeatureId {
        let element = self.studio_id();
        let id = self.session.apply(ModelOp::InsertFeature {
            element,
            at: InsertAt::End,
            name: name.to_string(),
            feature,
        });
        self.regenerate();
        FeatureId(id)
    }

    /// Undo; returns false when there is nothing to undo.
    pub fn undo(&mut self) -> bool {
        let r = self.session.undo();
        if r {
            self.regenerate();
        }
        r
    }

    /// Redo; returns false when there is nothing to redo.
    pub fn redo(&mut self) -> bool {
        let r = self.session.redo();
        if r {
            self.regenerate();
        }
        r
    }

    /// A fresh name like "Extrude 3" for a new feature of this type.
    pub fn next_name(&self, type_name: &str) -> String {
        let n = self
            .studio()
            .features
            .iter()
            .filter(|f| f.feature.type_name() == type_name)
            .count();
        format!("{type_name} {}", n + 1)
    }

    /// Regenerates all parts and rebuilds display meshes.
    pub fn regenerate(&mut self) {
        let studio = self.studio().clone();
        self.result = regenerate(&self.kernel, &studio);
        let mut meshes = Vec::new();
        for (i, b) in self.result.bodies.iter().enumerate() {
            match self.kernel.tessellate(&b.body, &self.tessellation) {
                Ok(m) => meshes.push(to_mesh_data(&m, PART_COLORS[i % PART_COLORS.len()])),
                Err(e) => {
                    if let Some(s) = self
                        .result
                        .statuses
                        .iter_mut()
                        .find(|(f, _)| *f == b.feature)
                    {
                        s.1 = FeatureStatus::Failed(format!("display: {e}"));
                    }
                }
            }
        }
        self.meshes = meshes;
        self.generation += 1;
    }
}
