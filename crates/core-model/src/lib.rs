//! The OpenDatum document model.
//!
//! A [`Document`] is derived state: it is computed by replaying the effective
//! [`ModelOp`]s of an [`OpLog`] through a deterministic reducer
//! ([`Document::apply`]). [`Session`] bundles a log with its cached document
//! and is what applications edit.
//!
//! The model is independent of the geometry kernel: features are described
//! with plain data (planes, profile curves, dimensions in SI units) and are
//! evaluated by the `features` crate.

pub mod sketch;
pub mod units;

pub use glam::{DVec2, DVec3};
pub use opendatum_oplog::{Op, OpId, OpKind, OpLog, OpLogError, ReplicaId};
use serde::{Deserialize, Serialize};
pub use sketch::{RegionKey, SketchData, SketchDef, SketchEdit, SketchId};

/// Identifies an element (tab) of a document; derived from the creating op.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ElementId(pub OpId);

/// Identifies a feature; derived from the creating op.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct FeatureId(pub OpId);

/// A sketch/profile plane.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum PlaneRef {
    /// World XY (normal +Z).
    Top,
    /// World XZ (normal −Y).
    Front,
    /// World YZ (normal +X).
    Right,
    /// An explicit plane.
    Custom {
        /// Origin.
        origin: DVec3,
        /// Normal.
        normal: DVec3,
        /// In-plane x direction hint.
        x_dir: DVec3,
    },
}

/// One curve of a feature profile, in plane coordinates (metres/radians).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum CurveDef {
    /// Straight segment.
    Line {
        /// Start.
        start: DVec2,
        /// End.
        end: DVec2,
    },
    /// Circular arc; positive sweep is counter-clockwise.
    Arc {
        /// Centre.
        center: DVec2,
        /// Radius.
        radius: f64,
        /// Start angle.
        start_angle: f64,
        /// Signed sweep.
        sweep_angle: f64,
    },
    /// Full circle.
    Circle {
        /// Centre.
        center: DVec2,
        /// Radius.
        radius: f64,
    },
    /// Interpolating spline through points.
    Spline {
        /// Points to pass through (first and last equal for a closed spline).
        through: Vec<DVec2>,
    },
}

/// A planar profile made of closed loops of curves.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProfileDef {
    /// Plane.
    pub plane: PlaneRef,
    /// Curves.
    pub curves: Vec<CurveDef>,
}

/// Where a sweep's profile comes from.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ProfileSource {
    /// Curves given directly.
    Inline(ProfileDef),
    /// Regions of an earlier sketch feature; the profile is the boundary of
    /// their union, in the sketch's plane.
    Sketch {
        /// The sketch feature.
        sketch: FeatureId,
        /// The selected regions.
        regions: Vec<RegionKey>,
    },
}

/// A revolve axis, in the profile's plane coordinates.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum RevolveAxis {
    /// A point and a direction.
    Inline {
        /// Point on the axis.
        origin: DVec2,
        /// Direction.
        dir: DVec2,
    },
    /// A line of the profile's sketch (typically construction geometry),
    /// from its start towards its end.
    SketchLine {
        /// The line entity.
        line: SketchId,
    },
}

/// A modeling feature: pure data, evaluated by the `features` crate.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Feature {
    /// A 2D sketch: entities and constraints on a plane.
    Sketch(SketchDef),
    /// An axis-aligned box primitive with one corner at `corner`.
    Box {
        /// Corner position.
        corner: DVec3,
        /// Extents along x, y, z (metres; may be negative).
        size: DVec3,
    },
    /// Linear extrusion of a profile.
    Extrude {
        /// Profile.
        profile: ProfileSource,
        /// Signed distance along the plane normal.
        distance: f64,
        /// Extrude half to each side.
        symmetric: bool,
    },
    /// Rotational sweep of a profile about an in-plane axis.
    Revolve {
        /// Profile.
        profile: ProfileSource,
        /// Axis.
        axis: RevolveAxis,
        /// Signed angle in radians.
        angle: f64,
    },
}

impl Feature {
    /// A short type name for display.
    pub fn type_name(&self) -> &'static str {
        match self {
            Feature::Sketch(_) => "Sketch",
            Feature::Box { .. } => "Box",
            Feature::Extrude { .. } => "Extrude",
            Feature::Revolve { .. } => "Revolve",
        }
    }
}

/// A feature in a part studio's feature list.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FeatureEntry {
    /// Stable id.
    pub id: FeatureId,
    /// Display name.
    pub name: String,
    /// Suppressed features are skipped by regeneration.
    pub suppressed: bool,
    /// Definition.
    pub feature: Feature,
}

/// A Part Studio: an ordered feature list producing parts.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PartStudio {
    /// Features in regeneration order.
    pub features: Vec<FeatureEntry>,
    /// The rollback bar (ADR 0013): features after it are not regenerated,
    /// and the application inserts new features at it. The reducer keeps it
    /// normalized: [`InsertAt::After`] names an existing feature other than
    /// the last, [`InsertAt::Start`] is used only when there are features
    /// after it, and a bar after every feature is [`InsertAt::End`].
    #[serde(default)]
    pub rollback: InsertAt,
}

impl PartStudio {
    /// Index of a feature.
    pub fn index_of(&self, id: FeatureId) -> Option<usize> {
        self.features.iter().position(|f| f.id == id)
    }

    /// The number of features before the rollback bar: `features[..i]` are
    /// regenerated, the rest are rolled back.
    pub fn rollback_index(&self) -> usize {
        match self.rollback {
            InsertAt::End => self.features.len(),
            InsertAt::Start => 0,
            InsertAt::After(f) => self.index_of(f).map_or(self.features.len(), |i| i + 1),
        }
    }

    /// The features before the rollback bar, in regeneration order.
    pub fn active_features(&self) -> &[FeatureEntry] {
        &self.features[..self.rollback_index()]
    }

    /// True if the feature exists and is after the rollback bar.
    pub fn is_rolled_back(&self, id: FeatureId) -> bool {
        self.index_of(id)
            .is_some_and(|i| i >= self.rollback_index())
    }

    /// The bar's canonical form for its current position (see
    /// [`PartStudio::rollback`]).
    fn normalize_rollback(&mut self) {
        let i = self.rollback_index();
        self.rollback = if i >= self.features.len() {
            InsertAt::End
        } else if i == 0 {
            InsertAt::Start
        } else {
            InsertAt::After(self.features[i - 1].id)
        };
    }

    /// Before the feature at `index` is removed: if it anchors the bar, the
    /// bar moves to the feature before it (or the start), so the same
    /// features stay active.
    fn release_anchor(&mut self, index: usize) {
        if self.rollback == InsertAt::After(self.features[index].id) {
            self.rollback = match index {
                0 => InsertAt::Start,
                i => InsertAt::After(self.features[i - 1].id),
            };
        }
    }

    /// Inserts `entry` at `pos`. A feature placed exactly at a rollback bar
    /// that is not at the end goes above it (the bar advances past it).
    fn insert_at_bar_rule(&mut self, pos: usize, entry: FeatureEntry) {
        let bar = self.rollback_index();
        let id = entry.id;
        self.features.insert(pos, entry);
        if self.rollback != InsertAt::End && pos == bar {
            self.rollback = InsertAt::After(id);
        }
        self.normalize_rollback();
    }
}

/// Kinds of elements a document can hold. More (assemblies, drawings,
/// variable studios, scripts, attachments) arrive in later phases.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElementKind {
    /// Part Studio.
    PartStudio,
}

/// An element (tab) of a document.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Element {
    /// Stable id.
    pub id: ElementId,
    /// Tab name.
    pub name: String,
    /// Part Studio content.
    pub studio: PartStudio,
}

/// Where to insert a feature; also the position of a rollback bar.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum InsertAt {
    /// At the end of the list.
    #[default]
    End,
    /// At the start of the list.
    Start,
    /// Directly after the given feature.
    After(FeatureId),
}

/// An edit to the document; the payload of the op log.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ModelOp {
    /// Creates an element; its id is the op id.
    CreateElement {
        /// Name.
        name: String,
        /// Kind.
        kind: ElementKind,
    },
    /// Deletes an element. Refused for a document's last element.
    DeleteElement {
        /// Element.
        element: ElementId,
    },
    /// Renames an element.
    RenameElement {
        /// Element.
        element: ElementId,
        /// New name.
        name: String,
    },
    /// Inserts a feature; its id is the op id.
    InsertFeature {
        /// Element.
        element: ElementId,
        /// Position.
        at: InsertAt,
        /// Display name.
        name: String,
        /// Definition.
        feature: Feature,
    },
    /// Replaces a feature's definition.
    EditFeature {
        /// Element.
        element: ElementId,
        /// Feature.
        feature: FeatureId,
        /// New definition.
        definition: Feature,
    },
    /// Renames a feature.
    RenameFeature {
        /// Element.
        element: ElementId,
        /// Feature.
        feature: FeatureId,
        /// New name.
        name: String,
    },
    /// Suppresses or unsuppresses a feature.
    SetSuppressed {
        /// Element.
        element: ElementId,
        /// Feature.
        feature: FeatureId,
        /// New state.
        suppressed: bool,
    },
    /// Moves a feature.
    MoveFeature {
        /// Element.
        element: ElementId,
        /// Feature.
        feature: FeatureId,
        /// New position.
        to: InsertAt,
    },
    /// Deletes a feature.
    DeleteFeature {
        /// Element.
        element: ElementId,
        /// Feature.
        feature: FeatureId,
    },
    /// Moves a Part Studio's rollback bar (ADR 0013).
    SetRollback {
        /// Element.
        element: ElementId,
        /// New position: the bar goes after this feature, or at the start
        /// or end of the list.
        at: InsertAt,
    },
    /// Edits a sketch feature: a batch applied atomically (ADR 0011).
    EditSketch {
        /// Element.
        element: ElementId,
        /// The sketch feature.
        feature: FeatureId,
        /// Edits, in order.
        edits: Vec<SketchEdit>,
    },
}

/// A replayed op that could not be applied (e.g. an edit to a feature that
/// a concurrent op deleted). Recorded, never fatal.
#[derive(Clone, Debug, PartialEq)]
pub struct Conflict {
    /// The op.
    pub op: OpId,
    /// Why it was skipped.
    pub reason: String,
}

/// The derived document state.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Document {
    /// Elements in creation order.
    pub elements: Vec<Element>,
    /// Ops skipped during replay.
    pub conflicts: Vec<Conflict>,
}

impl Document {
    /// Replays a log's effective ops.
    pub fn replay(log: &OpLog<ModelOp>) -> Document {
        Self::from_ops(log.effective())
    }

    /// Replays the log as of the given heads (history/versions).
    pub fn replay_at(log: &OpLog<ModelOp>, heads: &[OpId]) -> Document {
        Self::from_ops(log.effective_at(heads))
    }

    fn from_ops(ops: Vec<(OpId, &ModelOp)>) -> Document {
        let mut d = Document::default();
        for (id, op) in ops {
            if let Err(reason) = d.apply(id, op) {
                d.conflicts.push(Conflict { op: id, reason });
            }
        }
        d
    }

    /// Looks up an element.
    pub fn element(&self, id: ElementId) -> Option<&Element> {
        self.elements.iter().find(|e| e.id == id)
    }

    fn element_mut(&mut self, id: ElementId) -> Result<&mut Element, String> {
        self.elements
            .iter_mut()
            .find(|e| e.id == id)
            .ok_or_else(|| "element no longer exists".to_string())
    }

    /// The deterministic reducer: applies one op.
    pub fn apply(&mut self, id: OpId, op: &ModelOp) -> Result<(), String> {
        match op {
            ModelOp::CreateElement { name, kind } => {
                match kind {
                    ElementKind::PartStudio => self.elements.push(Element {
                        id: ElementId(id),
                        name: name.clone(),
                        studio: PartStudio::default(),
                    }),
                }
                Ok(())
            }
            ModelOp::DeleteElement { element } => {
                let i = self
                    .elements
                    .iter()
                    .position(|e| e.id == *element)
                    .ok_or("element no longer exists")?;
                if self.elements.len() == 1 {
                    return Err("a document keeps at least one element".into());
                }
                self.elements.remove(i);
                Ok(())
            }
            ModelOp::RenameElement { element, name } => {
                self.element_mut(*element)?.name = name.clone();
                Ok(())
            }
            ModelOp::InsertFeature {
                element,
                at,
                name,
                feature,
            } => {
                let studio = &mut self.element_mut(*element)?.studio;
                let entry = FeatureEntry {
                    id: FeatureId(id),
                    name: name.clone(),
                    suppressed: false,
                    feature: feature.clone(),
                };
                let pos = position(studio, *at)?;
                studio.insert_at_bar_rule(pos, entry);
                Ok(())
            }
            ModelOp::EditFeature {
                element,
                feature,
                definition,
            } => {
                feature_mut(&mut self.element_mut(*element)?.studio, *feature)?.feature =
                    definition.clone();
                Ok(())
            }
            ModelOp::RenameFeature {
                element,
                feature,
                name,
            } => {
                feature_mut(&mut self.element_mut(*element)?.studio, *feature)?.name = name.clone();
                Ok(())
            }
            ModelOp::SetSuppressed {
                element,
                feature,
                suppressed,
            } => {
                feature_mut(&mut self.element_mut(*element)?.studio, *feature)?.suppressed =
                    *suppressed;
                Ok(())
            }
            ModelOp::MoveFeature {
                element,
                feature,
                to,
            } => {
                let studio = &mut self.element_mut(*element)?.studio;
                if *to == InsertAt::After(*feature) {
                    return Ok(());
                }
                let from = studio
                    .index_of(*feature)
                    .ok_or("feature no longer exists")?;
                if let InsertAt::After(anchor) = to {
                    studio
                        .index_of(*anchor)
                        .ok_or("anchor feature no longer exists")?;
                }
                studio.release_anchor(from);
                let entry = studio.features.remove(from);
                studio.normalize_rollback();
                let p = position(studio, *to).expect("the anchor was checked");
                studio.insert_at_bar_rule(p, entry);
                Ok(())
            }
            ModelOp::DeleteFeature { element, feature } => {
                let studio = &mut self.element_mut(*element)?.studio;
                let i = studio
                    .index_of(*feature)
                    .ok_or("feature no longer exists")?;
                studio.release_anchor(i);
                studio.features.remove(i);
                studio.normalize_rollback();
                Ok(())
            }
            ModelOp::SetRollback { element, at } => {
                let studio = &mut self.element_mut(*element)?.studio;
                if let InsertAt::After(f) = at {
                    studio.index_of(*f).ok_or("feature no longer exists")?;
                }
                studio.rollback = *at;
                studio.normalize_rollback();
                Ok(())
            }
            ModelOp::EditSketch {
                element,
                feature,
                edits,
            } => match &mut feature_mut(&mut self.element_mut(*element)?.studio, *feature)?.feature
            {
                Feature::Sketch(def) => def.apply_edits(id, edits),
                _ => Err("feature is not a sketch".into()),
            },
        }
    }
}

fn position(studio: &PartStudio, at: InsertAt) -> Result<usize, String> {
    match at {
        InsertAt::End => Ok(studio.features.len()),
        InsertAt::Start => Ok(0),
        InsertAt::After(f) => studio
            .index_of(f)
            .map(|i| i + 1)
            .ok_or_else(|| "anchor feature no longer exists".to_string()),
    }
}

fn feature_mut(studio: &mut PartStudio, id: FeatureId) -> Result<&mut FeatureEntry, String> {
    studio
        .features
        .iter_mut()
        .find(|f| f.id == id)
        .ok_or_else(|| "feature no longer exists".to_string())
}

/// The replica that writes genesis ops (a new document's first Part Studio).
/// Its ops never appear on any editing session's undo stack.
pub const GENESIS: ReplicaId = ReplicaId(0);

/// An editing session: an op log plus its derived document.
#[derive(Clone, Debug)]
pub struct Session {
    log: OpLog<ModelOp>,
    doc: Document,
}

impl Session {
    /// A new document with one empty Part Studio, edited by `replica`.
    pub fn new_document(replica: ReplicaId) -> Session {
        let mut log = OpLog::new(replica);
        let genesis = Op {
            id: OpId {
                lamport: 1,
                replica: GENESIS,
            },
            deps: vec![],
            kind: OpKind::Apply(ModelOp::CreateElement {
                name: "Part Studio 1".into(),
                kind: ElementKind::PartStudio,
            }),
        };
        log.integrate(genesis).expect("genesis op is valid");
        let doc = Document::replay(&log);
        Session { log, doc }
    }

    /// Opens a serialized log.
    pub fn open(bytes: &[u8], replica: ReplicaId) -> Result<Session, OpLogError> {
        let log = OpLog::decode(bytes, replica)?;
        let doc = Document::replay(&log);
        Ok(Session { log, doc })
    }

    /// Serializes the log (the native format).
    pub fn save(&self) -> Result<Vec<u8>, OpLogError> {
        self.log.encode()
    }

    /// The derived document.
    pub fn document(&self) -> &Document {
        &self.doc
    }

    /// The underlying log.
    pub fn log(&self) -> &OpLog<ModelOp> {
        &self.log
    }

    /// The first Part Studio's id.
    pub fn first_studio(&self) -> Option<ElementId> {
        self.doc.elements.first().map(|e| e.id)
    }

    /// Applies a user edit. The op is appended even if it conflicts (it is
    /// recorded); returns its id.
    pub fn apply(&mut self, op: ModelOp) -> OpId {
        let id = self.log.append(op);
        // Appended ops sort last, so incremental application equals replay.
        let op = match &self.log.get(id).expect("just appended").kind {
            OpKind::Apply(p) => p.clone(),
            OpKind::Undo { .. } => unreachable!(),
        };
        if let Err(reason) = self.doc.apply(id, &op) {
            self.doc.conflicts.push(Conflict { op: id, reason });
        }
        id
    }

    /// Undoes this session's last action. Returns false if nothing to undo.
    pub fn undo(&mut self) -> bool {
        let r = self.log.undo().is_some();
        if r {
            self.doc = Document::replay(&self.log);
        }
        r
    }

    /// Redoes the last undone action. Returns false if nothing to redo.
    pub fn redo(&mut self) -> bool {
        let r = self.log.redo().is_some();
        if r {
            self.doc = Document::replay(&self.log);
        }
        r
    }

    /// How many of this session's actions can be undone; undoing back to an
    /// earlier depth reverts everything done since.
    pub fn undo_depth(&self) -> usize {
        self.log.undo_depth()
    }

    /// True if there is something to undo.
    pub fn can_undo(&self) -> bool {
        self.log.can_undo()
    }

    /// True if there is something to redo.
    pub fn can_redo(&self) -> bool {
        self.log.can_redo()
    }

    /// Integrates ops from another replica (sync) and re-derives the document.
    pub fn integrate(&mut self, ops: Vec<Op<ModelOp>>) -> Result<(), OpLogError> {
        for op in ops {
            self.log.integrate(op)?;
        }
        self.doc = Document::replay(&self.log);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn box_feature(s: f64) -> Feature {
        Feature::Box {
            corner: DVec3::ZERO,
            size: DVec3::splat(s),
        }
    }

    #[test]
    fn insert_edit_suppress_move_delete_undo() {
        let mut s = Session::new_document(ReplicaId(1));
        assert!(!s.can_undo(), "genesis is not undoable");
        let el = s.first_studio().unwrap();
        let a = FeatureId(s.apply(ModelOp::InsertFeature {
            element: el,
            at: InsertAt::End,
            name: "A".into(),
            feature: box_feature(1.0),
        }));
        let b = FeatureId(s.apply(ModelOp::InsertFeature {
            element: el,
            at: InsertAt::End,
            name: "B".into(),
            feature: box_feature(2.0),
        }));
        let names = |s: &Session| -> Vec<String> {
            s.document().elements[0]
                .studio
                .features
                .iter()
                .map(|f| f.name.clone())
                .collect()
        };
        assert_eq!(names(&s), ["A", "B"]);
        s.apply(ModelOp::MoveFeature {
            element: el,
            feature: b,
            to: InsertAt::Start,
        });
        assert_eq!(names(&s), ["B", "A"]);
        s.apply(ModelOp::SetSuppressed {
            element: el,
            feature: a,
            suppressed: true,
        });
        assert!(s.document().elements[0].studio.features[1].suppressed);
        s.apply(ModelOp::EditFeature {
            element: el,
            feature: a,
            definition: box_feature(5.0),
        });
        s.apply(ModelOp::DeleteFeature {
            element: el,
            feature: b,
        });
        assert_eq!(names(&s), ["A"]);
        // Undo everything step by step, then redo everything.
        let mut states = vec![s.document().clone()];
        while s.undo() {
            states.push(s.document().clone());
        }
        assert!(s.document().elements[0].studio.features.is_empty());
        assert_eq!(s.document().elements.len(), 1);
        states.pop();
        while s.redo() {
            assert_eq!(Some(s.document()), states.last());
            states.pop();
        }
        assert!(s.document().conflicts.is_empty());
        // Save/load round trip.
        let bytes = s.save().unwrap();
        let back = Session::open(&bytes, ReplicaId(1)).unwrap();
        assert_eq!(back.document(), s.document());
    }

    #[test]
    fn concurrent_edit_of_deleted_feature_is_a_conflict_not_a_crash() {
        let mut a = Session::new_document(ReplicaId(1));
        let el = a.first_studio().unwrap();
        let f = FeatureId(a.apply(ModelOp::InsertFeature {
            element: el,
            at: InsertAt::End,
            name: "A".into(),
            feature: box_feature(1.0),
        }));
        let mut b = Session::open(&a.save().unwrap(), ReplicaId(2)).unwrap();
        a.apply(ModelOp::DeleteFeature {
            element: el,
            feature: f,
        });
        b.apply(ModelOp::EditFeature {
            element: el,
            feature: f,
            definition: box_feature(3.0),
        });
        let to_b = a.log().missing_for(&b.log().heads());
        let to_a = b.log().missing_for(&a.log().heads());
        a.integrate(to_a).unwrap();
        b.integrate(to_b).unwrap();
        assert_eq!(a.document(), b.document());
        assert!(a.document().elements[0].studio.features.is_empty());
        assert_eq!(a.document().conflicts.len(), 1);
    }

    #[test]
    fn model_ops_serialize() {
        let op = ModelOp::InsertFeature {
            element: ElementId(OpId {
                lamport: 1,
                replica: GENESIS,
            }),
            at: InsertAt::End,
            name: "Revolve 1".into(),
            feature: Feature::Revolve {
                profile: ProfileSource::Inline(ProfileDef {
                    plane: PlaneRef::Front,
                    curves: vec![CurveDef::Circle {
                        center: DVec2::new(0.05, 0.0),
                        radius: 0.01,
                    }],
                }),
                axis: RevolveAxis::Inline {
                    origin: DVec2::ZERO,
                    dir: DVec2::Y,
                },
                angle: 1.0,
            },
        };
        let json = serde_json::to_string(&op).unwrap();
        let back: ModelOp = serde_json::from_str(&json).unwrap();
        assert_eq!(back, op);
    }
}
