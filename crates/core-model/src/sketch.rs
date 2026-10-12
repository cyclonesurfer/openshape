//! Sketches in the document (ADR 0011).
//!
//! A sketch feature stores its plane, its entities (points and curves) and
//! its constraints, keyed by [`SketchId`]s, together with the last solved
//! positions (the solver's starting point). Users edit sketches through
//! [`SketchEdit`] batches carried by `ModelOp::EditSketch`; a batch applies
//! atomically and is one undo step.

use std::collections::BTreeSet;

use glam::DVec2;
use openshape_oplog::OpId;
use serde::{Deserialize, Serialize};

pub use openshape_sketch_solver::{
    Constraint, ConstraintEntry, DistanceKind, Entity, Geometry, Sketch, TangentKind,
};

use crate::PlaneRef;

/// Identifies a sketch entity or constraint: the op that created it and its
/// index within that op's batch.
///
/// Inside a [`SketchEdit`] batch, `op == OpId::NULL` means "the op being
/// applied", so a batch can create items and refer to them; the reducer
/// substitutes the real op id.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct SketchId {
    /// Creating op.
    pub op: OpId,
    /// Index within the op's batch.
    pub index: u32,
}

impl SketchId {
    /// A new item's id within the batch being built (see the type docs).
    pub const fn new_in_batch(index: u32) -> SketchId {
        SketchId {
            op: OpId::NULL,
            index,
        }
    }

    /// True for ids created by the batch being applied.
    pub fn is_new_in_batch(&self) -> bool {
        self.op == OpId::NULL
    }

    fn resolved(self, op: OpId) -> SketchId {
        if self.is_new_in_batch() {
            SketchId {
                op,
                index: self.index,
            }
        } else {
            self
        }
    }
}

/// The sketch data of the document model.
pub type SketchData = Sketch<SketchId>;

/// A sketch feature's definition.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SketchDef {
    /// The sketch plane.
    pub plane: PlaneRef,
    /// Entities and constraints, with their last solved positions.
    pub sketch: SketchData,
}

/// One edit within a sketch edit batch.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum SketchEdit {
    /// Adds an entity. Points must be added before curves that use them.
    AddEntity(Entity<SketchId>),
    /// Adds a constraint or dimension.
    AddConstraint(ConstraintEntry<SketchId>),
    /// Changes a dimension's value.
    SetValue {
        /// The dimension.
        constraint: SketchId,
        /// New value (metres or radians).
        value: f64,
    },
    /// Deletes entities and constraints, together with everything that
    /// references them (curves using a deleted point, constraints on a
    /// deleted entity; a mirror, offset or pattern loses only the pairs,
    /// pieces or copies of deleted entities, `Constraint::without`). Items
    /// that no longer exist are ignored, so concurrent deletes merge.
    Delete {
        /// Items to delete.
        items: Vec<SketchId>,
    },
    /// Marks an entity as construction geometry or not.
    SetConstruction {
        /// The entity.
        entity: SketchId,
        /// New flag.
        construction: bool,
    },
    /// Moves points (usually the solver's result after an edit).
    MovePoints {
        /// Points and their new positions.
        points: Vec<(SketchId, DVec2)>,
    },
    /// Changes a circle's radius (usually the solver's result).
    SetRadius {
        /// The circle.
        entity: SketchId,
        /// New radius (> 0).
        radius: f64,
    },
    /// Changes an interpolated spline's parameters (re-parameterization
    /// after its points moved).
    SetSplineParams {
        /// The spline.
        entity: SketchId,
        /// New parameters: strictly increasing, first 0, last 1.
        params: Vec<f64>,
    },
    /// Moves the sketch to another plane.
    SetPlane(PlaneRef),
}

impl SketchDef {
    /// Applies a batch created by op `op`, atomically: on error the sketch
    /// is unchanged and the error names the failing edit.
    pub fn apply_edits(&mut self, op: OpId, edits: &[SketchEdit]) -> Result<(), String> {
        let mut next = self.clone();
        for (i, edit) in edits.iter().enumerate() {
            next.apply_one(op, edit)
                .map_err(|e| format!("sketch edit {}: {e}", i + 1))?;
        }
        *self = next;
        Ok(())
    }

    fn all_ids(&self) -> BTreeSet<SketchId> {
        self.sketch
            .entities
            .iter()
            .map(|e| e.id)
            .chain(self.sketch.constraints.iter().map(|c| c.id))
            .collect()
    }

    fn entity_mut(&mut self, id: SketchId) -> Result<&mut Entity<SketchId>, String> {
        self.sketch
            .entities
            .iter_mut()
            .find(|e| e.id == id)
            .ok_or_else(|| format!("entity {id:?} does not exist"))
    }

    fn apply_one(&mut self, op: OpId, edit: &SketchEdit) -> Result<(), String> {
        let r = |id: SketchId| id.resolved(op);
        match edit {
            SketchEdit::AddEntity(e) => {
                let e = resolve_entity(e, op);
                if self.all_ids().contains(&e.id) {
                    return Err(format!("id {:?} is already used", e.id));
                }
                check_geometry(&e.geometry)?;
                for p in e.geometry.referenced_points() {
                    match self.sketch.entity(p) {
                        Some(Entity {
                            geometry: Geometry::Point { .. },
                            ..
                        }) => {}
                        Some(_) => return Err(format!("{p:?} is not a point")),
                        None => return Err(format!("point {p:?} does not exist")),
                    }
                }
                self.sketch.entities.push(e);
                Ok(())
            }
            SketchEdit::AddConstraint(c) => {
                let c = resolve_constraint(c, op);
                if self.all_ids().contains(&c.id) {
                    return Err(format!("id {:?} is already used", c.id));
                }
                if let Some(v) = c.constraint.value() {
                    if !v.is_finite() {
                        return Err("dimension value is not finite".into());
                    }
                }
                if let Constraint::Fix { at, .. } = c.constraint {
                    if !at.is_finite() {
                        return Err("fixed position is not finite".into());
                    }
                }
                for e in c.constraint.referenced() {
                    if self.sketch.entity(e).is_none() {
                        return Err(format!("entity {e:?} does not exist"));
                    }
                }
                self.sketch.constraints.push(c);
                Ok(())
            }
            SketchEdit::SetValue { constraint, value } => {
                let id = r(*constraint);
                if !value.is_finite() {
                    return Err("dimension value is not finite".into());
                }
                let c = self
                    .sketch
                    .constraints
                    .iter_mut()
                    .find(|c| c.id == id)
                    .ok_or_else(|| format!("constraint {id:?} does not exist"))?;
                if c.constraint.set_value(*value) {
                    Ok(())
                } else {
                    Err(format!("{id:?} is not a dimension"))
                }
            }
            SketchEdit::Delete { items } => {
                let mut gone: BTreeSet<SketchId> = items.iter().map(|i| r(*i)).collect();
                // Close over references: curves using deleted points, then
                // anything referencing deleted entities.
                loop {
                    let before = gone.len();
                    for e in &self.sketch.entities {
                        if !gone.contains(&e.id)
                            && e.geometry
                                .referenced_points()
                                .iter()
                                .any(|p| gone.contains(p))
                        {
                            gone.insert(e.id);
                        }
                    }
                    if gone.len() == before {
                        break;
                    }
                }
                self.sketch.entities.retain(|e| !gone.contains(&e.id));
                self.sketch.constraints = std::mem::take(&mut self.sketch.constraints)
                    .into_iter()
                    .filter(|c| !gone.contains(&c.id))
                    .filter_map(|c| {
                        let constraint = c.constraint.without(|e| gone.contains(&e))?;
                        Some(ConstraintEntry {
                            id: c.id,
                            constraint,
                        })
                    })
                    .collect();
                Ok(())
            }
            SketchEdit::SetConstruction {
                entity,
                construction,
            } => {
                self.entity_mut(r(*entity))?.construction = *construction;
                Ok(())
            }
            SketchEdit::MovePoints { points } => {
                for (id, to) in points {
                    if !to.is_finite() {
                        return Err("point position is not finite".into());
                    }
                    let id = r(*id);
                    match &mut self.entity_mut(id)?.geometry {
                        Geometry::Point { at } => *at = *to,
                        _ => return Err(format!("{id:?} is not a point")),
                    }
                }
                Ok(())
            }
            SketchEdit::SetRadius { entity, radius } => {
                if !(radius.is_finite() && *radius > 0.0) {
                    return Err("radius must be positive".into());
                }
                let id = r(*entity);
                match &mut self.entity_mut(id)?.geometry {
                    Geometry::Circle { radius: r0, .. } => {
                        *r0 = *radius;
                        Ok(())
                    }
                    _ => Err(format!("{id:?} is not a circle")),
                }
            }
            SketchEdit::SetSplineParams { entity, params } => {
                let id = r(*entity);
                match &mut self.entity_mut(id)?.geometry {
                    Geometry::InterpolatedSpline { through, params: p } => {
                        check_params(params, through.len())?;
                        *p = params.clone();
                        Ok(())
                    }
                    _ => Err(format!("{id:?} is not an interpolated spline")),
                }
            }
            SketchEdit::SetPlane(plane) => {
                self.plane = *plane;
                Ok(())
            }
        }
    }
}

fn resolve_entity(e: &Entity<SketchId>, op: OpId) -> Entity<SketchId> {
    Entity {
        id: e.id.resolved(op),
        geometry: e.geometry.map_keys(|id| id.resolved(op)),
        construction: e.construction,
    }
}

fn resolve_constraint(c: &ConstraintEntry<SketchId>, op: OpId) -> ConstraintEntry<SketchId> {
    ConstraintEntry {
        id: c.id.resolved(op),
        constraint: c.constraint.map_keys(|id| id.resolved(op)),
    }
}

fn check_params(params: &[f64], points: usize) -> Result<(), String> {
    let ok = params.len() == points
        && points >= 2
        && params.first() == Some(&0.0)
        && params.last() == Some(&1.0)
        && params.windows(2).all(|w| w[1] > w[0])
        && params.iter().all(|p| p.is_finite());
    if ok {
        Ok(())
    } else {
        Err("spline parameters must increase strictly from 0 to 1, one per point".into())
    }
}

/// Structural checks that need no other entity.
fn check_geometry(g: &Geometry<SketchId>) -> Result<(), String> {
    match g {
        Geometry::Point { at } if !at.is_finite() => Err("point position is not finite".into()),
        Geometry::Line { start, end } if start == end => Err("a line needs two points".into()),
        Geometry::Circle { radius, .. } if !(radius.is_finite() && *radius > 0.0) => {
            Err("radius must be positive".into())
        }
        Geometry::Arc { center, start, end }
            if center == start || center == end || start == end =>
        {
            Err("an arc needs three distinct points".into())
        }
        Geometry::ControlSpline { control, degree } => {
            if !(1..=5).contains(degree) {
                Err("spline degree must be 1 to 5".into())
            } else if control.len() < *degree as usize + 1 {
                Err(format!(
                    "a degree {degree} spline needs {} control points",
                    degree + 1
                ))
            } else {
                Ok(())
            }
        }
        Geometry::InterpolatedSpline { through, params } => check_params(params, through.len()),
        _ => Ok(()),
    }
}

/// Names a region of a sketch's planar arrangement by the sketch items on
/// its boundary (design note 0006 §5), with an interior point as a
/// geometric tie-break.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RegionKey {
    /// Boundary loops: the outer loop first, then holes in a canonical order;
    /// each loop rotated to start at its smallest edge key.
    pub loops: Vec<Vec<EdgeKey>>,
    /// A point inside the region when it was selected (plane coordinates).
    pub interior: DVec2,
}

/// One boundary piece of a [`RegionKey`].
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EdgeKey {
    /// The sketch curve.
    pub entity: SketchId,
    /// True when the region boundary runs along the curve's own direction.
    pub forward: bool,
    /// Where the piece starts.
    pub from: VertexKey,
    /// Where the piece ends.
    pub to: VertexKey,
}

/// A vertex of the arrangement, named by the curves meeting there.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct VertexKey {
    /// The sketch curves through the vertex, sorted.
    pub entities: Vec<SketchId>,
    /// Distinguishes several vertices with the same curves, in order along
    /// the first curve.
    pub ordinal: u32,
}
