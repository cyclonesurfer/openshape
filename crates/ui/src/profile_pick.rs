//! Picking sketch regions (and a revolve axis) for extrude and revolve.
//!
//! Like [`crate::sketch_mode`], this is UI-toolkit independent: the app
//! shows a small window with the parameters and forwards viewport clicks
//! here.

use std::collections::BTreeMap;

use egui::Pos2;
use openshape_core_model::units::{LengthUnit, Quantity};
use openshape_core_model::{
    DVec2, Feature, FeatureEntry, FeatureId, ProfileSource, RegionKey, RevolveAxis, SketchId,
};
use openshape_features::regions::resolve_region;
use openshape_features::SketchResult;
use openshape_kernel::Kernel;
use openshape_sketch_solver::ResolvedCurve;

use crate::sketch_mode::{sample_curve, ScreenMap, SketchDisplay, CURVE_PICK};
use crate::Model;

/// Which sweep is being set up.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SweepKind {
    /// Linear extrusion.
    Extrude,
    /// Revolve about a sketch line.
    Revolve,
}

/// An open extrude or revolve dialog that takes its profile from sketch
/// regions picked in the viewport.
#[derive(Clone, Debug)]
pub struct ProfilePick {
    /// Extrude or revolve.
    pub kind: SweepKind,
    /// The feature being edited, or `None` when creating.
    pub editing: Option<FeatureId>,
    /// The sketch the regions belong to.
    pub sketch: Option<FeatureId>,
    /// Picked regions.
    pub regions: Vec<RegionKey>,
    /// Revolve axis: a line of the sketch.
    pub axis: Option<SketchId>,
    /// The next viewport click picks the axis instead of a region.
    pub picking_axis: bool,
    /// Extrude distance as typed (bare numbers are millimetres).
    pub distance: String,
    /// Extrude to both sides.
    pub symmetric: bool,
    /// Revolve angle as typed (bare numbers are degrees).
    pub angle: String,
    /// Why the last OK failed.
    pub error: Option<String>,
    /// The region under the pointer.
    pub hover: Option<(FeatureId, usize)>,
    outlines: BTreeMap<(FeatureId, usize), Vec<Vec<DVec2>>>,
    outline_generation: u64,
}

impl ProfilePick {
    /// A new, empty dialog.
    pub fn new(kind: SweepKind) -> ProfilePick {
        ProfilePick {
            kind,
            editing: None,
            sketch: None,
            regions: Vec::new(),
            axis: None,
            picking_axis: false,
            distance: "10".into(),
            symmetric: false,
            angle: "360".into(),
            error: None,
            hover: None,
            outlines: BTreeMap::new(),
            outline_generation: u64::MAX,
        }
    }

    /// The dialog for editing a feature whose profile comes from a sketch;
    /// `None` for other features.
    pub fn edit(entry: &FeatureEntry) -> Option<ProfilePick> {
        let mm = |v: f64| trim(v / LengthUnit::Millimetre.metres());
        let mut p = match &entry.feature {
            Feature::Extrude {
                profile: ProfileSource::Sketch { sketch, regions },
                distance,
                symmetric,
            } => ProfilePick {
                sketch: Some(*sketch),
                regions: regions.clone(),
                distance: mm(*distance),
                symmetric: *symmetric,
                ..ProfilePick::new(SweepKind::Extrude)
            },
            Feature::Revolve {
                profile: ProfileSource::Sketch { sketch, regions },
                axis,
                angle,
            } => ProfilePick {
                sketch: Some(*sketch),
                regions: regions.clone(),
                axis: match axis {
                    RevolveAxis::SketchLine { line } => Some(*line),
                    RevolveAxis::Inline { .. } => None,
                },
                angle: trim(angle.to_degrees()),
                ..ProfilePick::new(SweepKind::Revolve)
            },
            _ => return None,
        };
        p.editing = Some(entry.id);
        Some(p)
    }

    /// Window title.
    pub fn title(&self) -> &'static str {
        match (self.kind, self.editing.is_some()) {
            (SweepKind::Extrude, false) => "Extrude",
            (SweepKind::Extrude, true) => "Edit extrude",
            (SweepKind::Revolve, false) => "Revolve",
            (SweepKind::Revolve, true) => "Edit revolve",
        }
    }

    /// What the next viewport click does.
    pub fn hint(&self) -> &'static str {
        if self.picking_axis {
            "Click a line of the sketch to revolve about"
        } else {
            "Click sketch regions to add or remove them"
        }
    }

    /// Indices of the picked regions in the current regeneration (regions
    /// that no longer exist are left out).
    pub fn picked(&self, model: &Model) -> Vec<usize> {
        let Some(regions) = self
            .sketch
            .and_then(|f| model.sketch_result(f))
            .and_then(|s| s.regions.as_ref().ok())
        else {
            return Vec::new();
        };
        self.regions
            .iter()
            .filter_map(|k| resolve_region(model.kernel(), regions, k))
            .collect()
    }

    /// Updates the region under the pointer.
    pub fn hover(&mut self, model: &Model, view: &ScreenMap, pos: Option<Pos2>) {
        self.hover = match pos {
            Some(pos) if !self.picking_axis => region_at(model, view, pos),
            _ => None,
        };
    }

    /// A viewport click: toggles the region under `pos`, or picks the axis.
    pub fn click(&mut self, model: &Model, view: &ScreenMap, pos: Pos2) {
        self.error = None;
        if self.picking_axis {
            match self.sketch.and_then(|f| line_at(model, view, f, pos)) {
                Some(line) => {
                    self.axis = Some(line);
                    self.picking_axis = false;
                }
                None => self.error = Some("Click a line of the profile's sketch".into()),
            }
            return;
        }
        let Some((feature, index)) = region_at(model, view, pos) else {
            return;
        };
        let Some(key) = model
            .sketch_result(feature)
            .and_then(|s| s.regions.as_ref().ok())
            .and_then(|r| r.keys.get(index))
            .cloned()
        else {
            return;
        };
        if self.sketch != Some(feature) {
            // A profile comes from one sketch.
            self.sketch = Some(feature);
            self.regions.clear();
            self.axis = None;
        }
        // Clicking a picked region removes it.
        let regions = model
            .sketch_result(feature)
            .and_then(|s| s.regions.as_ref().ok());
        let existing = regions.and_then(|r| {
            self.regions
                .iter()
                .position(|k| resolve_region(model.kernel(), r, k) == Some(index))
        });
        match existing {
            Some(i) => {
                self.regions.remove(i);
            }
            None => self.regions.push(key),
        }
    }

    /// The feature the dialog describes.
    pub fn build(&self) -> Result<Feature, String> {
        let sketch = self.sketch.ok_or("Pick at least one sketch region")?;
        if self.regions.is_empty() {
            return Err("Pick at least one sketch region".into());
        }
        let profile = ProfileSource::Sketch {
            sketch,
            regions: self.regions.clone(),
        };
        match self.kind {
            SweepKind::Extrude => {
                let distance = Quantity::parse_length(&self.distance, LengthUnit::Millimetre)
                    .map_err(|e| format!("Distance: {e}"))?;
                if distance == 0.0 {
                    return Err("Distance: must not be zero".into());
                }
                Ok(Feature::Extrude {
                    profile,
                    distance,
                    symmetric: self.symmetric,
                })
            }
            SweepKind::Revolve => {
                let line = self.axis.ok_or("Pick an axis line")?;
                let angle =
                    Quantity::parse_angle(&self.angle).map_err(|e| format!("Angle: {e}"))?;
                if angle == 0.0 || angle.abs() > std::f64::consts::TAU + 1e-12 {
                    return Err("Angle: must be non-zero and at most 360°".into());
                }
                Ok(Feature::Revolve {
                    profile,
                    axis: RevolveAxis::SketchLine { line },
                    angle,
                })
            }
        }
    }

    /// Boundary polylines of a region (plane coordinates), cached per
    /// regeneration.
    pub fn outline(&mut self, model: &Model, feature: FeatureId, region: usize) -> &[Vec<DVec2>] {
        if self.outline_generation != model.generation {
            self.outlines.clear();
            self.outline_generation = model.generation;
        }
        self.outlines
            .entry((feature, region))
            .or_insert_with(|| region_outline(model, feature, region))
    }
}

fn trim(x: f64) -> String {
    let s = format!("{x:.6}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Boundary polylines of region `region` of a sketch.
pub fn region_outline(model: &Model, feature: FeatureId, region: usize) -> Vec<Vec<DVec2>> {
    let Some(regions) = model
        .sketch_result(feature)
        .and_then(|s| s.regions.as_ref().ok())
    else {
        return Vec::new();
    };
    match model
        .kernel()
        .region_profile(&regions.arrangement, &regions.curves, &[region])
    {
        Ok(curves) => curves.iter().map(sample_curve).collect(),
        Err(_) => Vec::new(),
    }
}

/// The region of any sketch under `pos`, nearest the viewer.
pub fn region_at(model: &Model, view: &ScreenMap, pos: Pos2) -> Option<(FeatureId, usize)> {
    model
        .sketches()
        .iter()
        .filter_map(|s: &SketchResult| {
            let regions = s.regions.as_ref().ok()?;
            let (p, t) = view.with_frame(s.frame).hit(pos)?;
            let i = model
                .kernel()
                .region_containing(&regions.arrangement, &regions.curves, p)?;
            Some((t, s.feature, i))
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, f, i)| (f, i))
}

/// The line of sketch `feature` (construction lines included) under `pos`.
pub fn line_at(model: &Model, view: &ScreenMap, feature: FeatureId, pos: Pos2) -> Option<SketchId> {
    let s = model.sketch_result(feature)?;
    let d = SketchDisplay::from_result(s);
    let map = view.with_frame(s.frame);
    d.curves
        .iter()
        .filter(|c| matches!(c.curve, ResolvedCurve::Line { .. }))
        .filter_map(|c| {
            let dist = map
                .polyline(&c.polyline)
                .iter()
                .flat_map(|run| run.windows(2).map(|w| segment_distance(pos, w[0], w[1])))
                .fold(f32::INFINITY, f32::min);
            (dist <= CURVE_PICK).then_some((c.id, dist))
        })
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(id, _)| id)
}

fn segment_distance(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let l2 = ab.length_sq();
    let t = if l2 > 0.0 {
        ((p - a).dot(ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance(a + ab * t)
}
