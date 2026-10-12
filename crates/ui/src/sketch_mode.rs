//! Sketch mode: drawing, constraining, dimensioning and dragging sketch
//! geometry.
//!
//! Nothing here uses egui's widgets or painter, only its plain geometry
//! types: [`SketchEditor`] receives pointer and key events in screen points,
//! maps them to the sketch plane with a [`ScreenMap`], and edits the document
//! through [`Model::edit_sketch`]. `app.rs` routes egui input here and paints
//! the [`SketchDisplay`] and its [`Annotation`]s; `tests/sketch_mode.rs`
//! drives the same editor headlessly.

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use egui::{Pos2, Rect, Vec2};
use glam::DMat4;
use openshape_core_model::sketch::{
    Constraint, ConstraintEntry, DistanceKind, Entity, Geometry, TangentKind,
};
use openshape_core_model::units::{LengthUnit, Quantity};
use openshape_core_model::{DVec2, DVec3, FeatureId, OpId, SketchData, SketchEdit, SketchId};
use openshape_features::{kernel_curve, EntityState, SketchResult, Solution};
use openshape_kernel::{Frame, ProfileCurve};
use openshape_render::{Camera, Projection};
use openshape_sketch_solver::{
    conic, measure, resolve_curves, solve, ResolvedCurve, SolveError, SolveOptions, SolveStatus,
    RHO_MIN,
};
use openshape_text::Font;

use crate::sketch_ops::fillet::fillet;
use crate::sketch_ops::shapes::{self, TangentStart};
use crate::sketch_ops::transform::{self, Similarity};
use crate::sketch_ops::trim::{self, Cuts, Plan};
use crate::sketch_ops::{curve_ending_at, deletion, leaving_direction, Batch};
use crate::sketch_ops::{mirror, offset, pattern, text};
use crate::{geometry_edits, Model};

/// Pick radius for points, in screen points.
pub const POINT_PICK: f32 = 7.0;
/// Pick radius for curves, in screen points.
pub const CURVE_PICK: f32 = 5.0;
/// Snap radius for the sketch origin, in screen points.
pub const ORIGIN_PICK: f32 = 8.0;
/// Pointer travel (screen points) that turns a press into a drag.
pub const DRAG_THRESHOLD: f32 = 4.0;
/// A new line within this angle of horizontal or vertical is made so.
const INFER_ANGLE: f64 = 3.0 * PI / 180.0;
/// Positions closer than this are the same click (metres).
const SAME_POINT: f64 = 1e-7;

// ---------------------------------------------------------------------------
// Screen mapping
// ---------------------------------------------------------------------------

/// Maps between a plane and screen points for one frame of the viewport.
#[derive(Clone, Debug)]
pub struct ScreenMap {
    camera: Camera,
    frame: Frame,
    rect: Rect,
    view_proj: DMat4,
}

impl ScreenMap {
    /// The map for `camera` drawing into `rect` (screen points), with plane
    /// coordinates in `frame`.
    pub fn new(camera: &Camera, frame: Frame, rect: Rect) -> ScreenMap {
        let aspect = rect.width() / rect.height().max(1.0);
        // Screen x and y do not depend on the clip planes, so the scene
        // bounds the renderer uses for them are not needed here.
        let view_proj = camera.proj_matrix_f64(aspect, None) * camera.view_matrix_f64();
        ScreenMap {
            camera: camera.clone(),
            frame,
            rect,
            view_proj,
        }
    }

    /// The same view with another plane.
    pub fn with_frame(&self, frame: Frame) -> ScreenMap {
        ScreenMap {
            frame,
            ..self.clone()
        }
    }

    /// The plane.
    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    /// The viewport rectangle.
    pub fn rect(&self) -> Rect {
        self.rect
    }

    /// Projects a world point; `None` behind the eye.
    pub fn world_to_screen(&self, p: DVec3) -> Option<Pos2> {
        let c = self.view_proj * p.extend(1.0);
        if !(c.w > 1e-12 && c.is_finite()) {
            return None;
        }
        let (x, y) = (c.x / c.w, c.y / c.w);
        Some(Pos2::new(
            self.rect.min.x + ((x + 1.0) * 0.5 * self.rect.width() as f64) as f32,
            self.rect.min.y + ((1.0 - y) * 0.5 * self.rect.height() as f64) as f32,
        ))
    }

    /// Projects a plane point.
    pub fn to_screen(&self, p: DVec2) -> Option<Pos2> {
        self.world_to_screen(self.frame.point2(p))
    }

    /// The plane point under screen point `s` and its distance along the
    /// view ray; `None` when the plane is seen edge-on or lies behind the
    /// eye.
    pub fn hit(&self, s: Pos2) -> Option<(DVec2, f64)> {
        let (o, d) = self.camera.pixel_ray(
            s.x - self.rect.min.x,
            s.y - self.rect.min.y,
            self.rect.width(),
            self.rect.height(),
        );
        let n = self.frame.z;
        let denom = d.dot(n);
        if denom.abs() < 1e-6 {
            return None;
        }
        let t = (self.frame.origin - o).dot(n) / denom;
        if matches!(self.camera.projection, Projection::Perspective { .. }) && t <= 0.0 {
            return None;
        }
        let l = self.frame.to_local(o + d * t);
        let p = DVec2::new(l.x, l.y);
        p.is_finite().then_some((p, t))
    }

    /// The plane point under screen point `s`.
    pub fn to_plane(&self, s: Pos2) -> Option<DVec2> {
        self.hit(s).map(|h| h.0)
    }

    /// Plane length of one screen point near `s` (0 if unknown).
    pub fn plane_per_point(&self, s: Pos2) -> f64 {
        match (
            self.to_plane(s),
            self.to_plane(s + Vec2::new(1.0, 0.0)),
            self.to_plane(s + Vec2::new(0.0, 1.0)),
        ) {
            (Some(a), Some(b), Some(c)) => (b - a).length().max((c - a).length()),
            _ => 0.0,
        }
    }

    /// A plane polyline projected to the screen, split where it passes
    /// behind the eye.
    pub fn polyline(&self, points: &[DVec2]) -> Vec<Vec<Pos2>> {
        let mut out = Vec::new();
        let mut run = Vec::new();
        for p in points {
            match self.to_screen(*p) {
                Some(s) => run.push(s),
                None => {
                    if run.len() > 1 {
                        out.push(std::mem::take(&mut run));
                    }
                    run.clear();
                }
            }
        }
        if run.len() > 1 {
            out.push(run);
        }
        out
    }
}

/// Camera yaw and pitch that look straight at a plane, with the plane's x
/// axis to the right where possible.
pub fn normal_view(frame: &Frame) -> (f64, f64) {
    let z = frame.z;
    let pitch = z.z.clamp(-1.0, 1.0).asin();
    let yaw = if z.x.hypot(z.y) < 1e-9 {
        // Looking straight down or up: choose the yaw that keeps the
        // plane's x axis pointing right.
        frame.x.y.atan2(frame.x.x)
    } else {
        z.x.atan2(-z.y)
    };
    (yaw, pitch)
}

// ---------------------------------------------------------------------------
// Display geometry
// ---------------------------------------------------------------------------

/// Number of segments used to draw a curve.
fn segments(c: &ProfileCurve) -> usize {
    match c {
        ProfileCurve::Line { .. } => 1,
        ProfileCurve::Circle { .. } => 96,
        ProfileCurve::Arc { sweep_angle, .. } => {
            ((sweep_angle.abs() / (3.75f64).to_radians()).ceil() as usize).clamp(4, 96)
        }
        ProfileCurve::Spline { points, .. } => (points.len() * 16).clamp(32, 1024),
    }
}

/// Points along a curve for drawing and picking (empty if the curve is
/// invalid).
pub fn sample_curve(c: &ProfileCurve) -> Vec<DVec2> {
    c.sample(segments(c)).unwrap_or_default()
}

/// Points along a resolved sketch curve (empty if it has no kernel curve,
/// e.g. a degenerate ellipse).
pub fn sample_resolved(c: &ResolvedCurve) -> Vec<DVec2> {
    if let ResolvedCurve::Conic {
        start,
        control,
        end,
        weight,
    } = *c
    {
        // Denser where it bends: a full conic turns sharply near its
        // control point (design note 0009 §4.2).
        let p = [start, control, end];
        return conic::samples(p, weight)
            .into_iter()
            .map(|t| conic::point(p, weight, t))
            .collect();
    }
    kernel_curve(c)
        .map(|k| sample_curve(&k))
        .unwrap_or_default()
}

/// One sketch curve as drawn.
#[derive(Clone, Debug)]
pub struct DisplayCurve {
    /// The entity.
    pub id: SketchId,
    /// Construction geometry.
    pub construction: bool,
    /// Explicit geometry.
    pub curve: ResolvedCurve,
    /// Points along the curve, in plane coordinates.
    pub polyline: Vec<DVec2>,
}

/// A sketch's geometry and diagnostics, ready to draw and pick.
#[derive(Clone, Debug)]
pub struct SketchDisplay {
    /// The sketch feature.
    pub feature: FeatureId,
    /// The sketch plane.
    pub frame: Frame,
    /// The geometry shown (solved, or the solver's input when it failed).
    pub sketch: SketchData,
    /// Curves.
    pub curves: Vec<DisplayCurve>,
    /// Free or constrained, per entity.
    pub states: BTreeMap<SketchId, EntityState>,
    /// Constraints that cannot be satisfied together.
    pub conflicting: BTreeSet<SketchId>,
    /// Entities referenced by conflicting constraints.
    pub conflicted: BTreeSet<SketchId>,
    /// Constraints implied by others.
    pub redundant: BTreeSet<SketchId>,
    /// Solver outcome.
    pub status: SolveStatus,
    /// Remaining degrees of freedom.
    pub dof: usize,
    /// Placed dimension labels: offsets from each dimension's anchor
    /// (`SketchDef::labels`, or the label being dragged).
    pub labels: BTreeMap<SketchId, DVec2>,
}

impl SketchDisplay {
    /// The display of a regenerated sketch.
    pub fn from_result(r: &SketchResult) -> SketchDisplay {
        Self::build(r.feature, r.frame, &r.solution, Some(&r.curves))
    }

    /// The display of a solution (a drag preview).
    pub fn from_solution(
        feature: FeatureId,
        frame: Frame,
        solution: &Solution<SketchId>,
    ) -> SketchDisplay {
        Self::build(feature, frame, solution, None)
    }

    fn build(
        feature: FeatureId,
        frame: Frame,
        solution: &Solution<SketchId>,
        curves: Option<&[(SketchId, bool, ResolvedCurve)]>,
    ) -> SketchDisplay {
        let resolved;
        let curves = match curves {
            Some(c) => c,
            None => {
                resolved = resolve_curves(&solution.sketch).unwrap_or_default();
                &resolved
            }
        };
        let curves = curves
            .iter()
            .map(|(id, construction, c)| DisplayCurve {
                id: *id,
                construction: *construction,
                polyline: sample_resolved(c),
                curve: c.clone(),
            })
            .collect();
        let conflicting: BTreeSet<SketchId> = solution.conflicting.iter().copied().collect();
        let conflicted = solution
            .sketch
            .constraints
            .iter()
            .filter(|c| conflicting.contains(&c.id))
            .flat_map(|c| c.constraint.referenced())
            .collect();
        SketchDisplay {
            feature,
            frame,
            sketch: solution.sketch.clone(),
            curves,
            states: solution.entities.clone(),
            conflicting,
            conflicted,
            redundant: solution.redundant.iter().copied().collect(),
            status: solution.status.clone(),
            dof: solution.dof,
            labels: BTreeMap::new(),
        }
    }

    /// Point entities and their positions.
    pub fn points(&self) -> impl Iterator<Item = (SketchId, DVec2)> + '_ {
        self.sketch
            .entities
            .iter()
            .filter_map(|e| match e.geometry {
                Geometry::Point { at } => Some((e.id, at)),
                _ => None,
            })
    }

    /// The drawn curve of an entity.
    pub fn curve(&self, id: SketchId) -> Option<&DisplayCurve> {
        self.curves.iter().find(|c| c.id == id)
    }

    /// Where a dimension's label offset is measured from: the middle of
    /// what it measures (design note 0007 §12). `None` for constraints
    /// that are not dimensions.
    pub fn dimension_anchor(&self, dimension: SketchId) -> Option<DVec2> {
        let points: BTreeMap<SketchId, DVec2> = self.points().collect();
        let c = self.sketch.constraint(dimension)?;
        Some(measured(&self.sketch, &points, &c.constraint)?.anchor())
    }
}

// ---------------------------------------------------------------------------
// Constraint glyphs and dimension labels
// ---------------------------------------------------------------------------

/// A constraint glyph or dimension label, in screen space.
#[derive(Clone, Debug, PartialEq)]
pub struct Annotation {
    /// The constraint.
    pub constraint: SketchId,
    /// Text.
    pub text: String,
    /// Where the text box is.
    pub rect: Rect,
    /// Dimension and leader lines.
    pub lines: Vec<[Pos2; 2]>,
    /// True for dimensions.
    pub dimension: bool,
}

/// The glyph drawn for a (non-dimension) constraint.
pub fn glyph(c: &Constraint<SketchId>) -> &'static str {
    match c {
        Constraint::Coincident { .. } => "•",
        Constraint::PointOnCurve { .. } => "○",
        Constraint::Horizontal { .. } | Constraint::HorizontalPoints { .. } => "H",
        Constraint::Vertical { .. } | Constraint::VerticalPoints { .. } => "V",
        Constraint::Parallel { .. } => "‖",
        Constraint::Perpendicular { .. } => "⟂",
        Constraint::Equal { .. } => "=",
        Constraint::Midpoint { .. } => "M",
        Constraint::Concentric { .. } => "◎",
        Constraint::Symmetric { .. } => "S",
        Constraint::Tangent { .. } => "T",
        Constraint::Curvature { .. } => "C",
        Constraint::Normal { .. } => "N",
        Constraint::Fix { .. } => "F",
        Constraint::Mirror { .. } => "↔",
        Constraint::LinearPattern { .. } => "⇉",
        Constraint::Distance { .. }
        | Constraint::Length { .. }
        | Constraint::Angle { .. }
        | Constraint::Radius { .. }
        | Constraint::Diameter { .. }
        | Constraint::Offset { .. }
        | Constraint::CircularPattern { .. }
        | Constraint::Rho { .. } => "",
    }
}

/// Entities a glyph is drawn next to.
fn glyph_anchors(c: &Constraint<SketchId>) -> Vec<SketchId> {
    match c {
        Constraint::Coincident { a, .. }
        | Constraint::HorizontalPoints { a, .. }
        | Constraint::VerticalPoints { a, .. } => vec![*a],
        Constraint::PointOnCurve { point, .. } | Constraint::Midpoint { point, .. } => {
            vec![*point]
        }
        Constraint::Fix { point, .. } => vec![*point],
        Constraint::Horizontal { line } | Constraint::Vertical { line } => vec![*line],
        Constraint::Normal { line, .. } => vec![*line],
        Constraint::Mirror { axis, .. } => vec![*axis],
        Constraint::LinearPattern { direction, .. } => vec![*direction],
        Constraint::Parallel { a, b }
        | Constraint::Perpendicular { a, b }
        | Constraint::Equal { a, b }
        | Constraint::Concentric { a, b }
        | Constraint::Symmetric { a, b, .. }
        | Constraint::Tangent { a, b, .. }
        | Constraint::Curvature { a, b, .. } => vec![*a, *b],
        _ => Vec::new(),
    }
}

/// Approximate size of a label (screen points), matching the 13-point font
/// `app.rs` draws with.
pub fn label_size(text: &str) -> Vec2 {
    Vec2::new(7.5 * text.chars().count() as f32 + 8.0, 17.0)
}

/// Formats a length for display.
pub fn format_length(metres: f64) -> String {
    LengthUnit::Millimetre.format(metres)
}

/// Formats an angle (shown unsigned, in degrees).
pub fn format_angle(radians: f64) -> String {
    format!("{}°", trim_number(radians.abs().to_degrees()))
}

fn trim_number(x: f64) -> String {
    let s = format!("{x:.4}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" {
        "0".into()
    } else {
        s.to_string()
    }
}

/// Text of a dimension.
pub fn dimension_text(c: &Constraint<SketchId>) -> String {
    match c {
        Constraint::Distance { value, .. } | Constraint::Length { value, .. } => {
            format_length(*value)
        }
        Constraint::Radius { value, .. } => format!("R {}", format_length(*value)),
        Constraint::Diameter { value, .. } => format!("Ø {}", format_length(*value)),
        Constraint::Angle { value, .. } | Constraint::CircularPattern { angle: value, .. } => {
            format_angle(*value)
        }
        Constraint::Offset { distance, .. } => format_length(*distance),
        Constraint::Rho { value, .. } => format!("ρ {}", trim_number(*value)),
        _ => String::new(),
    }
}

fn perp(v: DVec2) -> DVec2 {
    DVec2::new(-v.y, v.x)
}

/// Screen midpoint (by length) of a projected polyline, with its unit
/// direction there.
fn screen_midpoint(runs: &[Vec<Pos2>]) -> Option<(Pos2, Vec2)> {
    let run = runs.iter().max_by(|a, b| a.len().cmp(&b.len()))?;
    let total: f32 = run.windows(2).map(|w| w[0].distance(w[1])).sum();
    let mut acc = 0.0;
    for w in run.windows(2) {
        let l = w[0].distance(w[1]);
        if acc + l >= total * 0.5 && l > 0.0 {
            let f = (total * 0.5 - acc) / l;
            return Some((w[0] + (w[1] - w[0]) * f, (w[1] - w[0]) / l));
        }
        acc += l;
    }
    Some((run[0], Vec2::X))
}

/// Constraint glyphs and dimension labels of a sketch.
pub fn annotations(d: &SketchDisplay, map: &ScreenMap) -> Vec<Annotation> {
    let mut slots: BTreeMap<SketchId, usize> = BTreeMap::new();
    let points: BTreeMap<SketchId, DVec2> = d.points().collect();
    let centre = (!points.is_empty())
        .then(|| points.values().copied().sum::<DVec2>() / points.len() as f64)
        .and_then(|c| map.to_screen(c));
    let mut out = Vec::new();
    for c in &d.sketch.constraints {
        if c.constraint.value().is_some() {
            out.extend(dimension_annotation(d, &points, map, c));
            continue;
        }
        let text = glyph(&c.constraint);
        for anchor in glyph_anchors(&c.constraint) {
            let slot = slots.entry(anchor).or_default();
            let k = *slot as f32;
            *slot += 1;
            let center = if let Some(p) = points.get(&anchor) {
                map.to_screen(*p)
                    .map(|s| s + Vec2::new(12.0 + 15.0 * k, -12.0))
            } else if let Some(curve) = d.curve(anchor) {
                // Inside, towards the sketch's centre (dimensions go
                // outside); above when the curve passes through it.
                screen_midpoint(&map.polyline(&curve.polyline)).map(|(m, t)| {
                    let mut n = Vec2::new(-t.y, t.x);
                    let inward = centre.map_or(0.0, |c| (c - m).dot(n));
                    if inward.abs() > 0.5 {
                        if inward < 0.0 {
                            n = -n;
                        }
                    } else if n.y > 0.0 {
                        n = -n;
                    }
                    m + n * 13.0 + t * (15.0 * k)
                })
            } else {
                None
            };
            if let Some(center) = center {
                out.push(Annotation {
                    constraint: c.id,
                    text: text.to_string(),
                    rect: Rect::from_center_size(center, label_size(text)),
                    lines: Vec::new(),
                    dimension: false,
                });
            }
        }
    }
    out
}

fn line_ends(
    sketch: &SketchData,
    points: &BTreeMap<SketchId, DVec2>,
    id: SketchId,
) -> Option<(DVec2, DVec2)> {
    match sketch.entity(id)?.geometry {
        Geometry::Line { start, end } => Some((*points.get(&start)?, *points.get(&end)?)),
        _ => None,
    }
}

/// Centre, radius and a direction to put a radial label at.
fn circular(
    sketch: &SketchData,
    points: &BTreeMap<SketchId, DVec2>,
    id: SketchId,
) -> Option<(DVec2, f64, DVec2)> {
    match sketch.entity(id)?.geometry {
        Geometry::Circle { center, radius } => Some((
            *points.get(&center)?,
            radius,
            DVec2::new(1.0, 1.0).normalize(),
        )),
        Geometry::Arc { center, start, end } => {
            let (c, s, e) = (
                *points.get(&center)?,
                *points.get(&start)?,
                *points.get(&end)?,
            );
            let a0 = (s - c).to_angle();
            let sweep = ((e - c).to_angle() - a0).rem_euclid(TAU);
            let mid = a0 + sweep * 0.5;
            Some((c, (s - c).length(), DVec2::from_angle(mid)))
        }
        _ => None,
    }
}

/// A conic's start, control point and end, and its ρ.
fn conic_of(
    sketch: &SketchData,
    points: &BTreeMap<SketchId, DVec2>,
    id: SketchId,
) -> Option<([DVec2; 3], f64)> {
    match sketch.entity(id)?.geometry {
        Geometry::Conic {
            start,
            control,
            end,
            rho,
        } => Some((
            [
                *points.get(&start)?,
                *points.get(&control)?,
                *points.get(&end)?,
            ],
            rho,
        )),
        _ => None,
    }
}

/// What a dimension measures, in plane coordinates: enough to draw it and to
/// find its anchor.
#[derive(Clone, Copy, Debug)]
enum Measured {
    /// Between two points along their direction (a length, an aligned
    /// distance).
    Aligned(DVec2, DVec2),
    /// Between two points along x.
    Horizontal(DVec2, DVec2),
    /// Between two points along y.
    Vertical(DVec2, DVec2),
    /// From a point (or a line's midpoint) to its foot on a line, whose
    /// nearer end is `near`.
    Foot { p: DVec2, foot: DVec2, near: DVec2 },
    /// A circle's or arc's radius or diameter; `dir` is where an automatic
    /// label goes.
    Radial {
        centre: DVec2,
        radius: f64,
        dir: DVec2,
        diameter: bool,
    },
    /// The angle between two lines meeting at `vertex`, turning `sweep`
    /// from direction angle `start`; `a` and `b` are the lines.
    Angle {
        vertex: DVec2,
        start: f64,
        sweep: f64,
        a: (DVec2, DVec2),
        b: (DVec2, DVec2),
    },
    /// An offset's distance, from a source to its image.
    Across { from: DVec2, to: DVec2 },
    /// A circular pattern's angle, about `centre` from `source` (a point of
    /// its first source).
    Turn {
        centre: DVec2,
        source: DVec2,
        sweep: f64,
    },
    /// A conic's ρ: its shoulder point, and the unit direction of the
    /// median from the chord's midpoint through it (design note 0009 §7).
    Shoulder { at: DVec2, dir: DVec2 },
}

impl Measured {
    /// The middle of what is measured: a placed label's offset is from here
    /// (design note 0007 §12).
    fn anchor(&self) -> DVec2 {
        match *self {
            Measured::Aligned(a, b) | Measured::Horizontal(a, b) | Measured::Vertical(a, b) => {
                (a + b) * 0.5
            }
            Measured::Foot { p, foot, .. } => (p + foot) * 0.5,
            Measured::Radial { centre, .. } | Measured::Turn { centre, .. } => centre,
            Measured::Angle { vertex, .. } => vertex,
            Measured::Across { from, to } => (from + to) * 0.5,
            Measured::Shoulder { at, .. } => at,
        }
    }
}

fn measured(
    sketch: &SketchData,
    points: &BTreeMap<SketchId, DVec2>,
    c: &Constraint<SketchId>,
) -> Option<Measured> {
    Some(match c {
        Constraint::Length { line, .. } => {
            let (a, b) = line_ends(sketch, points, *line)?;
            Measured::Aligned(a, b)
        }
        Constraint::Distance { a, b, kind, .. } => {
            match (points.get(a).copied(), points.get(b).copied()) {
                (Some(pa), Some(pb)) => match kind {
                    DistanceKind::Aligned => Measured::Aligned(pa, pb),
                    DistanceKind::Horizontal => Measured::Horizontal(pa, pb),
                    DistanceKind::Vertical => Measured::Vertical(pa, pb),
                },
                _ => {
                    // Point to line, or line to line: from a point (or the
                    // first line's midpoint) to its foot on the other line.
                    let (p, line) = match (points.get(a), points.get(b)) {
                        (Some(p), None) => (*p, *b),
                        (None, Some(p)) => (*p, *a),
                        _ => {
                            let (a0, a1) = line_ends(sketch, points, *a)?;
                            ((a0 + a1) * 0.5, *b)
                        }
                    };
                    let (l0, l1) = line_ends(sketch, points, line)?;
                    let dir = (l1 - l0).normalize_or_zero();
                    let foot = l0 + dir * (p - l0).dot(dir);
                    let near = if (foot - l0).length() < (foot - l1).length() {
                        l0
                    } else {
                        l1
                    };
                    Measured::Foot { p, foot, near }
                }
            }
        }
        Constraint::Radius { curve, .. } | Constraint::Diameter { curve, .. } => {
            let (centre, radius, dir) = circular(sketch, points, *curve)?;
            Measured::Radial {
                centre,
                radius,
                dir,
                diameter: matches!(c, Constraint::Diameter { .. }),
            }
        }
        Constraint::Angle { a, b, value } => {
            let (a0, a1) = line_ends(sketch, points, *a)?;
            let (b0, b1) = line_ends(sketch, points, *b)?;
            let (da, db) = (a1 - a0, b1 - b0);
            let cross = da.perp_dot(db);
            let vertex = if cross.abs() > 1e-12 * da.length() * db.length() {
                a0 + da * ((b0 - a0).perp_dot(db) / cross)
            } else {
                (a0 + a1) * 0.5
            };
            Measured::Angle {
                vertex,
                start: da.to_angle(),
                sweep: *value,
                a: (a0, a1),
                b: (b0, b1),
            }
        }
        Constraint::Offset { pieces, .. } => {
            // From the middle of the first source square across to its
            // image.
            let &(src, img, _) = pieces.first()?;
            let (from, to) = match (
                line_ends(sketch, points, src),
                line_ends(sketch, points, img),
            ) {
                (Some((a, b)), Some((p, q))) => {
                    let m = (a + b) * 0.5;
                    let u = (q - p).normalize_or_zero();
                    (m, p + u * (m - p).dot(u))
                }
                _ => {
                    let (c, r, dir) = circular(sketch, points, src)?;
                    let (_, r2, _) = circular(sketch, points, img)?;
                    (c + dir * r, c + dir * r2)
                }
            };
            Measured::Across { from, to }
        }
        Constraint::CircularPattern {
            centre,
            angle,
            copies,
        } => {
            let c = *points.get(centre)?;
            let &(src, _, _) = copies.first()?;
            let source = match sketch.entity(src)?.geometry {
                Geometry::Point { .. } => *points.get(&src)?,
                Geometry::Line { .. } => {
                    let (a, b) = line_ends(sketch, points, src)?;
                    (a + b) * 0.5
                }
                Geometry::Circle { center, .. } => *points.get(&center)?,
                Geometry::Arc { .. } => {
                    let (cc, r, dir) = circular(sketch, points, src)?;
                    cc + dir * r
                }
                ref g => {
                    let pts = g.referenced_points();
                    let sum: DVec2 = pts.iter().filter_map(|q| points.get(q)).sum();
                    sum / pts.len().max(1) as f64
                }
            };
            Measured::Turn {
                centre: c,
                source,
                sweep: *angle,
            }
        }
        Constraint::Rho { curve, .. } => {
            let (p, rho) = conic_of(sketch, points, *curve)?;
            let m = (p[0] + p[2]) * 0.5;
            Measured::Shoulder {
                at: m + (p[1] - m) * rho,
                dir: (p[1] - m).normalize_or(DVec2::Y),
            }
        }
        _ => return None,
    })
}

/// `x` wrapped into (−π, π].
fn wrap_angle(x: f64) -> f64 {
    let w = (x + PI).rem_euclid(TAU) - PI;
    if w <= -PI {
        w + TAU
    } else {
        w
    }
}

/// The closest point of segment `ab` to `p` (plane coordinates).
fn closest_on_segment(p: DVec2, a: DVec2, b: DVec2) -> DVec2 {
    let ab = b - a;
    let t = if ab.length_squared() > 0.0 {
        ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0)
    } else {
        0.0
    };
    a + ab * t
}

/// An arc about `centre` of radius `r`, turning `sweep` from angle `start`,
/// widened to reach the direction of `reach` (a placed label) when given.
fn dimension_arc(
    centre: DVec2,
    r: f64,
    start: f64,
    sweep: f64,
    reach: Option<DVec2>,
) -> Vec<[DVec2; 2]> {
    let (mut lo, mut hi) = (sweep.min(0.0), sweep.max(0.0));
    if let Some(l) = reach {
        let rel = sweep * 0.5 + wrap_angle((l - centre).to_angle() - start - sweep * 0.5);
        lo = lo.min(rel);
        hi = hi.max(rel);
    }
    let n = ((hi - lo).abs() / (TAU / 96.0)).ceil().clamp(8.0, 96.0) as usize;
    let arc: Vec<DVec2> = (0..=n)
        .map(|i| centre + DVec2::from_angle(start + lo + (hi - lo) * i as f64 / n as f64) * r)
        .collect();
    arc.windows(2).map(|w| [w[0], w[1]]).collect()
}

/// An extension line along the ray from `vertex` in direction `e` out to
/// radius `r`, from where the segment `seg` stops reaching along it; `None`
/// when the segment already reaches the arc.
fn ray_extension(vertex: DVec2, e: DVec2, seg: (DVec2, DVec2), r: f64) -> Option<[DVec2; 2]> {
    let reach = (seg.0 - vertex)
        .dot(e)
        .max((seg.1 - vertex).dot(e))
        .max(0.0);
    (reach < r).then(|| [vertex + e * reach, vertex + e * r * 1.08])
}

/// A dimension's lines and label centre (plane coordinates). `placed` is
/// where the user put the label; without it the label goes outside, on the
/// side away from `centre` (the sketch's centroid), `off` from what is
/// measured.
fn dimension_geometry(
    m: Measured,
    placed: Option<DVec2>,
    off: f64,
    centre: DVec2,
) -> (Vec<[DVec2; 2]>, DVec2) {
    let tiny = 1e-9 * (1.0 + centre.length());
    let leader = |a: DVec2, b: DVec2, l: DVec2| {
        let q = closest_on_segment(l, a, b);
        ((l - q).length() > tiny).then_some([q, l])
    };
    match (m, placed) {
        (Measured::Aligned(a, b), None) => {
            let mut n = perp((b - a).normalize_or_zero());
            let away = ((a + b) * 0.5 - centre).dot(n);
            if away.abs() > tiny {
                if away < 0.0 {
                    n = -n;
                }
            } else if n.y < 0.0 || (n.y == 0.0 && n.x > 0.0) {
                n = -n;
            }
            let (a2, b2) = (a + n * off, b + n * off);
            (
                vec![[a, a + n * off * 1.25], [b, b + n * off * 1.25], [a2, b2]],
                // The label sits on the dimension line, covering its middle.
                (a2 + b2) * 0.5,
            )
        }
        (Measured::Aligned(a, b), Some(l)) => {
            // The dimension line runs parallel to the measured segment
            // through the label, reaching out to it past either end.
            let u = (b - a).normalize_or(DVec2::X);
            let n = perp(u);
            let h = (l - a).dot(n);
            let t = (l - a).dot(u);
            let len = (b - a).dot(u);
            let over = n * h.signum() * off * 0.25;
            let (a2, b2) = (a + n * h, b + n * h);
            (
                vec![
                    [a, a2 + over],
                    [b, b2 + over],
                    [a2 + u * t.min(0.0), a2 + u * t.max(len)],
                ],
                l,
            )
        }
        (Measured::Horizontal(pa, pb), placement) => {
            let y = match placement {
                Some(l) => l.y,
                None => {
                    let below = (pa.y + pb.y) * 0.5 < centre.y - tiny;
                    if below {
                        pa.y.min(pb.y) - off
                    } else {
                        pa.y.max(pb.y) + off
                    }
                }
            };
            let (a2, b2) = (DVec2::new(pa.x, y), DVec2::new(pb.x, y));
            let over = |p: DVec2| DVec2::new(0.0, (y - p.y).signum() * off * 0.25);
            let (x0, x1) = match placement {
                Some(l) => (pa.x.min(pb.x).min(l.x), pa.x.max(pb.x).max(l.x)),
                None => (pa.x.min(pb.x), pa.x.max(pb.x)),
            };
            (
                vec![
                    [pa, a2 + over(pa)],
                    [pb, b2 + over(pb)],
                    [DVec2::new(x0, y), DVec2::new(x1, y)],
                ],
                placement.unwrap_or((a2 + b2) * 0.5),
            )
        }
        (Measured::Vertical(pa, pb), placement) => {
            let x = match placement {
                Some(l) => l.x,
                None => {
                    let left = (pa.x + pb.x) * 0.5 < centre.x - tiny;
                    if left {
                        pa.x.min(pb.x) - off
                    } else {
                        pa.x.max(pb.x) + off
                    }
                }
            };
            let (a2, b2) = (DVec2::new(x, pa.y), DVec2::new(x, pb.y));
            let over = |p: DVec2| DVec2::new((x - p.x).signum() * off * 0.25, 0.0);
            let (y0, y1) = match placement {
                Some(l) => (pa.y.min(pb.y).min(l.y), pa.y.max(pb.y).max(l.y)),
                None => (pa.y.min(pb.y), pa.y.max(pb.y)),
            };
            (
                vec![
                    [pa, a2 + over(pa)],
                    [pb, b2 + over(pb)],
                    [DVec2::new(x, y0), DVec2::new(x, y1)],
                ],
                placement.unwrap_or((a2 + b2) * 0.5),
            )
        }
        (Measured::Foot { p, foot, near }, placement) => {
            let mut lines = vec![[p, foot], [foot, near]];
            let at = match placement {
                Some(l) => {
                    lines.extend(leader(p, foot, l));
                    l
                }
                None => {
                    let n = perp((foot - p).normalize_or_zero());
                    (p + foot) * 0.5 + n * off * 0.6
                }
            };
            (lines, at)
        }
        (
            Measured::Radial {
                centre: c,
                radius: r,
                dir,
                diameter,
            },
            placement,
        ) => {
            // Along the direction of the label, out to it when it is
            // outside the circle.
            let (dir, out, at) = match placement {
                Some(l) => {
                    let dir = (l - c).try_normalize().unwrap_or(dir);
                    (dir, (l - c).length().max(r), l)
                }
                None => (dir, r + off, c + dir * (r + off * 1.7)),
            };
            let from = if diameter { c - dir * r } else { c };
            (vec![[from, c + dir * out]], at)
        }
        (
            Measured::Angle {
                vertex,
                start,
                sweep,
                a,
                b,
            },
            placement,
        ) => {
            // A placed label picks which of the two opposite angles (equal
            // to each other) the arc spans: the one on its side.
            let flip = placement
                .is_some_and(|l| (l - vertex).dot(DVec2::from_angle(start + sweep * 0.5)) < 0.0);
            let start = if flip { start + PI } else { start };
            let r = match placement {
                Some(l) => (l - vertex).length(),
                None => off * 1.6,
            };
            let mut lines = dimension_arc(vertex, r, start, sweep, placement);
            let (ea, eb) = (DVec2::from_angle(start), DVec2::from_angle(start + sweep));
            lines.extend(ray_extension(vertex, ea, a, r));
            lines.extend(ray_extension(vertex, eb, b, r));
            let at = placement.unwrap_or_else(|| {
                vertex + DVec2::from_angle(start + sweep * 0.5) * (r + off * 0.9)
            });
            (lines, at)
        }
        (Measured::Across { from, to }, placement) => {
            let mut lines = vec![[from, to]];
            let at = match placement {
                Some(l) => {
                    lines.extend(leader(from, to, l));
                    l
                }
                None => to + (to - from).normalize_or_zero() * off * 0.8,
            };
            (lines, at)
        }
        (
            Measured::Turn {
                centre: c,
                source,
                sweep,
            },
            placement,
        ) => {
            // An arc about the centre from the source to its first copy,
            // with the two radii; a placed label sets the arc's radius.
            let v = source - c;
            let start = if v.length() > off { v.to_angle() } else { 0.0 };
            let r = match placement {
                Some(l) => (l - c).length(),
                None if v.length() > off => v.length(),
                None => off * 1.6,
            };
            let mut lines = dimension_arc(c, r, start, sweep, placement);
            lines.push([c, c + DVec2::from_angle(start) * r]);
            lines.push([c, c + DVec2::from_angle(start + sweep) * r]);
            let at = placement
                .unwrap_or_else(|| c + DVec2::from_angle(start + sweep * 0.5) * (r + off * 0.9));
            (lines, at)
        }
        (Measured::Shoulder { at, dir }, placement) => {
            // A leader from the shoulder point to the label, which goes on
            // along the median, outside the curve.
            let label = placement.unwrap_or(at + dir * off * 1.4);
            let lines = if (label - at).length() > tiny {
                vec![[at, label]]
            } else {
                Vec::new()
            };
            (lines, label)
        }
    }
}

fn dimension_annotation(
    d: &SketchDisplay,
    points: &BTreeMap<SketchId, DVec2>,
    map: &ScreenMap,
    c: &ConstraintEntry<SketchId>,
) -> Option<Annotation> {
    let m = measured(&d.sketch, points, &c.constraint)?;
    // Automatic offsets are a fixed screen distance, scaled where the
    // dimension's first item is (an angle's vertex, between nearly
    // parallel lines, can be far off screen).
    let sketch = &d.sketch;
    let near = c.constraint.referenced().iter().find_map(|e| {
        points
            .get(e)
            .copied()
            .or_else(|| line_ends(sketch, points, *e).map(|(a, _)| a))
            .or_else(|| circular(sketch, points, *e).map(|(c, _, _)| c))
            .or_else(|| conic_of(sketch, points, *e).map(|(p, _)| p[0]))
    })?;
    let ppp = map.plane_per_point(map.to_screen(near)?);
    let anchor = m.anchor();
    let centre = if points.is_empty() {
        DVec2::ZERO
    } else {
        points.values().copied().sum::<DVec2>() / points.len() as f64
    };
    let placed = d.labels.get(&c.id).map(|o| anchor + *o);
    let (lines, at) = dimension_geometry(m, placed, 24.0 * ppp, centre);
    let text = dimension_text(&c.constraint);
    let center = map.to_screen(at)?;
    let lines = lines
        .into_iter()
        .filter_map(|[a, b]| Some([map.to_screen(a)?, map.to_screen(b)?]))
        .collect();
    Some(Annotation {
        constraint: c.id,
        rect: Rect::from_center_size(center, label_size(&text)),
        text,
        lines,
        dimension: true,
    })
}

// ---------------------------------------------------------------------------
// Picking and snapping
// ---------------------------------------------------------------------------

/// Something under the pointer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pick {
    /// A point or curve.
    Entity(SketchId),
    /// A constraint glyph or dimension label.
    Constraint(SketchId),
}

impl Pick {
    /// The picked item.
    pub fn id(self) -> SketchId {
        match self {
            Pick::Entity(id) | Pick::Constraint(id) => id,
        }
    }
}

/// Distance from `p` to segment `ab`, and the closest point.
fn segment_distance(p: Pos2, a: Pos2, b: Pos2) -> (f32, Pos2) {
    let ab = b - a;
    let l2 = ab.length_sq();
    let t = if l2 > 0.0 {
        ((p - a).dot(ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let q = a + ab * t;
    (p.distance(q), q)
}

/// The nearest point entity within `radius` of `pos`.
fn nearest_point(
    d: &SketchDisplay,
    map: &ScreenMap,
    pos: Pos2,
    radius: f32,
    exclude: &[SketchId],
) -> Option<(SketchId, DVec2)> {
    d.points()
        .filter(|(id, _)| !exclude.contains(id))
        .filter_map(|(id, at)| Some((id, at, map.to_screen(at)?.distance(pos))))
        .filter(|(_, _, dist)| *dist <= radius)
        .min_by(|a, b| a.2.total_cmp(&b.2))
        .map(|(id, at, _)| (id, at))
}

/// The nearest curve within `radius` of `pos`, with the closest screen
/// point on it.
fn nearest_curve(
    d: &SketchDisplay,
    map: &ScreenMap,
    pos: Pos2,
    radius: f32,
) -> Option<(SketchId, Pos2)> {
    nearest_curve_where(d, map, pos, radius, |_| true)
}

/// [`nearest_curve`] among the curves for which `keep` is true.
fn nearest_curve_where(
    d: &SketchDisplay,
    map: &ScreenMap,
    pos: Pos2,
    radius: f32,
    keep: impl Fn(&DisplayCurve) -> bool,
) -> Option<(SketchId, Pos2)> {
    let mut best: Option<(SketchId, Pos2, f32)> = None;
    for c in d.curves.iter().filter(|c| keep(c)) {
        for run in map.polyline(&c.polyline) {
            for w in run.windows(2) {
                let (dist, q) = segment_distance(pos, w[0], w[1]);
                if dist <= radius && best.is_none_or(|b| dist < b.2) {
                    best = Some((c.id, q, dist));
                }
            }
        }
    }
    best.map(|(id, q, _)| (id, q))
}

/// What is under `pos`: labels first, then points, then curves.
pub fn pick(
    d: &SketchDisplay,
    annotations: &[Annotation],
    map: &ScreenMap,
    pos: Pos2,
) -> Option<Pick> {
    if let Some(a) = annotations.iter().rev().find(|a| a.rect.contains(pos)) {
        return Some(Pick::Constraint(a.constraint));
    }
    if let Some((id, _)) = nearest_point(d, map, pos, POINT_PICK, &[]) {
        return Some(Pick::Entity(id));
    }
    nearest_curve(d, map, pos, CURVE_PICK).map(|(id, _)| Pick::Entity(id))
}

/// What a click snapped to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SnapTarget {
    /// Nothing: a free position.
    Free,
    /// An existing point (the new point is made coincident with it).
    Point(SketchId),
    /// An existing curve (the new point is constrained onto it).
    Curve(SketchId),
    /// The sketch origin (the new point is fixed there).
    Origin,
}

/// A snapped pointer position.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Snap {
    /// Position in plane coordinates.
    pub at: DVec2,
    /// What it snapped to.
    pub target: SnapTarget,
}

/// The exact point of a curve nearest `p` for lines, circulars, ellipses and
/// conics, else `p`.
fn onto_curve(d: &SketchDisplay, id: SketchId, p: DVec2) -> DVec2 {
    match d.curve(id).map(|c| &c.curve) {
        Some(ResolvedCurve::Line { start, end }) => {
            let dir = (*end - *start).normalize_or_zero();
            *start + dir * (p - *start).dot(dir)
        }
        Some(ResolvedCurve::Circle { center, radius })
        | Some(ResolvedCurve::Arc { center, radius, .. }) => {
            let v = p - *center;
            if v.length() > 0.0 {
                *center + v.normalize() * *radius
            } else {
                p
            }
        }
        Some(ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        })
        | Some(ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            ..
        }) => closest_on_ellipse(*center, *major, *minor, p),
        Some(&ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        }) => {
            let c = [start, control, end];
            conic::point(c, weight, conic::closest_param(c, weight, p))
        }
        _ => p,
    }
}

/// The point of an ellipse closest to `p`, for a `p` near it (a point of
/// its drawn polyline): Newton's method on the distance's derivative from
/// `p`'s eccentric anomaly, which is already close.
fn closest_on_ellipse(c: DVec2, a: DVec2, b: DVec2, p: DVec2) -> DVec2 {
    let mut t = shapes::eccentric_anomaly(c, a, b, p);
    for _ in 0..8 {
        let (sin, cos) = t.sin_cos();
        let e = c + a * cos + b * sin - p;
        let d1 = b * cos - a * sin;
        let d2 = -(a * cos + b * sin);
        let (f, fp) = (d1.dot(e), d2.dot(e) + d1.length_squared());
        if fp <= 0.0 || fp.is_nan() {
            break;
        }
        let step = f / fp;
        t -= step;
        if step.abs() <= 1e-15 {
            break;
        }
    }
    shapes::ellipse_point(c, a, b, t)
}

/// Snaps `pos` to a point, the origin or a curve of `d`.
pub fn snap(d: &SketchDisplay, map: &ScreenMap, pos: Pos2, exclude: &[SketchId]) -> Option<Snap> {
    let free = map.to_plane(pos)?;
    if let Some((id, at)) = nearest_point(d, map, pos, POINT_PICK, exclude) {
        return Some(Snap {
            at,
            target: SnapTarget::Point(id),
        });
    }
    if map
        .to_screen(DVec2::ZERO)
        .is_some_and(|o| o.distance(pos) <= ORIGIN_PICK)
    {
        return Some(Snap {
            at: DVec2::ZERO,
            target: SnapTarget::Origin,
        });
    }
    // Not onto a text's letters: a point cannot be held on them (ADR 0015).
    let texts = text::texts(&d.sketch);
    if let Some((id, q)) = nearest_curve_where(d, map, pos, CURVE_PICK, |c| !texts.contains(&c.id))
    {
        let at = onto_curve(d, id, map.to_plane(q).unwrap_or(free));
        return Some(Snap {
            at,
            target: SnapTarget::Curve(id),
        });
    }
    Some(Snap {
        at: free,
        target: SnapTarget::Free,
    })
}

// ---------------------------------------------------------------------------
// Building edit batches
// ---------------------------------------------------------------------------

/// True when two clicks are different places (not the same point, nor
/// snapped to the same point entity).
fn distinct_clicks(a: &Snap, b: &Snap) -> bool {
    (a.at - b.at).length() > SAME_POINT
        && !(a.target == b.target && matches!(a.target, SnapTarget::Point(_)))
}

/// The degree of a control-point spline the tool makes with `count`
/// control points: cubic, or lower when there are too few points.
pub fn control_degree(count: usize) -> u32 {
    (count.saturating_sub(1)).clamp(1, 3) as u32
}

/// Chord-length parameters for a through-point spline through `points`,
/// from 0 to exactly 1; `None` when the points are all in one place.
pub fn chord_params(points: &[DVec2]) -> Option<Vec<f64>> {
    let mut lengths = vec![0.0];
    for w in points.windows(2) {
        let last = *lengths.last().unwrap_or(&0.0);
        lengths.push(last + (w[1] - w[0]).length());
    }
    let total = *lengths.last().unwrap_or(&0.0);
    if total <= SAME_POINT {
        return None;
    }
    let mut params: Vec<f64> = lengths.iter().map(|l| l / total).collect();
    if let Some(last) = params.last_mut() {
        *last = 1.0;
    }
    Some(params)
}

/// The smooth closed spline (ADR 0016) the spline tools make when the
/// clicks at `points` close on the first one again (`points` ends on its
/// first point): control-point (`control`) or through-point.
pub fn closed_spline(points: &[DVec2], control: bool) -> Option<ResolvedCurve> {
    let n = points.len().checked_sub(1)?;
    if control {
        openshape_sketch_solver::control_spline_curve(points, control_degree(n), true)
    } else {
        let params = chord_params(points)?;
        openshape_sketch_solver::interpolated_spline_curve(points, &params, true)
    }
}

/// The control-point spline (clamped, uniform knots) on `control`, from
/// three points; `None` with fewer.
pub fn control_spline(control: &[DVec2]) -> Option<ProfileCurve> {
    if control.len() < 3 {
        return None;
    }
    let degree = control_degree(control.len()) as usize;
    Some(ProfileCurve::Spline {
        degree,
        knots: openshape_sketch_solver::clamped_uniform_knots(control.len(), degree),
        points: control.to_vec(),
        weights: None,
    })
}

/// Where a tangent arc clicked at `s` starts: `s` must have snapped to the
/// end point of a line, arc or spline.
fn tangent_start(d: &SketchDisplay, s: &Snap) -> Option<TangentStart> {
    let SnapTarget::Point(point) = s.target else {
        return None;
    };
    let (curve, at_end) = curve_ending_at(&d.sketch, point)?;
    let direction = leaving_direction(&d.curve(curve)?.curve, at_end)?;
    Some(TangentStart {
        curve,
        point,
        at: d.sketch.point(point)?,
        direction,
    })
}

/// Centre and radius of the circle through three points.
pub fn circumcircle(a: DVec2, b: DVec2, c: DVec2) -> Option<(DVec2, f64)> {
    let (b, c) = (b - a, c - a);
    let d = 2.0 * b.perp_dot(c);
    // Reject (nearly) collinear points: the circle would be enormous.
    if d.abs() <= 2e-4 * b.length() * c.length() || d == 0.0 {
        return None;
    }
    let u = DVec2::new(
        c.y * b.length_squared() - b.y * c.length_squared(),
        b.x * c.length_squared() - c.x * b.length_squared(),
    ) / d;
    Some((a + u, u.length()))
}

/// Points of an arc, `sweep` signed.
fn arc_points(center: DVec2, radius: f64, start: f64, sweep: f64) -> Vec<DVec2> {
    let n = ((sweep.abs() / (3.75f64).to_radians()).ceil() as usize).clamp(4, 96);
    (0..=n)
        .map(|i| center + DVec2::from_angle(start + sweep * i as f64 / n as f64) * radius)
        .collect()
}

/// Counter-clockwise sweep from angle `a` to angle `b`, in [0, 2π).
fn ccw(a: f64, b: f64) -> f64 {
    (b - a).rem_euclid(TAU)
}

// ---------------------------------------------------------------------------
// Tools
// ---------------------------------------------------------------------------

/// A sketch tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Tool {
    /// Select, drag, delete.
    Select,
    /// Lines, chained end to end.
    Line,
    /// Axis-aligned rectangle from two corners.
    Rectangle,
    /// Axis-aligned rectangle from its centre and a corner.
    CenterRectangle,
    /// Rectangle at any angle: two corners along one side, then the width.
    AlignedRectangle,
    /// Circle from its centre and a point on it.
    Circle,
    /// Arc through three points (start, end, a point on the arc).
    Arc,
    /// Arc from its centre, start and end.
    CenterArc,
    /// Arc continuing a curve tangentially from its end.
    TangentArc,
    /// Slot: two arc centres, then the width.
    Slot,
    /// Regular polygon from its centre and a vertex.
    Polygon,
    /// A point.
    Point,
    /// Spline through clicked points (finish with Enter or a double-click).
    Spline,
    /// Removes the piece of a curve between the curves crossing it.
    Trim,
    /// Extends a line or arc to the next curve.
    Extend,
    /// Splits a line or arc in two at a point.
    Split,
    /// Rounds the corner where two lines or arcs meet.
    SketchFillet,
    /// Mirror images of curves and points about a line.
    Mirror,
    /// Offsets a chain of lines and arcs, or a circle, to one side.
    Offset,
    /// Copies of the selection along a line.
    LinearPattern,
    /// Copies of the selection about a centre.
    CircularPattern,
    /// A B-spline shaped by clicked control points.
    ControlSpline,
    /// Ellipse from its centre, the end of one axis and the other axis's
    /// extent.
    Ellipse,
    /// Elliptical arc: the ellipse's three clicks, then its start and end.
    EllipticalArc,
    /// Conic arc: its start, its end, where its end tangents meet, then
    /// its shoulder (design note 0009 §7).
    Conic,
    /// Moves the selection (or a copy) from one point to another (design
    /// note 0007 §13).
    Move,
    /// Turns the selection (or a copy) about a centre.
    Rotate,
    /// Scales the selection (or a copy) about a base point.
    Scale,
    /// A line of text in a box: a corner, then the height (ADR 0015).
    Text,
}

impl Tool {
    /// Every tool, in toolbar order.
    pub const ALL: [Tool; 29] = [
        Tool::Select,
        Tool::Line,
        Tool::Rectangle,
        Tool::CenterRectangle,
        Tool::AlignedRectangle,
        Tool::Circle,
        Tool::Ellipse,
        Tool::Arc,
        Tool::CenterArc,
        Tool::TangentArc,
        Tool::EllipticalArc,
        Tool::Slot,
        Tool::Polygon,
        Tool::Point,
        Tool::Spline,
        Tool::ControlSpline,
        Tool::Conic,
        Tool::Text,
        Tool::Trim,
        Tool::Extend,
        Tool::Split,
        Tool::SketchFillet,
        Tool::Mirror,
        Tool::Offset,
        Tool::LinearPattern,
        Tool::CircularPattern,
        Tool::Move,
        Tool::Rotate,
        Tool::Scale,
    ];

    /// The tools sharing one toolbar button with this one (shown in its
    /// drop-down), in menu order; the first is the group's default.
    pub fn group(self) -> &'static [Tool] {
        match self {
            Tool::Rectangle | Tool::CenterRectangle | Tool::AlignedRectangle => &[
                Tool::Rectangle,
                Tool::CenterRectangle,
                Tool::AlignedRectangle,
            ],
            Tool::Arc | Tool::CenterArc | Tool::TangentArc | Tool::EllipticalArc => &[
                Tool::Arc,
                Tool::CenterArc,
                Tool::TangentArc,
                Tool::EllipticalArc,
            ],
            Tool::Trim | Tool::Extend | Tool::Split => &[Tool::Trim, Tool::Extend, Tool::Split],
            Tool::Select => &[Tool::Select],
            Tool::Line => &[Tool::Line],
            Tool::Circle | Tool::Ellipse => &[Tool::Circle, Tool::Ellipse],
            Tool::Slot => &[Tool::Slot],
            Tool::Polygon => &[Tool::Polygon],
            Tool::Point => &[Tool::Point],
            Tool::Spline | Tool::ControlSpline | Tool::Conic => {
                &[Tool::Spline, Tool::ControlSpline, Tool::Conic]
            }
            Tool::SketchFillet => &[Tool::SketchFillet],
            Tool::Mirror => &[Tool::Mirror],
            Tool::Offset => &[Tool::Offset],
            Tool::LinearPattern | Tool::CircularPattern => {
                &[Tool::LinearPattern, Tool::CircularPattern]
            }
            Tool::Move | Tool::Rotate | Tool::Scale => &[Tool::Move, Tool::Rotate, Tool::Scale],
            Tool::Text => &[Tool::Text],
        }
    }

    /// Toolbar label.
    pub fn label(self) -> &'static str {
        match self {
            Tool::Select => "Select",
            Tool::Line => "Line",
            Tool::Rectangle => "Corner rectangle",
            Tool::CenterRectangle => "Centre-point rectangle",
            Tool::AlignedRectangle => "Aligned rectangle",
            Tool::Circle => "Circle",
            Tool::Arc => "3-point arc",
            Tool::CenterArc => "Centre arc",
            Tool::TangentArc => "Tangent arc",
            Tool::Slot => "Slot",
            Tool::Polygon => "Polygon",
            Tool::Point => "Point",
            Tool::Spline => "Spline",
            Tool::ControlSpline => "Control-point spline",
            Tool::Trim => "Trim",
            Tool::Extend => "Extend",
            Tool::Split => "Split",
            Tool::SketchFillet => "Sketch fillet",
            Tool::Mirror => "Mirror",
            Tool::Offset => "Offset",
            Tool::LinearPattern => "Linear pattern",
            Tool::CircularPattern => "Circular pattern",
            Tool::Ellipse => "Ellipse",
            Tool::EllipticalArc => "Elliptical arc",
            Tool::Conic => "Conic",
            Tool::Move => "Move",
            Tool::Rotate => "Rotate",
            Tool::Scale => "Scale",
            Tool::Text => "Text",
        }
    }

    /// True for the tools that edit existing curves (trim, extend, split,
    /// fillet).
    pub fn edits_curves(self) -> bool {
        matches!(
            self,
            Tool::Trim | Tool::Extend | Tool::Split | Tool::SketchFillet
        )
    }

    /// True for the pattern tools.
    pub fn is_pattern(self) -> bool {
        matches!(self, Tool::LinearPattern | Tool::CircularPattern)
    }

    /// True for the transform tools (move, rotate, scale).
    pub fn is_transform(self) -> bool {
        matches!(self, Tool::Move | Tool::Rotate | Tool::Scale)
    }

    /// Keyboard shortcut.
    pub fn key(self) -> Option<egui::Key> {
        Some(match self {
            Tool::Line => egui::Key::L,
            Tool::Rectangle => egui::Key::R,
            Tool::Circle => egui::Key::C,
            Tool::Arc => egui::Key::A,
            Tool::Point => egui::Key::P,
            Tool::Spline => egui::Key::S,
            Tool::Trim => egui::Key::T,
            Tool::Select
            | Tool::CenterRectangle
            | Tool::AlignedRectangle
            | Tool::CenterArc
            | Tool::TangentArc
            | Tool::Slot
            | Tool::Polygon
            | Tool::Extend
            | Tool::Split
            | Tool::SketchFillet
            | Tool::Mirror
            | Tool::Offset
            | Tool::LinearPattern
            | Tool::CircularPattern
            | Tool::ControlSpline
            | Tool::Ellipse
            | Tool::EllipticalArc
            | Tool::Conic
            | Tool::Move
            | Tool::Rotate
            | Tool::Scale
            | Tool::Text => return None,
        })
    }

    /// What to do next, given the clicks so far.
    pub fn hint(self, clicks: usize) -> &'static str {
        match (self, clicks) {
            (Tool::Select, _) => {
                "Click to select (Shift adds), drag to move, box-select or place a dimension label, Delete removes"
            }
            (Tool::Line, 0) => "Click the start point",
            (Tool::Line, _) => "Click the end point (double-click or Esc ends the chain)",
            (Tool::Rectangle, 0) => "Click the first corner",
            (Tool::Rectangle, _) => "Click the opposite corner",
            (Tool::CenterRectangle, 0) => "Click the centre",
            (Tool::CenterRectangle, _) => "Click a corner",
            (Tool::AlignedRectangle, 0) => "Click the first corner",
            (Tool::AlignedRectangle, 1) => "Click the next corner along one side",
            (Tool::AlignedRectangle, _) => "Click to set the width",
            (Tool::Circle, 0) => "Click the centre",
            (Tool::Circle, _) => "Click a point on the circle",
            (Tool::Arc, 0) => "Click the start point",
            (Tool::Arc, 1) => "Click the end point",
            (Tool::Arc, _) => "Click a point on the arc",
            (Tool::CenterArc, 0) => "Click the centre",
            (Tool::CenterArc, 1) => "Click the start point",
            (Tool::CenterArc, _) => "Move around the centre and click the end point",
            (Tool::TangentArc, 0) => "Click the end of a line, arc or spline",
            (Tool::TangentArc, _) => "Click the end of the arc",
            (Tool::Slot, 0) => "Click the first end's centre",
            (Tool::Slot, 1) => "Click the second end's centre",
            (Tool::Slot, _) => "Click to set the width",
            (Tool::Polygon, 0) => "Click the centre",
            (Tool::Polygon, _) => "Click a vertex",
            (Tool::Point, _) => "Click to place a point",
            (Tool::Spline, 0) => "Click the first point",
            (Tool::Spline, 1 | 2) => "Click more points; Enter or double-click finishes",
            (Tool::Spline, _) => {
                "Click more points; Enter or double-click finishes, the first point closes it smoothly"
            }
            (Tool::ControlSpline, 0) => "Click the first control point (the curve starts there)",
            (Tool::ControlSpline, 1 | 2) => "Click more control points",
            (Tool::ControlSpline, _) => {
                "Click more control points; Enter or double-click finishes, the first point closes it smoothly"
            }
            (Tool::Trim, _) => {
                "Click the piece of a curve to remove (a curve nothing crosses is deleted)"
            }
            (Tool::Extend, _) => {
                "Click near the end of a line or arc to extend it to the next curve"
            }
            (Tool::Split, _) => "Click a line or arc where it should be split",
            (Tool::SketchFillet, _) => "Click the corner point where two lines or arcs meet",
            (Tool::Mirror, _) => "Click the mirror line, then the curves to mirror",
            (Tool::Offset, _) => {
                "Click a line, arc or circle to offset it with the curves joined to it"
            }
            (Tool::LinearPattern, 0) => {
                "Select curves or points, then click where the step between copies starts"
            }
            (Tool::LinearPattern, _) => "Click where the step to the first copy ends",
            (Tool::CircularPattern, _) => "Select curves or points, then click the centre",
            (Tool::Move, _) => "Select curves or points, then click where to move them from",
            (Tool::Rotate, _) => "Select curves or points, then click the centre to turn them about",
            (Tool::Scale, _) => "Select curves or points, then click the base point to scale them about",
            (Tool::Ellipse | Tool::EllipticalArc, 0) => "Click the centre",
            (Tool::Ellipse | Tool::EllipticalArc, 1) => "Click the end of the first axis",
            (Tool::Ellipse, _) | (Tool::EllipticalArc, 2) => "Click to set the second axis",
            (Tool::EllipticalArc, 3) => "Click the start point on the ellipse",
            (Tool::EllipticalArc, _) => "Move around the ellipse and click the end point",
            (Tool::Text, 0) => "Type the text and choose a font, then click a corner of its box",
            (Tool::Text, _) => "Click to set the text's height (the box grows towards the pointer)",
            (Tool::Conic, 0) => "Click the start point",
            (Tool::Conic, 1) => "Click the end point",
            (Tool::Conic, 2) => "Click where the tangents at the two ends meet",
            (Tool::Conic, _) => {
                "Click the shoulder: how far the curve reaches towards the corner sets its ρ"
            }
        }
    }
}

/// A relation inferred for the line being drawn (design note 0007 §5),
/// added as a constraint with it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Inference {
    /// Along the sketch x axis.
    Horizontal,
    /// Along the sketch y axis.
    Vertical,
    /// Parallel to a line.
    Parallel(SketchId),
    /// Perpendicular to a line.
    Perpendicular(SketchId),
    /// Tangent to the curve whose end the line starts at.
    Tangent(SketchId),
}

impl Inference {
    /// The glyph shown next to the pointer (the constraint's glyph).
    pub fn glyph(self) -> &'static str {
        match self {
            Inference::Horizontal => "H",
            Inference::Vertical => "V",
            Inference::Parallel(_) => "‖",
            Inference::Perpendicular(_) => "⟂",
            Inference::Tangent(_) => "T",
        }
    }
}

/// A geometric constraint applied to the selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ConstraintTool {
    /// Two points coincide, or a point lies on a curve.
    Coincident,
    /// Lines horizontal, or two points level.
    Horizontal,
    /// Lines vertical, or two points aligned vertically.
    Vertical,
    /// Lines parallel.
    Parallel,
    /// Two lines perpendicular.
    Perpendicular,
    /// Two curves tangent.
    Tangent,
    /// A spline and a curve curvature continuous at a shared end.
    Curvature,
    /// Equal lengths or radii.
    Equal,
    /// A point at a line's midpoint.
    Midpoint,
    /// Two circles or arcs share a centre.
    Concentric,
    /// Two points mirror each other about a line.
    Symmetric,
    /// A line through a circle's centre.
    Normal,
    /// Points fixed where they are.
    Fix,
}

impl ConstraintTool {
    /// Every constraint tool, in toolbar order.
    pub const ALL: [ConstraintTool; 13] = [
        ConstraintTool::Coincident,
        ConstraintTool::Horizontal,
        ConstraintTool::Vertical,
        ConstraintTool::Parallel,
        ConstraintTool::Perpendicular,
        ConstraintTool::Tangent,
        ConstraintTool::Curvature,
        ConstraintTool::Equal,
        ConstraintTool::Midpoint,
        ConstraintTool::Concentric,
        ConstraintTool::Symmetric,
        ConstraintTool::Normal,
        ConstraintTool::Fix,
    ];

    /// Toolbar label.
    pub fn label(self) -> &'static str {
        match self {
            ConstraintTool::Coincident => "Coincident",
            ConstraintTool::Horizontal => "Horizontal",
            ConstraintTool::Vertical => "Vertical",
            ConstraintTool::Parallel => "Parallel",
            ConstraintTool::Perpendicular => "Perpendicular",
            ConstraintTool::Tangent => "Tangent",
            ConstraintTool::Curvature => "Curvature",
            ConstraintTool::Equal => "Equal",
            ConstraintTool::Midpoint => "Midpoint",
            ConstraintTool::Concentric => "Concentric",
            ConstraintTool::Symmetric => "Symmetric",
            ConstraintTool::Normal => "Normal",
            ConstraintTool::Fix => "Fix",
        }
    }

    /// What the tool needs selected.
    pub fn needs(self) -> &'static str {
        match self {
            ConstraintTool::Coincident => "two points, or a point and a curve",
            ConstraintTool::Horizontal | ConstraintTool::Vertical => {
                "one or more lines, or two points"
            }
            ConstraintTool::Parallel => "two or more lines",
            ConstraintTool::Perpendicular => "two lines",
            ConstraintTool::Tangent => "two curves (not two lines)",
            ConstraintTool::Curvature => "a spline and a curve sharing an end point",
            ConstraintTool::Equal => "two or more lines, or two or more circles and arcs",
            ConstraintTool::Midpoint => "a point and a line",
            ConstraintTool::Concentric => "two circles or arcs",
            ConstraintTool::Symmetric => "two points and a line",
            ConstraintTool::Normal => "a line and a circle or arc",
            ConstraintTool::Fix => "one or more points",
        }
    }
}

/// What kind of entity an id is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Point,
    Line,
    Circle,
    Arc,
    Spline,
    /// An ellipse or an elliptical arc.
    Ellipse,
    /// A conic arc.
    Conic,
    /// Text (ADR 0015).
    Text,
}

impl Kind {
    fn circular(self) -> bool {
        matches!(self, Kind::Circle | Kind::Arc)
    }
}

fn kind(sketch: &SketchData, id: SketchId) -> Option<Kind> {
    Some(match sketch.entity(id)?.geometry {
        Geometry::Point { .. } => Kind::Point,
        Geometry::Line { .. } => Kind::Line,
        Geometry::Circle { .. } => Kind::Circle,
        Geometry::Arc { .. } => Kind::Arc,
        Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. } => Kind::Spline,
        Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } => Kind::Ellipse,
        Geometry::Conic { .. } => Kind::Conic,
        Geometry::Text { .. } => Kind::Text,
    })
}

/// Why tangency and curvature continuity refuse a smooth closed spline.
const CLOSED_HAS_NO_ENDS: &str = "A smooth closed spline has no end points: tangency and \
     curvature continuity join curves at a shared end (ADR 0016)";

/// True when `id` is a smooth closed spline (ADR 0016).
fn is_closed_spline(sketch: &SketchData, id: SketchId) -> bool {
    sketch
        .entity(id)
        .is_some_and(|e| e.geometry.is_closed_spline())
}

/// The end points of a curve (none for a circle, an ellipse or a smooth
/// closed spline).
fn end_points(sketch: &SketchData, id: SketchId) -> Vec<SketchId> {
    sketch
        .entity(id)
        .and_then(|e| e.geometry.end_points())
        .map_or_else(Vec::new, |(s, e)| vec![s, e])
}

/// An end point of `a` that `b` shares: the same point, a coincident one,
/// or one at the same position.
fn shared_end(sketch: &SketchData, a: SketchId, b: SketchId) -> Option<SketchId> {
    let (ea, eb) = (end_points(sketch, a), end_points(sketch, b));
    let coincident = |p: SketchId, q: SketchId| {
        p == q
            || sketch.constraints.iter().any(|c| {
                matches!(c.constraint, Constraint::Coincident { a, b }
                    if (a == p && b == q) || (a == q && b == p))
            })
            || matches!((sketch.point(p), sketch.point(q)),
                (Some(x), Some(y)) if (x - y).length() <= SAME_POINT)
    };
    ea.iter()
        .find(|p| eb.iter().any(|q| coincident(**p, *q)))
        .copied()
}

/// The centre of a circle, arc, ellipse or elliptical arc.
fn centre_of(sketch: &SketchData, id: SketchId) -> Option<DVec2> {
    match sketch.entity(id)?.geometry {
        Geometry::Circle { center, .. }
        | Geometry::Arc { center, .. }
        | Geometry::Ellipse { center, .. }
        | Geometry::EllipticalArc { center, .. } => sketch.point(center),
        _ => None,
    }
}

fn center_radius(sketch: &SketchData, id: SketchId) -> Option<(DVec2, f64)> {
    match sketch.entity(id)?.geometry {
        Geometry::Circle { center, radius } => Some((sketch.point(center)?, radius)),
        Geometry::Arc { center, start, .. } => {
            let c = sketch.point(center)?;
            Some((c, (sketch.point(start)? - c).length()))
        }
        _ => None,
    }
}

/// The constraints a tool adds for a selection (entities in selection
/// order), or why it does not apply.
pub fn constraints_for(
    tool: ConstraintTool,
    selection: &[SketchId],
    sketch: &SketchData,
) -> Result<Vec<Constraint<SketchId>>, String> {
    if selection.iter().any(|id| text::is_text(sketch, *id)) {
        return Err(text::CONSTRAIN_THE_BOX.into());
    }
    let items: Vec<(SketchId, Kind)> = selection
        .iter()
        .filter_map(|id| Some((*id, kind(sketch, *id)?)))
        .collect();
    let kinds: Vec<Kind> = items.iter().map(|(_, k)| *k).collect();
    let ids: Vec<SketchId> = items.iter().map(|(id, _)| *id).collect();
    let all = |k: Kind| !kinds.is_empty() && kinds.iter().all(|x| *x == k);
    let find = |want: Kind| items.iter().find(|(_, k)| *k == want).map(|(id, _)| *id);
    let need = || Err(format!("{} needs {}", tool.label(), tool.needs()));
    use Constraint as C;
    Ok(match tool {
        ConstraintTool::Coincident => match kinds.as_slice() {
            [Kind::Point, Kind::Point] => vec![C::Coincident {
                a: ids[0],
                b: ids[1],
            }],
            [Kind::Point, _] => vec![C::PointOnCurve {
                point: ids[0],
                curve: ids[1],
            }],
            [_, Kind::Point] => vec![C::PointOnCurve {
                point: ids[1],
                curve: ids[0],
            }],
            _ => return need(),
        },
        ConstraintTool::Horizontal | ConstraintTool::Vertical => {
            let h = tool == ConstraintTool::Horizontal;
            if all(Kind::Line) {
                ids.iter()
                    .map(|&line| {
                        if h {
                            C::Horizontal { line }
                        } else {
                            C::Vertical { line }
                        }
                    })
                    .collect()
            } else if kinds == [Kind::Point, Kind::Point] {
                let (a, b) = (ids[0], ids[1]);
                vec![if h {
                    C::HorizontalPoints { a, b }
                } else {
                    C::VerticalPoints { a, b }
                }]
            } else {
                return need();
            }
        }
        ConstraintTool::Parallel if all(Kind::Line) && ids.len() >= 2 => ids[1..]
            .iter()
            .map(|&b| C::Parallel { a: ids[0], b })
            .collect(),
        ConstraintTool::Perpendicular if kinds == [Kind::Line, Kind::Line] => {
            vec![C::Perpendicular {
                a: ids[0],
                b: ids[1],
            }]
        }
        ConstraintTool::Equal
            if ids.len() >= 2 && (all(Kind::Line) || kinds.iter().all(|k| k.circular())) =>
        {
            ids[1..]
                .iter()
                .map(|&b| C::Equal { a: ids[0], b })
                .collect()
        }
        ConstraintTool::Midpoint if ids.len() == 2 => match (find(Kind::Point), find(Kind::Line)) {
            (Some(point), Some(line)) => vec![C::Midpoint { point, line }],
            _ => return need(),
        },
        ConstraintTool::Concentric if ids.len() == 2 && kinds.iter().all(|k| k.circular()) => {
            vec![C::Concentric {
                a: ids[0],
                b: ids[1],
            }]
        }
        ConstraintTool::Symmetric if ids.len() == 3 => {
            let pts: Vec<SketchId> = items
                .iter()
                .filter(|(_, k)| *k == Kind::Point)
                .map(|(id, _)| *id)
                .collect();
            match (pts.as_slice(), find(Kind::Line)) {
                ([a, b], Some(axis)) => vec![C::Symmetric { a: *a, b: *b, axis }],
                _ => return need(),
            }
        }
        ConstraintTool::Normal if ids.len() == 2 => {
            let curve = items.iter().find(|(_, k)| k.circular()).map(|(id, _)| *id);
            match (find(Kind::Line), curve) {
                (Some(line), Some(curve)) => vec![C::Normal { line, curve }],
                _ => return need(),
            }
        }
        ConstraintTool::Fix if all(Kind::Point) => ids
            .iter()
            .filter_map(|&point| {
                Some(C::Fix {
                    point,
                    at: sketch.point(point)?,
                })
            })
            .collect(),
        ConstraintTool::Tangent if ids.len() == 2 && !kinds.contains(&Kind::Point) => {
            let (a, b) = (ids[0], ids[1]);
            if kinds == [Kind::Line, Kind::Line] {
                return Err("Two lines can't be tangent; use Parallel".into());
            }
            if is_closed_spline(sketch, a) || is_closed_spline(sketch, b) {
                return Err(CLOSED_HAS_NO_ENDS.into());
            }
            let kind = if let Some(point) = shared_end(sketch, a, b) {
                TangentKind::AtPoint { point }
            } else if kinds.contains(&Kind::Spline) {
                return Err("Tangency with a spline needs a shared end point".into());
            } else if kinds.contains(&Kind::Conic) {
                return Err("Tangency with a conic needs a shared end point".into());
            } else if kinds.contains(&Kind::Line) {
                let (line, circ) = if kinds[0] == Kind::Line {
                    (a, b)
                } else {
                    (b, a)
                };
                let ends = end_points(sketch, line);
                let (s, e) = (
                    sketch.point(ends[0]).ok_or("missing point")?,
                    sketch.point(ends[1]).ok_or("missing point")?,
                );
                let c = centre_of(sketch, circ).ok_or("missing centre")?;
                TangentKind::LineCircle {
                    center_left: (e - s).perp_dot(c - s) > 0.0,
                }
            } else if kinds.contains(&Kind::Ellipse) {
                return Err(
                    "An ellipse can be tangent to a line, or to a curve at a shared end point"
                        .into(),
                );
            } else {
                let (c1, r1) = center_radius(sketch, a).ok_or("missing centre")?;
                let (c2, r2) = center_radius(sketch, b).ok_or("missing centre")?;
                TangentKind::CircleCircle {
                    internal: (c1 - c2).length() < r1.max(r2),
                }
            };
            vec![C::Tangent { a, b, kind }]
        }
        ConstraintTool::Curvature if ids.len() == 2 && !kinds.contains(&Kind::Point) => {
            let (a, b) = (ids[0], ids[1]);
            if !kinds.contains(&Kind::Spline) {
                return Err("Curvature continuity needs a spline on at least one side".into());
            }
            if is_closed_spline(sketch, a) || is_closed_spline(sketch, b) {
                return Err(CLOSED_HAS_NO_ENDS.into());
            }
            let ends = |id: SketchId| !end_points(sketch, id).is_empty();
            if !ends(a) || !ends(b) {
                return Err("A full circle or ellipse has no end to join".into());
            }
            // Degree 1: a control spline of degree 1, or an interpolated
            // one through two points (degree min(3, points − 1)).
            let straight = |id: SketchId| match sketch.entity(id).map(|e| &e.geometry) {
                Some(Geometry::ControlSpline { degree, .. }) => *degree == 1,
                Some(Geometry::InterpolatedSpline { through, .. }) => through.len() == 2,
                _ => false,
            };
            if straight(a) || straight(b) {
                return Err("A degree-1 spline is straight at its ends; use Tangent".into());
            }
            let Some(point) = shared_end(sketch, a, b) else {
                return Err("Curvature continuity needs the curves to share an end point".into());
            };
            vec![C::Curvature { a, b, point }]
        }
        _ => return need(),
    })
}

/// The tangencies at a shared point that new curvature constraints `cs`
/// make redundant: on the same two curves at the same or a coincident
/// point.
fn replaced_tangencies(sketch: &SketchData, cs: &[Constraint<SketchId>]) -> Vec<SketchId> {
    let same = |p: SketchId, q: SketchId| {
        p == q
            || sketch.constraints.iter().any(|c| {
                matches!(c.constraint, Constraint::Coincident { a, b }
                    if (a == p && b == q) || (a == q && b == p))
            })
    };
    sketch
        .constraints
        .iter()
        .filter(|e| {
            let Constraint::Tangent {
                a,
                b,
                kind: TangentKind::AtPoint { point },
            } = e.constraint
            else {
                return false;
            };
            cs.iter().any(|c| {
                matches!(*c, Constraint::Curvature { a: x, b: y, point: q }
                    if ((x == a && y == b) || (x == b && y == a)) && same(point, q))
            })
        })
        .map(|e| e.id)
        .collect()
}

/// The dimension a selection gets (value not yet measured), or why none.
pub fn dimension_for(
    selection: &[SketchId],
    sketch: &SketchData,
) -> Result<Constraint<SketchId>, String> {
    if selection.iter().any(|id| text::is_text(sketch, *id)) {
        return Err(text::CONSTRAIN_THE_BOX.into());
    }
    let items: Vec<(SketchId, Kind)> = selection
        .iter()
        .filter_map(|id| Some((*id, kind(sketch, *id)?)))
        .collect();
    use Constraint as C;
    Ok(match items.as_slice() {
        [(line, Kind::Line)] => C::Length {
            line: *line,
            value: 0.0,
        },
        [(curve, Kind::Circle)] => C::Diameter {
            curve: *curve,
            value: 0.0,
        },
        [(curve, Kind::Arc)] => C::Radius {
            curve: *curve,
            value: 0.0,
        },
        [(curve, Kind::Conic)] => C::Rho {
            curve: *curve,
            value: 0.5,
        },
        [(a, Kind::Point), (b, Kind::Point)] => C::Distance {
            a: *a,
            b: *b,
            kind: DistanceKind::Aligned,
            value: 0.0,
        },
        [(a, Kind::Point), (b, Kind::Line)] | [(b, Kind::Line), (a, Kind::Point)] => C::Distance {
            a: *a,
            b: *b,
            kind: DistanceKind::Aligned,
            value: 0.0,
        },
        [(a, Kind::Line), (b, Kind::Line)] => {
            let dir = |l: SketchId| -> Option<DVec2> {
                let e = end_points(sketch, l);
                Some(sketch.point(e[1])? - sketch.point(e[0])?)
            };
            let (da, db) = (
                dir(*a).ok_or("missing point")?,
                dir(*b).ok_or("missing point")?,
            );
            if da.perp_dot(db).abs() <= 1e-9 * da.length() * db.length() {
                C::Distance {
                    a: *a,
                    b: *b,
                    kind: DistanceKind::Aligned,
                    value: 0.0,
                }
            } else {
                C::Angle {
                    a: *a,
                    b: *b,
                    value: 0.0,
                }
            }
        }
        _ => {
            return Err(
                "Dimension needs a line, a circle or arc, a conic, two points, a point \
                        and a line, or two lines"
                    .into(),
            )
        }
    })
}

/// A dimension value being entered.
#[derive(Clone, Debug, PartialEq)]
pub struct DimensionPrompt {
    /// The dimension: new (not yet added) or existing.
    pub constraint: Constraint<SketchId>,
    /// The existing dimension being changed, if any.
    pub existing: Option<SketchId>,
    /// The value as typed (bare numbers are millimetres or degrees).
    pub text: String,
    /// Where to show the prompt.
    pub at: Pos2,
    /// Why the last entry was refused.
    pub error: Option<String>,
}

impl DimensionPrompt {
    /// True for angles.
    pub fn is_angle(&self) -> bool {
        matches!(
            self.constraint,
            Constraint::Angle { .. } | Constraint::CircularPattern { .. }
        )
    }

    /// True for a conic's ρ, a plain number.
    pub fn is_rho(&self) -> bool {
        matches!(self.constraint, Constraint::Rho { .. })
    }

    /// What the typed value is in, shown beside it: millimetres, degrees,
    /// or ρ.
    pub fn unit(&self) -> &'static str {
        if self.is_angle() {
            "°"
        } else if self.is_rho() {
            "ρ"
        } else {
            "mm"
        }
    }

    /// How a point-to-point distance is measured, if it is one.
    pub fn distance_kind(&self) -> Option<DistanceKind> {
        match self.constraint {
            Constraint::Distance { kind, .. } if self.existing.is_none() => Some(kind),
            _ => None,
        }
    }

    /// Parses the typed value into SI units, keeping an angle's sign.
    fn value(&self) -> Result<f64, String> {
        let current = self.constraint.value().unwrap_or(0.0);
        if self.is_rho() {
            let v: f64 = self
                .text
                .trim()
                .parse()
                .map_err(|_| "ρ is a number between 0 and 1".to_string())?;
            return if (RHO_MIN..=1.0 - RHO_MIN).contains(&v) {
                Ok(v)
            } else {
                Err(format!(
                    "ρ is between 0 and 1 ({RHO_MIN} to {})",
                    1.0 - RHO_MIN
                ))
            };
        }
        if self.is_angle() {
            let v = Quantity::parse_angle(&self.text).map_err(|e| e.to_string())?;
            let v = v.abs() * if current < 0.0 { -1.0 } else { 1.0 };
            if let Constraint::CircularPattern { copies, .. } = &self.constraint {
                return pattern_angle(v, copies);
            }
            if v.abs() > PI {
                return Err("an angle between lines is at most 180°".into());
            }
            Ok(v)
        } else {
            let v = Quantity::parse_length(&self.text, LengthUnit::Millimetre)
                .map_err(|e| e.to_string())?;
            if v <= 0.0 {
                return Err("a length must be positive".into());
            }
            Ok(v)
        }
    }
}

/// Checks a circular pattern's typed angle: more than nothing and less
/// than a turn, and no copy turned onto its source by a whole number of
/// turns (the solver's own test, design note 0007 §11.1).
fn pattern_angle(v: f64, copies: &[(SketchId, SketchId, u32)]) -> Result<f64, String> {
    if !(v.is_finite() && v != 0.0 && v.abs() < TAU) {
        return Err("a pattern angle is more than 0° and less than 360°".into());
    }
    let last = copies.iter().map(|&(_, _, k)| k).max().unwrap_or(1);
    for k in 1..=last {
        let turns = v * f64::from(k) / TAU;
        if (turns - turns.round()).abs() * TAU <= 1e-9 {
            return Err(format!(
                "copy {k} would be turned onto its source ({}°)",
                trim_number((v * f64::from(k)).abs().to_degrees())
            ));
        }
    }
    Ok(v)
}

/// Text a prompt starts with for a value (millimetres or degrees, no unit).
fn prompt_text(c: &Constraint<SketchId>) -> String {
    let v = c.value().unwrap_or(0.0);
    match c {
        Constraint::Angle { .. } | Constraint::CircularPattern { .. } => {
            trim_number(v.abs().to_degrees())
        }
        Constraint::Rho { .. } => trim_number(v),
        _ => trim_number(v / LengthUnit::Millimetre.metres()),
    }
}

// ---------------------------------------------------------------------------
// The editor
// ---------------------------------------------------------------------------

/// An in-progress pointer drag.
#[derive(Clone, Debug)]
enum Drag {
    /// Moving points; `points` are where they started.
    Move {
        grab: DVec2,
        points: Vec<(SketchId, DVec2)>,
    },
    /// Changing a circle's radius.
    Radius { circle: SketchId },
    /// Changing a conic's ρ: its shoulder follows the pointer along the
    /// median from where it was grabbed, when ρ was `from`.
    Rho {
        conic: SketchId,
        grab: DVec2,
        from: f64,
    },
    /// Box selection.
    Box { from: Pos2, to: Pos2 },
    /// Moving a dimension's label: its offset from the dimension's anchor
    /// was `from` when the pointer was at `grab` (plane coordinates).
    Label {
        dimension: SketchId,
        grab: DVec2,
        from: DVec2,
        offset: DVec2,
    },
    /// A drawing tool's press-drag-release (a click at each end).
    Draw { from: Pos2 },
}

/// The state of sketch mode for one sketch.
#[derive(Clone, Debug)]
pub struct SketchEditor {
    /// The sketch being edited.
    pub feature: FeatureId,
    /// Current tool.
    pub tool: Tool,
    /// New geometry is construction geometry.
    pub construction: bool,
    /// Selected entities and constraints, in selection order.
    pub selection: Vec<SketchId>,
    /// What is under the pointer (Select tool).
    pub hover: Option<Pick>,
    /// Where the pointer snaps (drawing tools).
    pub cursor: Option<Snap>,
    /// The relation inferred for the line being drawn.
    pub inferred: Option<Inference>,
    /// Clicks of the drawing tool so far.
    pub clicks: Vec<Snap>,
    /// Number of sides for the polygon tool.
    pub polygon_sides: u32,
    /// Radius of the sketch fillet tool (metres).
    pub fillet_radius: f64,
    /// The last line the pointer passed over while drawing lines (for
    /// parallel and perpendicular inference).
    last_line: Option<SketchId>,
    /// Where the tangent arc being drawn starts.
    tangent_start: Option<TangentStart>,
    /// Signed angle the pointer has turned around a centre-point arc's
    /// centre since its start point.
    turn: f64,
    last_angle: Option<f64>,
    drag: Option<Drag>,
    preview: Option<Solution<SketchId>>,
    /// A dimension value being entered.
    pub prompt: Option<DimensionPrompt>,
    /// Feedback for the status bar.
    pub message: Option<String>,
    /// What the trim, extend or split tool would do at the pointer.
    pub curve_edit: Option<CurveEdit>,
    /// The last cuts found for trim, with the curves they were found in.
    cuts: Option<(Vec<(SketchId, ResolvedCurve)>, Cuts)>,
    /// The mirror tool's line, once chosen.
    pub mirror_axis: Option<SketchId>,
    /// What the mirror tool mirrors about the first line clicked: the
    /// selection when the tool was chosen.
    pub mirror_sources: Vec<SketchId>,
    /// The offset tool's chain, once chosen, with the curve clicked.
    pub offset_chain: Option<(SketchId, offset::Chain)>,
    /// What the pattern tools copy: the selection when the tool was
    /// chosen.
    pub pattern_sources: Vec<SketchId>,
    /// Instances a pattern tool makes, the sources included.
    pub pattern_count: u32,
    /// What the transform tools move: the selection when the tool was
    /// chosen.
    pub transform_sources: Vec<SketchId>,
    /// The transform tools keep the originals and move copies.
    pub transform_copy: bool,
    /// The transform tools keep constraints between the selection and the
    /// rest of the sketch (ADR 0014); otherwise they are removed.
    pub transform_keep: bool,
    /// The Move tool's typed offset (metres).
    pub transform_offset: DVec2,
    /// The Rotate tool's typed angle (radians, counter-clockwise).
    pub transform_angle: f64,
    /// The Scale tool's typed factor.
    pub transform_factor: f64,
    /// The text the Text tool places.
    pub text_content: String,
    /// The font the Text tool places it in.
    pub text_font: Font,
    /// The selected text's characters and font as being edited in the
    /// dialog (applied with [`SketchEditor::apply_text`]).
    pub text_edit: Option<TextEdit>,
}

/// A selected text's characters and font, as typed in the dialog.
#[derive(Clone, Debug, PartialEq)]
pub struct TextEdit {
    /// The text.
    pub entity: SketchId,
    /// Its characters.
    pub content: String,
    /// Its font.
    pub font: Font,
    /// Its characters and font when loaded from the document.
    pub loaded: (String, Font),
}

/// What the trim, extend or split tool would do at the pointer.
#[derive(Clone, Debug, PartialEq)]
pub struct CurveEdit {
    /// The curve under the pointer.
    pub target: SketchId,
    /// The edit, or why it cannot be made there.
    pub plan: Result<Plan, String>,
}

impl SketchEditor {
    /// Starts editing `feature` with the Select tool.
    pub fn new(feature: FeatureId) -> SketchEditor {
        SketchEditor {
            feature,
            tool: Tool::Select,
            construction: false,
            selection: Vec::new(),
            hover: None,
            cursor: None,
            inferred: None,
            clicks: Vec::new(),
            polygon_sides: 6,
            fillet_radius: 0.005,
            last_line: None,
            tangent_start: None,
            turn: 0.0,
            last_angle: None,
            drag: None,
            preview: None,
            prompt: None,
            message: None,
            curve_edit: None,
            cuts: None,
            mirror_axis: None,
            mirror_sources: Vec::new(),
            offset_chain: None,
            pattern_sources: Vec::new(),
            pattern_count: 3,
            transform_sources: Vec::new(),
            transform_copy: false,
            transform_keep: true,
            transform_offset: DVec2::ZERO,
            transform_angle: FRAC_PI_2,
            transform_factor: 2.0,
            text_content: text::DEFAULT_TEXT.to_string(),
            text_font: Font::DEFAULT,
            text_edit: None,
        }
    }

    /// The geometry to draw: the drag preview, or the regenerated sketch,
    /// with its placed dimension labels.
    pub fn display(&self, model: &Model) -> Option<SketchDisplay> {
        let r = model.sketch_result(self.feature)?;
        let mut d = match &self.preview {
            Some(p) => SketchDisplay::from_solution(self.feature, r.frame, p),
            None => SketchDisplay::from_result(r),
        };
        if let Some(def) = model.sketch_def(self.feature) {
            d.labels = def.labels.iter().copied().collect();
        }
        if let Some(Drag::Label {
            dimension, offset, ..
        }) = self.drag
        {
            d.labels.insert(dimension, offset);
        }
        Some(d)
    }

    /// The sketch plane.
    pub fn frame(&self, model: &Model) -> Option<Frame> {
        model.sketch_result(self.feature).map(|r| r.frame)
    }

    /// Switches tool, abandoning any clicks of the previous one.
    pub fn set_tool(&mut self, tool: Tool) {
        self.tool = tool;
        self.reset_tool();
        self.last_line = None;
        self.hover = None;
        if tool == Tool::Mirror {
            self.mirror_sources = self.selection.clone();
        }
        self.pattern_sources = if tool.is_pattern() {
            self.selection.clone()
        } else {
            Vec::new()
        };
        self.transform_sources = if tool.is_transform() {
            self.selection.clone()
        } else {
            Vec::new()
        };
    }

    fn reset_tool(&mut self) {
        self.clicks.clear();
        self.curve_edit = None;
        self.mirror_axis = None;
        self.mirror_sources.clear();
        self.offset_chain = None;
        self.inferred = None;
        self.tangent_start = None;
        self.turn = 0.0;
        self.last_angle = None;
    }

    /// Escape: abandons the current clicks, else returns to Select, else
    /// clears the selection.
    pub fn escape(&mut self) {
        if self.prompt.take().is_some() {
            return;
        }
        if self.drag.take().is_some() {
            self.preview = None;
            return;
        }
        if !self.clicks.is_empty() || self.mirror_axis.is_some() || self.offset_chain.is_some() {
            self.reset_tool();
        } else if self.tool != Tool::Select {
            self.set_tool(Tool::Select);
        } else {
            self.selection.clear();
        }
    }

    fn commit(&mut self, model: &mut Model, edits: Vec<SketchEdit>) -> Option<OpId> {
        match model.edit_sketch(self.feature, edits) {
            Ok(op) => {
                self.message = self.solve_message(model);
                Some(op)
            }
            Err(e) => {
                self.message = Some(format!("Edit refused: {e}"));
                None
            }
        }
    }

    /// Commits a shape tool's batch and resets the tool, or shows why the
    /// shape could not be made (keeping the clicks so far).
    fn commit_shape(&mut self, model: &mut Model, made: Result<Vec<SketchEdit>, String>) {
        match made {
            Ok(edits) => {
                if self.commit(model, edits).is_some() {
                    self.reset_tool();
                }
            }
            Err(e) => self.message = Some(e),
        }
    }

    /// A warning when the sketch no longer solves.
    fn solve_message(&self, model: &Model) -> Option<String> {
        match &model.sketch_result(self.feature)?.solution.status {
            SolveStatus::Failed(SolveError::Conflicting) => {
                Some("Constraints conflict (shown in red): delete one, or undo with Ctrl+Z".into())
            }
            SolveStatus::Failed(e) => Some(format!("Sketch does not solve: {e}")),
            _ => None,
        }
    }

    /// What to do next with the current tool.
    pub fn hint(&self) -> String {
        if self.tool == Tool::Mirror {
            return match (self.mirror_axis, self.mirror_sources.len()) {
                (Some(_), _) => {
                    "Click curves or points to mirror about the chosen line; Esc chooses another line"
                        .into()
                }
                (None, 0) => self.tool.hint(0).into(),
                (None, 1) => "Click the line to mirror the selected item about".into(),
                (None, n) => format!("Click the line to mirror the {n} selected items about"),
            };
        }
        if self.tool == Tool::Offset && self.offset_chain.is_some() {
            return "Move to the side to offset to and click; type the distance after; Esc chooses another curve"
                .into();
        }
        if self.tool.is_transform() {
            let items = match self.transform_sources.len() {
                0 => {
                    return format!(
                        "Select the curves and points to {}, then choose the tool",
                        self.tool.label().to_lowercase()
                    )
                }
                1 => "the selected item".to_string(),
                n => format!("the {n} selected items"),
            };
            return match (self.tool, self.clicks.len()) {
                (Tool::Move, 0) => {
                    format!("Click a point to move {items} from, or type the offsets and Apply")
                }
                (Tool::Move, _) => "Click where that point goes; Esc starts again".into(),
                (Tool::Rotate, 0) => format!("Click the centre to turn {items} about"),
                (Tool::Rotate, 1) => "Click a reference point, or type the angle and Apply".into(),
                (Tool::Rotate, _) => {
                    "Click where the reference point turns to; Esc starts again".into()
                }
                (_, 0) => format!("Click the base point to scale {items} about"),
                (_, 1) => "Click a reference point, or type the factor and Apply".into(),
                _ => "Click where the reference point goes; Esc starts again".into(),
            };
        }
        if self.tool.is_pattern() {
            let items = match self.pattern_sources.len() {
                0 => return "Select the curves and points to pattern, then choose the tool".into(),
                1 => "the selected item".to_string(),
                n => format!("the {n} selected items"),
            };
            return match (self.tool, self.clicks.len()) {
                (Tool::LinearPattern, 0) => {
                    format!("Click where the step between copies of {items} starts")
                }
                (Tool::LinearPattern, _) => {
                    "Click where the step to the first copy ends; Esc starts it again".into()
                }
                _ => format!("Click the centre to turn {items} about"),
            };
        }
        self.tool.hint(self.clicks.len()).into()
    }

    /// The sketch's state for the status bar: degrees of freedom, redundant
    /// or conflicting constraints.
    pub fn status(&self, model: &Model) -> String {
        let Some(r) = model.sketch_result(self.feature) else {
            return "The sketch is suppressed or was deleted".into();
        };
        let s = &r.solution;
        let mut text = match &s.status {
            SolveStatus::Solved if s.dof == 0 => "Fully defined".to_string(),
            SolveStatus::Solved => format!("{} degrees of freedom", s.dof),
            SolveStatus::Redundant => format!(
                "{} redundant constraint{} (orange) · {} degrees of freedom",
                s.redundant.len(),
                if s.redundant.len() == 1 { "" } else { "s" },
                s.dof
            ),
            SolveStatus::Failed(SolveError::Conflicting) => {
                "Conflicting constraints (red)".to_string()
            }
            SolveStatus::Failed(e) => format!("Not solved: {e}"),
        };
        if let Err(e) = &r.regions {
            if !matches!(s.status, SolveStatus::Failed(_)) {
                text.push_str(&format!(" · regions: {e}"));
            }
        }
        text
    }

    /// Updates hover, snapping and inference for the pointer at `pos`.
    pub fn hover(&mut self, model: &Model, map: &ScreenMap, pos: Option<Pos2>) {
        self.hover = None;
        self.cursor = None;
        self.inferred = None;
        self.curve_edit = None;
        let (Some(pos), Some(d)) = (pos, self.display(model)) else {
            return;
        };
        if self.tool == Tool::Select {
            if self.drag.is_none() {
                self.hover = pick(&d, &annotations(&d, map), map, pos);
            }
            return;
        }
        if self.tool == Tool::Mirror {
            self.curve_edit = self.plan_mirror(&d, map, pos);
            return;
        }
        if self.tool == Tool::Offset {
            self.curve_edit = self.plan_offset(model, &d, map, pos);
            return;
        }
        if self.tool.is_transform() {
            let Some(s) = snap(&d, map, pos, &[]) else {
                return;
            };
            self.cursor = Some(s);
            if let Some(t) = self.transform_at(s) {
                self.curve_edit = self.plan_transform(model, &d, t);
            }
            return;
        }
        if self.tool == Tool::SketchFillet {
            if let Some((corner, _)) = nearest_point(&d, map, pos, POINT_PICK, &[]) {
                let curves: Vec<(SketchId, ResolvedCurve)> =
                    d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
                let plan = fillet(&d.sketch, &curves, corner, self.fillet_radius);
                self.curve_edit = Some(CurveEdit {
                    target: corner,
                    plan,
                });
            }
            return;
        }
        if self.tool.edits_curves() {
            if let (Some((target, _)), Some(at)) =
                (nearest_curve(&d, map, pos, CURVE_PICK), map.to_plane(pos))
            {
                let plan = self.plan_curve_edit(model, &d, target, at);
                self.curve_edit = Some(CurveEdit { target, plan });
            }
            return;
        }
        let Some(mut s) = snap(&d, map, pos, &[]) else {
            return;
        };
        if self.tool == Tool::Line {
            if let SnapTarget::Curve(id) = s.target {
                if kind(&d.sketch, id) == Some(Kind::Line) {
                    self.last_line = Some(id);
                }
            }
            if s.target == SnapTarget::Free {
                if let Some(start) = self.clicks.first().copied() {
                    if let Some((at, inference)) = self.infer(&d, &start, s.at) {
                        s.at = at;
                        self.inferred = Some(inference);
                    }
                }
            }
        }
        if self.tool == Tool::CenterArc && self.clicks.len() == 2 {
            let a = (s.at - self.clicks[0].at).to_angle();
            if let Some(last) = self.last_angle {
                // Wrap the step into (−π, π].
                let mut step = (a - last).rem_euclid(TAU);
                if step > PI {
                    step -= TAU;
                }
                self.turn += step;
            } else {
                let a0 = (self.clicks[1].at - self.clicks[0].at).to_angle();
                let mut step = (a - a0).rem_euclid(TAU);
                if step > PI {
                    step -= TAU;
                }
                self.turn = step;
            }
            self.last_angle = Some(a);
        }
        if self.tool == Tool::EllipticalArc && self.clicks.len() == 4 {
            if let Some((a, b)) =
                shapes::ellipse_axes(self.clicks[0].at, self.clicks[1].at, self.clicks[2].at)
            {
                let c = self.clicks[0].at;
                let t = shapes::eccentric_anomaly(c, a, b, s.at);
                let last = self
                    .last_angle
                    .unwrap_or_else(|| shapes::eccentric_anomaly(c, a, b, self.clicks[3].at));
                // Wrap the step into (−π, π].
                let mut step = (t - last).rem_euclid(TAU);
                if step > PI {
                    step -= TAU;
                }
                self.turn += step;
                self.last_angle = Some(t);
            }
        }
        if self.tool.is_pattern() {
            self.curve_edit = self.plan_pattern(model, &d, s);
        }
        // A control-point spline closes on its first point.
        if matches!(self.tool, Tool::Spline | Tool::ControlSpline) && self.clicks.len() >= 3 {
            let first = self.clicks[0];
            if map
                .to_screen(first.at)
                .is_some_and(|f| f.distance(pos) <= POINT_PICK)
            {
                s = first;
            }
        }
        self.cursor = Some(s);
    }

    /// The mirror tool's plan at `pos`: while choosing the line, the line
    /// under the pointer (with the selection's images about it, if any);
    /// once it is chosen, the image of the curve or point under the
    /// pointer.
    fn plan_mirror(&self, d: &SketchDisplay, map: &ScreenMap, pos: Pos2) -> Option<CurveEdit> {
        let curves: Vec<(SketchId, ResolvedCurve)> =
            d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
        match self.mirror_axis {
            None => {
                let (target, _) = nearest_curve(d, map, pos, CURVE_PICK)?;
                let sources = mirror::sources(&d.sketch, &self.mirror_sources);
                let plan = if kind(&d.sketch, target) != Some(Kind::Line) {
                    Err("The mirror line must be a line".to_string())
                } else if sources.is_empty() {
                    // Choosing the line: nothing to show but the line.
                    Ok(Plan {
                        edits: Vec::new(),
                        removed: 0,
                        preview: trim::Preview::Shapes {
                            curves: Vec::new(),
                            points: Vec::new(),
                        },
                    })
                } else {
                    mirror::mirror(&d.sketch, &curves, &sources, target, true).map(|m| m.plan)
                };
                Some(CurveEdit { target, plan })
            }
            Some(axis) => {
                // A lone point under the pointer, else a curve.
                let lone = nearest_point(d, map, pos, POINT_PICK, &[])
                    .map(|(p, _)| p)
                    .filter(|p| {
                        !d.sketch
                            .entities
                            .iter()
                            .any(|e| e.geometry.referenced_points().contains(p))
                    });
                let target =
                    lone.or_else(|| nearest_curve(d, map, pos, CURVE_PICK).map(|(c, _)| c))?;
                let plan = if target == axis {
                    Err("This is the mirror line".to_string())
                } else {
                    mirror::mirror(&d.sketch, &curves, &[target], axis, true).map(|m| m.plan)
                };
                Some(CurveEdit { target, plan })
            }
        }
    }

    /// Commits the mirror tool's click at `pos`: chooses the line, or
    /// mirrors the selection about it, or the item under the pointer about
    /// the chosen line.
    fn click_mirror(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        self.hover(model, map, Some(pos));
        let Some(d) = self.display(model) else {
            return;
        };
        let Some(edit) = self.curve_edit.take() else {
            self.message = Some(match self.mirror_axis {
                None => "Click a line to mirror about".into(),
                Some(_) => "Click a curve or point to mirror".into(),
            });
            return;
        };
        if let Err(e) = edit.plan {
            self.message = Some(e);
            return;
        }
        let (axis, sources) = match self.mirror_axis {
            Some(axis) => (axis, vec![edit.target]),
            None => {
                let sources = mirror::sources(&d.sketch, &self.mirror_sources);
                if sources.is_empty() {
                    self.mirror_axis = Some(edit.target);
                    self.message = Some("Mirror line chosen".into());
                    self.hover(model, map, Some(pos));
                    return;
                }
                (edit.target, sources)
            }
        };
        let curves: Vec<(SketchId, ResolvedCurve)> =
            d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
        // Images share the source's points on the line where they can;
        // where that would repeat what a constraint already says, they get
        // points of their own (which the mirror keeps on the same spot).
        let mut outcome = Err(String::new());
        for share in [true, false] {
            let planned = mirror::mirror(&d.sketch, &curves, &sources, axis, share);
            outcome = planned.and_then(|m| {
                model
                    .edit_sketch_checked(self.feature, m.plan.edits.clone())
                    .map(|_| m.message())
            });
            if outcome.is_ok() {
                break;
            }
        }
        match outcome {
            Ok(done) => {
                self.message = self.solve_message(model).or(Some(done));
                self.mirror_axis = Some(axis);
                self.mirror_sources.clear();
                self.selection.clear();
            }
            Err(e) => self.message = Some(format!("Not done: {e}")),
        }
        self.hover(model, map, Some(pos));
    }

    /// The offset of the chosen chain for the pointer at `pos`: the side
    /// it is on and its distance from the clicked curve, rounded to a step
    /// of about 1.5 screen points.
    fn offset_at(
        &self,
        d: &SketchDisplay,
        map: &ScreenMap,
        pos: Pos2,
    ) -> Option<Result<offset::Offsetted, String>> {
        let (clicked, chain) = self.offset_chain.as_ref()?;
        let curves: Vec<(SketchId, ResolvedCurve)> =
            d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
        let at = map.to_plane(pos)?;
        let curve = &d.curve(*clicked)?.curve;
        let (distance, left) = offset::side_and_distance(curve, at)?;
        let distance = round_distance(distance, map.plane_per_point(pos));
        Some(offset::offset(&d.sketch, &curves, chain, distance, left))
    }

    /// The offset tool's plan at `pos`: before a curve is chosen, the chain
    /// of the curve under the pointer; after, the chain's offset towards
    /// the pointer, checked by solving it.
    fn plan_offset(
        &self,
        model: &Model,
        d: &SketchDisplay,
        map: &ScreenMap,
        pos: Pos2,
    ) -> Option<CurveEdit> {
        match &self.offset_chain {
            None => {
                let (target, _) = nearest_curve(d, map, pos, CURVE_PICK)?;
                let plan = offset::chain(&d.sketch, target).map(|chain| Plan {
                    edits: Vec::new(),
                    removed: 0,
                    preview: trim::Preview::Shapes {
                        curves: chain
                            .pieces
                            .iter()
                            .filter_map(|(id, _)| d.curve(*id).map(|c| c.curve.clone()))
                            .collect(),
                        points: Vec::new(),
                    },
                });
                Some(CurveEdit { target, plan })
            }
            Some((clicked, _)) => {
                let plan = self.offset_at(d, map, pos)?.and_then(|o| {
                    model
                        .check_sketch_edit(self.feature, &o.plan.edits)
                        .map(|_| o.plan)
                });
                Some(CurveEdit {
                    target: *clicked,
                    plan,
                })
            }
        }
    }

    /// Commits the offset tool's click at `pos`: chooses the chain, or
    /// offsets it towards the pointer and opens its distance for typing.
    fn click_offset(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        self.hover(model, map, Some(pos));
        let Some(d) = self.display(model) else {
            return;
        };
        let Some(edit) = self.curve_edit.take() else {
            self.message = Some(match self.offset_chain {
                None => "Click a line, arc or circle to offset".into(),
                Some(_) => "Move the pointer to the side to offset to".into(),
            });
            return;
        };
        if self.offset_chain.is_none() {
            match edit
                .plan
                .and_then(|_| offset::chain(&d.sketch, edit.target))
            {
                Ok(chain) => {
                    self.message = Some(match chain.pieces.len() {
                        1 => "Curve chosen".into(),
                        n => format!("Chain of {n} curves chosen"),
                    });
                    self.offset_chain = Some((edit.target, chain));
                    self.hover(model, map, Some(pos));
                }
                Err(e) => self.message = Some(e),
            }
            return;
        }
        let outcome = match self.offset_at(&d, map, pos) {
            Some(planned) => planned.and_then(|o| {
                model
                    .edit_sketch_checked(self.feature, o.plan.edits.clone())
                    .map(|op| (o, op))
            }),
            None => Err("Move the pointer to the side to offset to".into()),
        };
        match outcome {
            Ok((o, op)) => {
                self.message = self.solve_message(model).or(Some(o.message()));
                self.offset_chain = None;
                self.selection.clear();
                let id = SketchId {
                    op,
                    index: o.constraint.index,
                };
                self.edit_dimension(model, id, pos);
            }
            Err(e) => self.message = Some(format!("Not done: {e}")),
        }
        self.hover(model, map, Some(pos));
    }

    /// The similarity the transform tool's clicks and the pointer snapped
    /// at `s` give: Move from the first click to `s`; Rotate about the
    /// first click from the second's direction to `s`'s; Scale about the
    /// first click by `s`'s distance over the second's.
    fn transform_at(&self, s: Snap) -> Option<Similarity> {
        let c = &self.clicks;
        match (self.tool, c.len()) {
            (Tool::Move, 1) => {
                let v = s.at - c[0].at;
                (v.length() > SAME_POINT).then(|| Similarity::translation(v))
            }
            (Tool::Rotate, 2) => {
                let (a, b) = (c[1].at - c[0].at, s.at - c[0].at);
                (a.length() > SAME_POINT && b.length() > SAME_POINT)
                    .then(|| Similarity::rotation(c[0].at, a.perp_dot(b).atan2(a.dot(b))))
            }
            (Tool::Scale, 2) => {
                let (a, b) = ((c[1].at - c[0].at).length(), (s.at - c[0].at).length());
                (a > SAME_POINT && b > SAME_POINT).then(|| Similarity::scaling(c[0].at, b / a))
            }
            _ => None,
        }
    }

    /// The transform tool's plan for `t` (not checked).
    fn transform_plan(
        &self,
        model: &Model,
        d: &SketchDisplay,
        t: Similarity,
    ) -> Result<transform::Transformed, String> {
        let curves: Vec<(SketchId, ResolvedCurve)> =
            d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
        let labels = model
            .sketch_def(self.feature)
            .map(|def| def.labels.clone())
            .unwrap_or_default();
        let sources = mirror::sources(&d.sketch, &self.transform_sources);
        transform::transform(
            &d.sketch,
            &labels,
            &curves,
            &sources,
            t,
            transform::Options {
                copy: self.transform_copy,
                keep: self.transform_keep,
            },
        )
    }

    /// The transform tools' preview for `t`, checked like the edit, with
    /// the first source as the target.
    fn plan_transform(&self, model: &Model, d: &SketchDisplay, t: Similarity) -> Option<CurveEdit> {
        let target = *self.transform_sources.first()?;
        let plan = self.transform_plan(model, d, t).and_then(|p| {
            model
                .check_sketch_edit(self.feature, &p.plan.edits)
                .map(|_| p.plan)
        });
        Some(CurveEdit { target, plan })
    }

    /// Commits a transform tool's click at `pos`: the clicks before the
    /// last set the points; the last transforms the selection.
    fn click_transform(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        self.hover(model, map, Some(pos));
        let Some(s) = self.cursor else {
            return;
        };
        if self.transform_sources.is_empty() {
            self.message = Some(format!(
                "Select the curves and points to {} first",
                self.tool.label().to_lowercase()
            ));
            return;
        }
        let needed = if self.tool == Tool::Move { 1 } else { 2 };
        if self.clicks.len() < needed {
            if self
                .clicks
                .first()
                .is_some_and(|c| (c.at - s.at).length() <= SAME_POINT)
            {
                self.message = Some("Click a point away from the centre".into());
                return;
            }
            self.clicks.push(s);
            self.hover(model, map, Some(pos));
            return;
        }
        match self.transform_at(s) {
            Some(t) => self.commit_transform(model, t),
            None => self.message = Some("Click a point away from the first one".into()),
        }
    }

    /// Keeps [`SketchEditor::text_edit`] on the selected text: loaded when
    /// a single text is selected (with the Select tool) or when the
    /// document's text changed under it (undo, say), cleared otherwise.
    pub fn sync_text_edit(&mut self, model: &Model) {
        let selected = match (self.tool, self.selection.as_slice()) {
            (Tool::Select, [id]) => model
                .sketch_def(self.feature)
                .and_then(|def| Some((*id, text::text_of(&def.sketch, *id)?))),
            _ => None,
        };
        let Some((entity, (content, font))) = selected else {
            self.text_edit = None;
            return;
        };
        let font = Font::from_name(&font).unwrap_or(Font::DEFAULT);
        let current = (content.clone(), font);
        if self
            .text_edit
            .as_ref()
            .is_none_or(|e| e.entity != entity || e.loaded != current)
        {
            self.text_edit = Some(TextEdit {
                entity,
                content,
                font,
                loaded: current,
            });
        }
    }

    /// Gives the selected text (and its mirror images and pattern copies)
    /// the characters and font typed in the dialog, in one undo step; the
    /// box keeps its corner, height and direction and takes the new width.
    pub fn apply_text(&mut self, model: &mut Model) {
        let Some(edit) = self.text_edit.clone() else {
            return;
        };
        let Some(def) = model.sketch_def(self.feature) else {
            return;
        };
        let done = text::set_text(&def.sketch, edit.entity, &edit.content, edit.font)
            .and_then(|edits| model.edit_sketch_checked(self.feature, edits));
        self.message = match done {
            Ok(_) => {
                // The document has the typed text now; an undo reloads.
                self.text_edit = Some(TextEdit {
                    loaded: (edit.content.clone(), edit.font),
                    ..edit
                });
                self.solve_message(model).or(Some("Text changed".into()))
            }
            Err(e) => Some(format!("Text not changed: {e}")),
        };
    }

    /// Applies the transform tool's typed amount: the offsets (Move), the
    /// angle about the clicked centre (Rotate) or the factor about the
    /// clicked base point (Scale), in one undo step.
    pub fn apply_typed_transform(&mut self, model: &mut Model) {
        if self.transform_sources.is_empty() {
            self.message = Some(format!(
                "Select the curves and points to {} first",
                self.tool.label().to_lowercase()
            ));
            return;
        }
        let t = match (self.tool, self.clicks.first()) {
            (Tool::Move, _) => Similarity::translation(self.transform_offset),
            (Tool::Rotate, Some(c)) => Similarity::rotation(c.at, self.transform_angle),
            (Tool::Scale, Some(c)) => Similarity::scaling(c.at, self.transform_factor),
            (Tool::Rotate, None) => {
                self.message = Some("Click the centre to turn about first".into());
                return;
            }
            (Tool::Scale, None) => {
                self.message = Some("Click the base point to scale about first".into());
                return;
            }
            _ => return,
        };
        self.commit_transform(model, t);
    }

    /// Transforms the selection by `t` through the checked edit, then
    /// returns to Select with the moved entities (or the copies) selected.
    fn commit_transform(&mut self, model: &mut Model, t: Similarity) {
        let Some(d) = self.display(model) else {
            return;
        };
        let outcome = self.transform_plan(model, &d, t).and_then(|p| {
            model
                .edit_sketch_checked(self.feature, p.plan.edits.clone())
                .map(|op| (p, op))
        });
        match outcome {
            Ok((p, op)) => {
                let verb = match (self.transform_copy, self.tool) {
                    (true, _) => "Copied",
                    (_, Tool::Rotate) => "Rotated",
                    (_, Tool::Scale) => "Scaled",
                    _ => "Moved",
                };
                self.message = self.solve_message(model).or(Some(p.message(verb)));
                let result: Vec<SketchId> = p
                    .result
                    .iter()
                    .map(|id| {
                        if id.is_new_in_batch() {
                            SketchId {
                                op,
                                index: id.index,
                            }
                        } else {
                            *id
                        }
                    })
                    .collect();
                self.set_tool(Tool::Select);
                self.selection = result;
            }
            Err(e) => self.message = Some(format!("Not done: {e}")),
        }
    }

    /// The pattern for the pointer snapped at `s`, once the anchor is
    /// known: a linear pattern's start clicked (and the pointer elsewhere),
    /// or a circular pattern's centre under the pointer.
    fn pattern_at(&self, d: &SketchDisplay, s: Snap) -> Option<Result<pattern::Patterned, String>> {
        let anchor = match (self.tool, self.clicks.first()) {
            (Tool::LinearPattern, Some(from)) if (s.at - from.at).length() > SAME_POINT => {
                pattern::Anchor::Linear { from: *from, to: s }
            }
            (Tool::CircularPattern, _) => pattern::Anchor::Circular { centre: s },
            _ => return None,
        };
        let curves: Vec<(SketchId, ResolvedCurve)> =
            d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
        let sources = mirror::sources(&d.sketch, &self.pattern_sources);
        Some(pattern::pattern(
            &d.sketch,
            &curves,
            &sources,
            anchor,
            self.pattern_count,
        ))
    }

    /// The pattern tools' plan for the pointer snapped at `s`: the copies,
    /// checked by solving them, with the first source as the target.
    fn plan_pattern(&self, model: &Model, d: &SketchDisplay, s: Snap) -> Option<CurveEdit> {
        let target = *self.pattern_sources.first()?;
        let plan = self.pattern_at(d, s)?.and_then(|p| {
            model
                .check_sketch_edit(self.feature, &p.plan.edits)
                .map(|_| p.plan)
        });
        Some(CurveEdit { target, plan })
    }

    /// Commits a pattern tool's click at `pos`: a linear pattern's first
    /// click starts the step; the next (or a circular pattern's centre)
    /// makes the copies, returns to the Select tool and opens the step's
    /// length or the angle for typing.
    fn click_pattern(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        self.hover(model, map, Some(pos));
        let Some(s) = self.cursor else {
            return;
        };
        if self.pattern_sources.is_empty() {
            self.message = Some("Select the curves and points to pattern first".into());
            return;
        }
        if self.tool == Tool::LinearPattern && self.clicks.is_empty() {
            self.clicks.push(s);
            self.hover(model, map, Some(pos));
            return;
        }
        let Some(d) = self.display(model) else {
            return;
        };
        let Some(planned) = self.pattern_at(&d, s) else {
            return;
        };
        let outcome = planned.and_then(|p| {
            model
                .edit_sketch_checked(self.feature, p.plan.edits.clone())
                .map(|op| (p, op))
        });
        match outcome {
            Ok((p, op)) => {
                self.message = self.solve_message(model).or(Some(p.message()));
                self.selection.clear();
                self.set_tool(Tool::Select);
                match p.line {
                    Some(line) => self.prompt_length(
                        model,
                        SketchId {
                            op,
                            index: line.index,
                        },
                        pos,
                    ),
                    None => self.edit_dimension(
                        model,
                        SketchId {
                            op,
                            index: p.constraint.index,
                        },
                        pos,
                    ),
                }
            }
            Err(e) => self.message = Some(format!("Not done: {e}")),
        }
        self.hover(model, map, Some(pos));
    }

    /// Opens the prompt for a new length dimension on `line`, at its
    /// measured length.
    fn prompt_length(&mut self, model: &Model, line: SketchId, at: Pos2) {
        let Some(d) = self.display(model) else {
            return;
        };
        let mut constraint = Constraint::Length { line, value: 0.0 };
        if let Ok(v) = measure(&d.sketch, &constraint) {
            constraint.set_value(v);
            self.prompt = Some(DimensionPrompt {
                text: prompt_text(&constraint),
                constraint,
                existing: None,
                at,
                error: None,
            });
        }
    }

    /// Plans the trim, extend or split of `target` at plane point `at`.
    fn plan_curve_edit(
        &mut self,
        model: &Model,
        d: &SketchDisplay,
        target: SketchId,
        at: DVec2,
    ) -> Result<Plan, String> {
        if text::is_text(&d.sketch, target) {
            return Err(text::NOT_EDITABLE.into());
        }
        let curves: Vec<(SketchId, ResolvedCurve)> =
            d.curves.iter().map(|c| (c.id, c.curve.clone())).collect();
        // Letters cut curves, but no end can be held on them.
        let texts = text::texts(&d.sketch);
        match self.tool {
            Tool::Trim => {
                // Finding cuts intersects the curve with every other one;
                // reuse them while the pointer stays on the same curve.
                let fresh = !matches!(&self.cuts,
                    Some((c, k)) if k.target == target && *c == curves);
                if fresh {
                    let k = trim::cuts(model.kernel(), &curves, target, &texts)?;
                    self.cuts = Some((curves, k));
                }
                let (_, k) = self.cuts.as_ref().expect("cuts were just found");
                trim::trim(&d.sketch, k, at)
            }
            Tool::Extend => trim::extend(model.kernel(), &d.sketch, &curves, target, at, &texts),
            Tool::Split => trim::split(&d.sketch, &curves, target, at),
            _ => Err("not a curve editing tool".into()),
        }
    }

    /// Commits the trim, extend or split planned at `pos`.
    fn edit_curve(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        self.hover(model, map, Some(pos));
        match self.curve_edit.take() {
            Some(CurveEdit { plan: Ok(plan), .. }) => {
                let verb = match self.tool {
                    Tool::Trim => "Trimmed",
                    Tool::Extend => "Extended",
                    Tool::SketchFillet => "Filleted",
                    _ => "Split",
                };
                match model.edit_sketch_checked(self.feature, plan.edits) {
                    Ok(_) => {
                        self.message = self.solve_message(model).or(Some(match plan.removed {
                            0 => verb.to_string(),
                            1 => format!("{verb}; removed 1 constraint that no longer applied"),
                            n => format!("{verb}; removed {n} constraints that no longer applied"),
                        }));
                        if let Some(def) = model.sketch_def(self.feature) {
                            let sketch = &def.sketch;
                            self.selection.retain(|id| {
                                sketch.entity(*id).is_some() || sketch.constraint(*id).is_some()
                            });
                        }
                    }
                    Err(e) => self.message = Some(format!("Not done: {e}")),
                }
            }
            Some(CurveEdit { plan: Err(e), .. }) => self.message = Some(e),
            None if self.tool == Tool::SketchFillet => {
                self.message = Some("Click the corner point where two lines or arcs meet".into())
            }
            None => self.message = Some("Click a curve".into()),
        }
        self.hover(model, map, Some(pos));
    }

    /// The relation a line from `start` to `to` is close to (design note
    /// 0007 §5), with `to` moved onto it.
    fn infer(&self, d: &SketchDisplay, start: &Snap, to: DVec2) -> Option<(DVec2, Inference)> {
        let v = to - start.at;
        if v.length() <= SAME_POINT {
            return None;
        }
        let mut candidates = vec![
            (DVec2::X, Inference::Horizontal),
            (DVec2::Y, Inference::Vertical),
        ];
        let line_dir = |id: SketchId| match d.curve(id).map(|c| &c.curve) {
            Some(ResolvedCurve::Line { start, end }) => {
                Some((*end - *start).normalize_or_zero()).filter(|w| *w != DVec2::ZERO)
            }
            _ => None,
        };
        if let SnapTarget::Point(p) = start.target {
            if let Some((curve, at_end)) = curve_ending_at(&d.sketch, p) {
                match kind(&d.sketch, curve) {
                    Some(Kind::Line) => {
                        if let Some(w) = line_dir(curve) {
                            candidates.push((perp(w), Inference::Perpendicular(curve)));
                            candidates.push((w, Inference::Parallel(curve)));
                        }
                    }
                    Some(Kind::Arc | Kind::Spline | Kind::Ellipse | Kind::Conic) => {
                        if let Some(w) = d
                            .curve(curve)
                            .and_then(|c| leaving_direction(&c.curve, at_end))
                        {
                            candidates.push((w, Inference::Tangent(curve)));
                        }
                    }
                    _ => {}
                }
            }
        }
        if let Some(l) = self.last_line {
            if let Some(w) = line_dir(l) {
                candidates.push((perp(w), Inference::Perpendicular(l)));
                candidates.push((w, Inference::Parallel(l)));
            }
        }
        candidates.into_iter().find_map(|(w, inference)| {
            let along = v.dot(w);
            (v.perp_dot(w).abs() <= INFER_ANGLE.tan() * along.abs() && along != 0.0)
                .then(|| (start.at + w * along, inference))
        })
    }

    /// Geometry of the tool in progress, as plane polylines.
    pub fn rubber_band(&self) -> Vec<Vec<DVec2>> {
        let Some(cur) = self.cursor else {
            return Vec::new();
        };
        let c = &self.clicks;
        match (self.tool, c.len()) {
            (Tool::Line, 1) => vec![vec![c[0].at, cur.at]],
            (Tool::Rectangle, 1) => {
                let (a, b) = (c[0].at, cur.at);
                vec![vec![a, DVec2::new(b.x, a.y), b, DVec2::new(a.x, b.y), a]]
            }
            (Tool::Circle, 1) => {
                vec![arc_points(c[0].at, (cur.at - c[0].at).length(), 0.0, TAU)]
            }
            (Tool::Arc, 1) => vec![vec![c[0].at, cur.at]],
            (Tool::Arc, 2) => match circumcircle(c[0].at, c[1].at, cur.at) {
                Some((center, r)) => {
                    let (a0, a1, am) = (
                        (c[0].at - center).to_angle(),
                        (c[1].at - center).to_angle(),
                        (cur.at - center).to_angle(),
                    );
                    let sweep = ccw(a0, a1);
                    let sweep = if ccw(a0, am) < sweep {
                        sweep
                    } else {
                        sweep - TAU
                    };
                    vec![arc_points(center, r, a0, sweep)]
                }
                None => vec![vec![c[0].at, c[1].at]],
            },
            (Tool::CenterArc, 1)
            | (Tool::LinearPattern, 1)
            | (Tool::Move, 1)
            | (Tool::Rotate, 1)
            | (Tool::Scale, 1) => vec![vec![c[0].at, cur.at]],
            (Tool::Rotate, 2) | (Tool::Scale, 2) => vec![vec![c[1].at, c[0].at, cur.at]],
            (Tool::CenterArc, 2) => {
                let r = (c[1].at - c[0].at).length();
                let a0 = (c[1].at - c[0].at).to_angle();
                vec![
                    arc_points(c[0].at, r, a0, self.turn.clamp(-TAU, TAU)),
                    vec![c[0].at, cur.at],
                ]
            }
            (Tool::Ellipse | Tool::EllipticalArc, 1) => vec![vec![c[0].at, cur.at]],
            (Tool::Ellipse | Tool::EllipticalArc, n @ 2..=4) => {
                let (centre, m) = (c[0].at, c[1].at);
                let extent = if n == 2 { cur.at } else { c[2].at };
                let Some((a, b)) = shapes::ellipse_axes(centre, m, extent) else {
                    return vec![vec![centre, m]];
                };
                let axes = vec![centre + b, centre, m];
                let along = |p: DVec2| {
                    shapes::ellipse_point(centre, a, b, shapes::eccentric_anomaly(centre, a, b, p))
                };
                match n {
                    2 | 3 => {
                        let mut out = vec![shapes::ellipse_polyline(centre, a, b, 0.0, TAU), axes];
                        if n == 3 {
                            // Where the arc will start.
                            out.push(vec![centre, along(cur.at)]);
                        }
                        out
                    }
                    _ => {
                        let t0 = shapes::eccentric_anomaly(centre, a, b, c[3].at);
                        vec![
                            shapes::ellipse_polyline(centre, a, b, t0, self.turn.clamp(-TAU, TAU)),
                            axes,
                            vec![along(c[3].at), centre, along(cur.at)],
                        ]
                    }
                }
            }
            (Tool::CenterRectangle, 1) => {
                let k = shapes::centre_rectangle_corners(c[0].at, cur.at);
                vec![vec![k[0], k[1], k[2], k[3], k[0]], vec![k[0], k[2]]]
            }
            (Tool::AlignedRectangle | Tool::Slot, 1) => vec![vec![c[0].at, cur.at]],
            (Tool::AlignedRectangle, 2) => {
                match shapes::aligned_rectangle_corners(c[0].at, c[1].at, cur.at) {
                    Some(k) => vec![vec![k[0], k[1], k[2], k[3], k[0]]],
                    None => vec![vec![c[0].at, c[1].at]],
                }
            }
            (Tool::Slot, 2) => match shapes::slot_geometry(c[0].at, c[1].at, cur.at) {
                Some(g) => vec![shapes::slot_outline(&g), vec![c[0].at, c[1].at]],
                None => vec![vec![c[0].at, c[1].at]],
            },
            (Tool::Polygon, 1) => {
                let mut v = shapes::polygon_vertices(c[0].at, cur.at, self.polygon_sides);
                if let Some(first) = v.first().copied() {
                    v.push(first);
                }
                vec![v]
            }
            (Tool::TangentArc, 1) => match &self.tangent_start {
                Some(t) => vec![shapes::tangent_arc_outline(t, cur.at)],
                None => Vec::new(),
            },
            (Tool::Spline | Tool::ControlSpline, n) if n >= 3 && !distinct_clicks(&c[0], &cur) => {
                // Back on the first point: the smooth closed curve it
                // would make (ADR 0016), with the control polygon.
                let mut pts: Vec<DVec2> = c.iter().map(|s| s.at).collect();
                pts.push(c[0].at);
                let control = self.tool == Tool::ControlSpline;
                let mut out = Vec::new();
                if control {
                    out.push(pts.clone());
                }
                match closed_spline(&pts, control) {
                    Some(curve) => out.push(sample_resolved(&curve)),
                    None if !control => out.push(pts),
                    None => {}
                }
                out
            }
            (Tool::Spline, n) if n >= 1 => {
                let mut pts: Vec<DVec2> = c.iter().map(|s| s.at).collect();
                pts.push(cur.at);
                pts.dedup_by(|a, b| (*a - *b).length() <= SAME_POINT);
                match ProfileCurve::interpolate(&pts) {
                    Ok(curve) if pts.len() >= 2 => vec![sample_curve(&curve)],
                    _ => vec![pts],
                }
            }
            (Tool::Text, 1) => match text::layout(&self.text_content, self.text_font) {
                Ok(layout) => text::preview(c[0].at, cur.at, &layout),
                Err(_) => vec![vec![c[0].at, cur.at]],
            },
            (Tool::Conic, 1) => vec![vec![c[0].at, cur.at]],
            (Tool::Conic, 2) => {
                // The control polygon to the pointer and the parabola
                // (ρ = ½) it would make.
                let (a, e) = (c[0].at, c[1].at);
                if shapes::conic_control_ok(a, e, cur.at) {
                    vec![
                        shapes::conic_polyline([a, cur.at, e], 0.5),
                        vec![a, cur.at, e],
                    ]
                } else {
                    vec![vec![a, e]]
                }
            }
            (Tool::Conic, 3) => {
                // The conic whose shoulder follows the pointer along the
                // median, with the control polygon and the median.
                let (a, k, e) = (c[0].at, c[2].at, c[1].at);
                let rho = shapes::conic_rho(a, k, e, cur.at);
                vec![
                    shapes::conic_polyline([a, k, e], rho),
                    vec![a, k, e],
                    vec![(a + e) * 0.5, k],
                ]
            }
            (Tool::ControlSpline, n) if n >= 1 => {
                // The control polygon and, from three points, the curve.
                let mut pts: Vec<DVec2> = c.iter().map(|s| s.at).collect();
                if distinct_clicks(&c[n - 1], &cur) {
                    pts.push(cur.at);
                }
                let mut out = vec![pts.clone()];
                if let Some(curve) = control_spline(&pts) {
                    out.push(sample_curve(&curve));
                }
                out
            }
            _ => Vec::new(),
        }
    }

    /// A click at `pos`. `add` (Shift) extends the selection.
    pub fn click(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2, add: bool) {
        let Some(d) = self.display(model) else {
            return;
        };
        if self.tool == Tool::Select {
            let picked = pick(&d, &annotations(&d, map), map, pos);
            match picked {
                Some(p) => {
                    let id = p.id();
                    if add {
                        if let Some(i) = self.selection.iter().position(|s| *s == id) {
                            self.selection.remove(i);
                        } else {
                            self.selection.push(id);
                        }
                    } else {
                        self.selection = vec![id];
                    }
                }
                None if !add => self.selection.clear(),
                None => {}
            }
            return;
        }
        if self.tool == Tool::Mirror {
            self.click_mirror(model, map, pos);
            return;
        }
        if self.tool == Tool::Offset {
            self.click_offset(model, map, pos);
            return;
        }
        if self.tool.is_pattern() {
            self.click_pattern(model, map, pos);
            return;
        }
        if self.tool.is_transform() {
            self.click_transform(model, map, pos);
            return;
        }
        if self.tool.edits_curves() {
            self.edit_curve(model, map, pos);
            return;
        }
        // Use the hover state when it is for this position, so a click
        // lands where the user saw the snap.
        self.hover(model, map, Some(pos));
        if let Some(s) = self.cursor {
            self.tool_click(model, s);
        }
    }

    fn tool_click(&mut self, model: &mut Model, s: Snap) {
        let distinct = distinct_clicks;
        let construction = self.construction;
        match self.tool {
            // Curve editing tools act on clicks directly (`edit_curve`).
            Tool::Select
            | Tool::Trim
            | Tool::Extend
            | Tool::Split
            | Tool::SketchFillet
            | Tool::Mirror
            | Tool::Offset
            | Tool::LinearPattern
            | Tool::CircularPattern
            | Tool::Move
            | Tool::Rotate
            | Tool::Scale => {}
            Tool::Text => {
                let Some(first) = self.clicks.first().copied() else {
                    // Say at once if the text cannot be laid out.
                    match text::layout(&self.text_content, self.text_font) {
                        Ok(_) => self.clicks.push(s),
                        Err(e) => self.message = Some(format!("Text not placed: {e}")),
                    }
                    return;
                };
                let mut b = Batch::default();
                let made = text::place(
                    &mut b,
                    &first,
                    s.at,
                    &self.text_content,
                    self.text_font,
                    construction,
                );
                self.commit_shape(model, made.map(|()| b.edits));
            }
            Tool::Point => {
                let mut b = Batch::default();
                let p = b.snapped_point(&s);
                if construction {
                    b.edits.push(SketchEdit::SetConstruction {
                        entity: p,
                        construction: true,
                    });
                }
                self.commit(model, b.edits);
            }
            Tool::Line => {
                let Some(start) = self.clicks.first().copied() else {
                    self.clicks.push(s);
                    return;
                };
                if !distinct(&start, &s) {
                    return;
                }
                let mut b = Batch::default();
                let p0 = b.snapped_point(&start);
                let p1 = b.snapped_point(&s);
                let line = b.entity(Geometry::Line { start: p0, end: p1 }, construction);
                match self.inferred {
                    Some(Inference::Horizontal) => {
                        b.constraint(Constraint::Horizontal { line });
                    }
                    Some(Inference::Vertical) => {
                        b.constraint(Constraint::Vertical { line });
                    }
                    Some(Inference::Parallel(m)) => {
                        b.constraint(Constraint::Parallel { a: m, b: line });
                    }
                    Some(Inference::Perpendicular(m)) => {
                        b.constraint(Constraint::Perpendicular { a: m, b: line });
                    }
                    Some(Inference::Tangent(c)) => {
                        b.tangent_at(c, line, p0);
                    }
                    None => {}
                }
                if let Some(op) = self.commit(model, b.edits) {
                    // Chain: the next line starts at this one's end.
                    let end = SketchId {
                        op,
                        index: p1.index,
                    };
                    let at = model
                        .sketch_result(self.feature)
                        .and_then(|r| r.solution.sketch.point(end))
                        .unwrap_or(s.at);
                    self.clicks = vec![Snap {
                        at,
                        target: SnapTarget::Point(end),
                    }];
                }
            }
            Tool::Rectangle => {
                let Some(a) = self.clicks.first().copied() else {
                    self.clicks.push(s);
                    return;
                };
                let (p, q) = (a.at, s.at);
                if (p.x - q.x).abs() <= SAME_POINT || (p.y - q.y).abs() <= SAME_POINT {
                    return;
                }
                let corners = [p, DVec2::new(q.x, p.y), q, DVec2::new(p.x, q.y)];
                let mut b = Batch::default();
                let mut ends = Vec::new();
                for i in 0..4 {
                    let s0 = b.point(corners[i]);
                    let s1 = b.point(corners[(i + 1) % 4]);
                    ends.push((s0, s1));
                }
                let mut lines = Vec::new();
                for (s0, s1) in &ends {
                    lines.push(b.entity(
                        Geometry::Line {
                            start: *s0,
                            end: *s1,
                        },
                        construction,
                    ));
                }
                for i in 0..4 {
                    b.constraint(Constraint::Coincident {
                        a: ends[i].1,
                        b: ends[(i + 1) % 4].0,
                    });
                }
                for (i, line) in lines.iter().enumerate() {
                    b.constraint(if i % 2 == 0 {
                        Constraint::Horizontal { line: *line }
                    } else {
                        Constraint::Vertical { line: *line }
                    });
                }
                b.attach(ends[0].0, &a);
                b.attach(ends[1].1, &s);
                if self.commit(model, b.edits).is_some() {
                    self.reset_tool();
                }
            }
            Tool::CenterRectangle | Tool::Polygon => {
                let Some(first) = self.clicks.first().copied() else {
                    self.clicks.push(s);
                    return;
                };
                if !distinct(&first, &s) {
                    return;
                }
                let mut b = Batch::default();
                let made = if self.tool == Tool::Polygon {
                    shapes::polygon(&mut b, &first, &s, self.polygon_sides, construction)
                } else {
                    shapes::centre_rectangle(&mut b, &first, &s, construction)
                };
                self.commit_shape(model, made.map(|()| b.edits));
            }
            Tool::AlignedRectangle | Tool::Slot => {
                if self.clicks.len() < 2 {
                    if self.clicks.first().is_none_or(|a| distinct(a, &s)) {
                        self.clicks.push(s);
                    }
                    return;
                }
                let (first, second) = (self.clicks[0], self.clicks[1]);
                let mut b = Batch::default();
                let made = if self.tool == Tool::Slot {
                    shapes::slot(&mut b, &first, &second, s.at, construction)
                } else {
                    shapes::aligned_rectangle(&mut b, &first, &second, s.at, construction)
                };
                self.commit_shape(model, made.map(|()| b.edits));
            }
            Tool::TangentArc => {
                let Some(start) = self.tangent_start else {
                    match self.display(model).and_then(|d| tangent_start(&d, &s)) {
                        Some(t) => {
                            self.tangent_start = Some(t);
                            self.clicks.push(s);
                        }
                        None => {
                            self.message = Some(
                                "Start a tangent arc at the end of a line, arc or spline".into(),
                            )
                        }
                    }
                    return;
                };
                if !distinct(&self.clicks[0], &s) {
                    return;
                }
                let mut b = Batch::default();
                let made = shapes::tangent_arc(&mut b, &start, &self.clicks[0], &s, construction);
                self.commit_shape(model, made.map(|()| b.edits));
            }
            Tool::Circle => {
                let Some(c) = self.clicks.first().copied() else {
                    self.clicks.push(s);
                    return;
                };
                let radius = (s.at - c.at).length();
                if radius <= SAME_POINT {
                    return;
                }
                let mut b = Batch::default();
                let center = b.snapped_point(&c);
                let circle = b.entity(Geometry::Circle { center, radius }, construction);
                if let SnapTarget::Point(q) = s.target {
                    b.constraint(Constraint::PointOnCurve {
                        point: q,
                        curve: circle,
                    });
                }
                if self.commit(model, b.edits).is_some() {
                    self.reset_tool();
                }
            }
            Tool::Arc => {
                if self.clicks.len() < 2 {
                    if self.clicks.first().is_none_or(|a| distinct(a, &s)) {
                        self.clicks.push(s);
                    }
                    return;
                }
                let (a, e) = (self.clicks[0], self.clicks[1]);
                let Some((center, _)) = circumcircle(a.at, e.at, s.at) else {
                    self.message = Some("The three points are in a line".into());
                    return;
                };
                let (a0, a1, am) = (
                    (a.at - center).to_angle(),
                    (e.at - center).to_angle(),
                    (s.at - center).to_angle(),
                );
                // The arc runs counter-clockwise from its start to its end
                // through the third point.
                let (first, second) = if ccw(a0, am) < ccw(a0, a1) {
                    (a, e)
                } else {
                    (e, a)
                };
                let mut b = Batch::default();
                let c = b.point(center);
                let p0 = b.snapped_point(&first);
                let p1 = b.snapped_point(&second);
                b.entity(
                    Geometry::Arc {
                        center: c,
                        start: p0,
                        end: p1,
                    },
                    construction,
                );
                if self.commit(model, b.edits).is_some() {
                    self.reset_tool();
                }
            }
            Tool::CenterArc => {
                if self.clicks.len() < 2 {
                    if self.clicks.first().is_none_or(|a| distinct(a, &s)) {
                        self.clicks.push(s);
                        self.turn = 0.0;
                        self.last_angle = None;
                    }
                    return;
                }
                let (c, start) = (self.clicks[0], self.clicks[1]);
                let r = (start.at - c.at).length();
                let dir = (s.at - c.at).normalize_or_zero();
                if dir == DVec2::ZERO {
                    return;
                }
                let end = Snap {
                    at: c.at + dir * r,
                    target: s.target,
                };
                let sweep = self.turn;
                if sweep.abs() < 1e-6 || (end.at - start.at).length() <= SAME_POINT {
                    return;
                }
                let (first, second) = if sweep > 0.0 {
                    (start, end)
                } else {
                    (end, start)
                };
                let mut b = Batch::default();
                let center = b.snapped_point(&c);
                let p0 = b.snapped_point(&first);
                let p1 = b.snapped_point(&second);
                b.entity(
                    Geometry::Arc {
                        center,
                        start: p0,
                        end: p1,
                    },
                    construction,
                );
                if self.commit(model, b.edits).is_some() {
                    self.reset_tool();
                }
            }
            Tool::Spline | Tool::ControlSpline => {
                if self.clicks.last().is_none_or(|a| distinct(a, &s)) {
                    let closes = self.clicks.len() >= 3 && !distinct(&self.clicks[0], &s);
                    self.clicks.push(s);
                    if closes {
                        self.finish(model);
                    }
                }
            }
            Tool::Conic => match self.clicks.len() {
                0 => self.clicks.push(s),
                1 => {
                    if distinct(&self.clicks[0], &s) {
                        self.clicks.push(s);
                    }
                }
                2 => {
                    if shapes::conic_control_ok(self.clicks[0].at, self.clicks[1].at, s.at) {
                        self.clicks.push(s);
                    } else {
                        self.message = Some(shapes::NO_CORNER.into());
                    }
                }
                _ => {
                    let c = &self.clicks;
                    let rho = shapes::conic_rho(c[0].at, c[2].at, c[1].at, s.at);
                    let mut b = Batch::default();
                    let made = shapes::conic(&mut b, &c[0], &c[1], &c[2], rho, construction);
                    self.commit_shape(model, made.map(|()| b.edits));
                }
            },
            Tool::Ellipse | Tool::EllipticalArc => {
                let arc = self.tool == Tool::EllipticalArc;
                match self.clicks.len() {
                    0 => self.clicks.push(s),
                    1 => {
                        if distinct(&self.clicks[0], &s) {
                            self.clicks.push(s);
                        }
                    }
                    2 if !arc => {
                        let mut b = Batch::default();
                        let made = shapes::ellipse(
                            &mut b,
                            &self.clicks[0],
                            &self.clicks[1],
                            s.at,
                            construction,
                        );
                        self.commit_shape(model, made.map(|()| b.edits));
                    }
                    2 => {
                        if shapes::ellipse_axes(self.clicks[0].at, self.clicks[1].at, s.at)
                            .is_some()
                        {
                            self.clicks.push(s);
                        } else {
                            self.message = Some(
                                "The second axis has no length: click away from the first axis"
                                    .into(),
                            );
                        }
                    }
                    3 => {
                        self.clicks.push(s);
                        self.turn = 0.0;
                        self.last_angle = None;
                    }
                    _ => {
                        let mut b = Batch::default();
                        let c = &self.clicks;
                        let made = shapes::elliptical_arc(
                            &mut b,
                            &c[0],
                            &c[1],
                            c[2].at,
                            &c[3],
                            &s,
                            self.turn,
                            construction,
                        );
                        self.commit_shape(model, made.map(|()| b.edits));
                    }
                }
            }
        }
    }

    /// Enter: finishes a spline; ends a line chain. A spline whose last
    /// click is its first point again closes smoothly (ADR 0016).
    pub fn finish(&mut self, model: &mut Model) {
        if matches!(self.tool, Tool::Spline | Tool::ControlSpline) {
            let clicks = std::mem::take(&mut self.clicks);
            let closed =
                clicks.len() >= 4 && !distinct_clicks(&clicks[0], &clicks[clicks.len() - 1]);
            let open = if closed {
                &clicks[..clicks.len() - 1]
            } else {
                &clicks[..]
            };
            let at: Vec<DVec2> = clicks.iter().map(|s| s.at).collect();
            let geometry = |points: Vec<SketchId>| match self.tool {
                Tool::ControlSpline => Some(Geometry::ControlSpline {
                    degree: control_degree(open.len()),
                    control: points,
                    periodic: closed,
                }),
                _ => Some(Geometry::InterpolatedSpline {
                    params: chord_params(&at)?,
                    through: points,
                    periodic: closed,
                }),
            };
            let too_few = match self.tool {
                Tool::ControlSpline => open.len() < 3,
                _ => open.len() < 2,
            };
            if too_few {
                if self.tool == Tool::ControlSpline && !clicks.is_empty() {
                    self.message = Some("A control-point spline needs at least 3 points".into());
                }
                self.reset_tool();
                return;
            }
            let mut b = Batch::default();
            let mut points: Vec<SketchId> = open.iter().map(|s| b.snapped_point(s)).collect();
            if closed {
                // The same point entity: the curve closes on it exactly.
                points.push(points[0]);
            }
            if let Some(g) = geometry(points) {
                b.entity(g, self.construction);
                self.commit(model, b.edits);
            }
        }
        self.reset_tool();
    }

    /// A double-click: finishes a spline or a line chain; on a dimension,
    /// edits its value.
    pub fn double_click(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        match self.tool {
            Tool::Spline | Tool::ControlSpline | Tool::Line => self.finish(model),
            Tool::Select => {
                let Some(d) = self.display(model) else {
                    return;
                };
                if let Some(Pick::Constraint(id)) = pick(&d, &annotations(&d, map), map, pos) {
                    self.edit_dimension(model, id, pos);
                }
            }
            _ => {}
        }
    }

    /// The pointer was pressed at `pos` and has started to move.
    pub fn drag_start(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        self.preview = None;
        if self.tool != Tool::Select {
            self.click(model, map, pos, false);
            self.drag = Some(Drag::Draw { from: pos });
            return;
        }
        let Some(d) = self.display(model) else {
            return;
        };
        let grab = map.to_plane(pos);
        let notes = annotations(&d, map);
        let picked = pick(&d, &notes, map, pos);
        self.drag = match (picked, grab) {
            (Some(Pick::Constraint(id)), Some(grab)) if d.dimension_anchor(id).is_some() => {
                // Grabbing a dimension's label moves it.
                let anchor = d.dimension_anchor(id);
                let label = notes
                    .iter()
                    .find(|a| a.constraint == id)
                    .and_then(|a| map.to_plane(a.rect.center()));
                match (anchor, label) {
                    (Some(anchor), Some(label)) => Some(Drag::Label {
                        dimension: id,
                        grab,
                        from: label - anchor,
                        offset: label - anchor,
                    }),
                    _ => None,
                }
            }
            (Some(Pick::Entity(id)), Some(grab)) => {
                let sketch = &d.sketch;
                match kind(sketch, id) {
                    Some(Kind::Circle) => Some(Drag::Radius { circle: id }),
                    Some(Kind::Conic) => match sketch.entity(id).map(|e| &e.geometry) {
                        Some(&Geometry::Conic { rho, .. }) => Some(Drag::Rho {
                            conic: id,
                            grab,
                            from: rho,
                        }),
                        _ => None,
                    },
                    Some(Kind::Point) => sketch.point(id).map(|at| Drag::Move {
                        grab,
                        points: vec![(id, at)],
                    }),
                    Some(_) => {
                        let ids = sketch
                            .entity(id)
                            .map(|e| e.geometry.referenced_points())
                            .unwrap_or_default();
                        Some(Drag::Move {
                            grab,
                            points: ids
                                .into_iter()
                                .filter_map(|p| Some((p, sketch.point(p)?)))
                                .collect(),
                        })
                    }
                    None => None,
                }
            }
            _ => Some(Drag::Box { from: pos, to: pos }),
        };
    }

    /// The pointer moved to `pos` during a drag.
    pub fn drag_to(&mut self, model: &Model, map: &ScreenMap, pos: Pos2) {
        let Some(drag) = self.drag.as_mut() else {
            return;
        };
        if let Drag::Label {
            grab, from, offset, ..
        } = drag
        {
            if let Some(p) = map.to_plane(pos) {
                *offset = *from + (p - *grab);
            }
            return;
        }
        let base = match &self.preview {
            Some(p) => p.sketch.clone(),
            None => match model.sketch_result(self.feature) {
                Some(r) => r.solution.sketch.clone(),
                None => return,
            },
        };
        let solution = match drag {
            Drag::Box { to, .. } => {
                *to = pos;
                return;
            }
            Drag::Label { .. } => return,
            Drag::Draw { .. } => {
                self.hover(model, map, Some(pos));
                return;
            }
            Drag::Move { grab, points } => {
                let Some(p) = map.to_plane(pos) else {
                    return;
                };
                let delta = p - *grab;
                let options = SolveOptions {
                    drag: points.iter().map(|(id, at)| (*id, *at + delta)).collect(),
                    ..SolveOptions::default()
                };
                solve(&base, &options)
            }
            Drag::Radius { circle } => {
                let Some(p) = map.to_plane(pos) else {
                    return;
                };
                let mut sketch = base;
                let center = match sketch.entity(*circle).map(|e| &e.geometry) {
                    Some(Geometry::Circle { center, .. }) => sketch.point(*center),
                    _ => None,
                };
                let Some(center) = center else {
                    return;
                };
                let r = (p - center).length();
                if r <= SAME_POINT {
                    return;
                }
                for e in &mut sketch.entities {
                    if e.id == *circle {
                        if let Geometry::Circle { radius, .. } = &mut e.geometry {
                            *radius = r;
                        }
                    }
                }
                solve(&sketch, &SolveOptions::default())
            }
            Drag::Rho { conic, grab, from } => {
                let Some(q) = map.to_plane(pos) else {
                    return;
                };
                let mut sketch = base;
                let points = match sketch.entity(*conic).map(|e| &e.geometry) {
                    Some(&Geometry::Conic {
                        start,
                        control,
                        end,
                        ..
                    }) => sketch
                        .point(start)
                        .zip(sketch.point(control))
                        .zip(sketch.point(end)),
                    _ => None,
                };
                let Some(((p0, p1), p2)) = points else {
                    return;
                };
                // The shoulder moves as far along the median as the
                // pointer has, within the tool's range (or the range
                // that already held it).
                let moved = shapes::median_position(p0, p1, p2, q)
                    - shapes::median_position(p0, p1, p2, *grab);
                let (lo, hi) = (shapes::TOOL_RHO.0.min(*from), shapes::TOOL_RHO.1.max(*from));
                let r = *from + moved;
                if !r.is_finite() {
                    return;
                }
                let r = r.clamp(lo, hi);
                for e in &mut sketch.entities {
                    if e.id == *conic {
                        if let Geometry::Conic { rho, .. } = &mut e.geometry {
                            *rho = r;
                        }
                    }
                }
                solve(&sketch, &SolveOptions::default())
            }
        };
        // Keep the last good state if this step does not solve.
        if !matches!(solution.status, SolveStatus::Failed(_)) {
            self.preview = Some(solution);
        }
    }

    /// The pointer was released at `pos`, ending a drag.
    pub fn drag_end(&mut self, model: &mut Model, map: &ScreenMap, pos: Pos2) {
        match self.drag.take() {
            Some(Drag::Move { .. }) | Some(Drag::Radius { .. }) | Some(Drag::Rho { .. }) => {
                let Some(preview) = self.preview.take() else {
                    self.message = Some("That geometry can't move that way".into());
                    return;
                };
                let edits = match model.sketch_def(self.feature) {
                    Some(def) => geometry_edits(&def.sketch, &preview.sketch),
                    None => return,
                };
                if !edits.is_empty() {
                    self.commit(model, edits);
                }
            }
            Some(Drag::Box { from, to }) => {
                let Some(d) = self.display(model) else {
                    return;
                };
                let r = Rect::from_two_pos(from, to);
                let inside = |pts: &[DVec2]| {
                    !pts.is_empty()
                        && pts
                            .iter()
                            .all(|p| map.to_screen(*p).is_some_and(|s| r.contains(s)))
                };
                let mut picked: Vec<SketchId> = d
                    .points()
                    .filter(|(_, at)| inside(&[*at]))
                    .map(|(id, _)| id)
                    .collect();
                picked.extend(
                    d.curves
                        .iter()
                        .filter(|c| inside(&c.polyline))
                        .map(|c| c.id),
                );
                for id in picked {
                    if !self.selection.contains(&id) {
                        self.selection.push(id);
                    }
                }
            }
            Some(Drag::Label {
                dimension,
                from,
                offset,
                ..
            }) => {
                if offset != from {
                    self.commit(
                        model,
                        vec![SketchEdit::PlaceLabel {
                            dimension,
                            offset: Some(offset),
                        }],
                    );
                }
            }
            Some(Drag::Draw { from }) if from.distance(pos) > DRAG_THRESHOLD => {
                self.click(model, map, pos, false);
            }
            Some(Drag::Draw { .. }) | None => {}
        }
        self.preview = None;
    }

    /// The box being dragged out, if any.
    pub fn selection_box(&self) -> Option<Rect> {
        match self.drag {
            Some(Drag::Box { from, to }) => Some(Rect::from_two_pos(from, to)),
            _ => None,
        }
    }

    /// Deletes the selection, with the end points and centres that only the
    /// deleted curves used.
    pub fn delete_selection(&mut self, model: &mut Model) {
        if self.selection.is_empty() {
            return;
        }
        let Some(def) = model.sketch_def(self.feature) else {
            return;
        };
        let selected: BTreeSet<SketchId> = self.selection.iter().copied().collect();
        let items = deletion(&def.sketch, &selected);
        self.selection.clear();
        self.commit(model, vec![SketchEdit::Delete { items }]);
    }

    /// The dimensions whose labels [`SketchEditor::reset_labels`] returns
    /// to automatic placement: the selected ones that were placed, or every
    /// placed one when no dimension is selected.
    pub fn labels_to_reset(&self, model: &Model) -> Vec<SketchId> {
        let Some(def) = model.sketch_def(self.feature) else {
            return Vec::new();
        };
        let is_dimension = |id: &SketchId| {
            def.sketch
                .constraint(*id)
                .is_some_and(|c| c.constraint.value().is_some())
        };
        let placed = def.labels.iter().map(|(id, _)| *id);
        if self.selection.iter().any(is_dimension) {
            placed.filter(|id| self.selection.contains(id)).collect()
        } else {
            placed.collect()
        }
    }

    /// Returns dimension labels to automatic placement (one undo step):
    /// see [`SketchEditor::labels_to_reset`].
    pub fn reset_labels(&mut self, model: &mut Model) {
        let edits: Vec<SketchEdit> = self
            .labels_to_reset(model)
            .into_iter()
            .map(|dimension| SketchEdit::PlaceLabel {
                dimension,
                offset: None,
            })
            .collect();
        if !edits.is_empty() {
            self.commit(model, edits);
        }
    }

    /// Toggles construction on the selected curves.
    pub fn toggle_construction(&mut self, model: &mut Model) {
        let Some(def) = model.sketch_def(self.feature) else {
            return;
        };
        let curves: Vec<&Entity<SketchId>> = self
            .selection
            .iter()
            .filter_map(|id| def.sketch.entity(*id))
            .filter(|e| e.geometry.is_curve())
            .collect();
        if curves.is_empty() {
            self.construction = !self.construction;
            return;
        }
        let to = !curves.iter().all(|e| e.construction);
        let edits = curves
            .iter()
            .map(|e| SketchEdit::SetConstruction {
                entity: e.id,
                construction: to,
            })
            .collect();
        self.commit(model, edits);
    }

    /// Whether a constraint tool applies to the selection.
    pub fn can_constrain(&self, model: &Model, tool: ConstraintTool) -> bool {
        model
            .sketch_def(self.feature)
            .is_some_and(|d| constraints_for(tool, &self.selection, &d.sketch).is_ok())
    }

    /// Applies a constraint tool to the selection.
    pub fn constrain(&mut self, model: &mut Model, tool: ConstraintTool) {
        let Some(d) = self.display(model) else {
            return;
        };
        match constraints_for(tool, &self.selection, &d.sketch) {
            Ok(cs) => {
                let mut b = Batch::default();
                // Curvature continuity implies the tangency it replaces
                // (design note 0010 §5).
                let replaced = replaced_tangencies(&d.sketch, &cs);
                if !replaced.is_empty() {
                    b.edits.push(SketchEdit::Delete {
                        items: replaced.clone(),
                    });
                }
                for c in cs {
                    b.constraint(c);
                }
                if self.commit(model, b.edits).is_some() {
                    self.selection.clear();
                    if !replaced.is_empty() {
                        let solved = self.message.take().unwrap_or_default();
                        let note =
                            "Replaced the tangency there, which curvature continuity implies.";
                        self.message = Some(format!("{note} {solved}").trim_end().to_string());
                    }
                }
            }
            Err(e) => self.message = Some(e),
        }
    }

    /// Starts dimensioning the selection; the value is entered in the
    /// prompt and applied by [`SketchEditor::commit_prompt`].
    pub fn dimension(&mut self, model: &Model, at: Pos2) {
        let Some(d) = self.display(model) else {
            return;
        };
        let made = dimension_for(&self.selection, &d.sketch).and_then(|mut c| {
            let v = measure(&d.sketch, &c).map_err(|e| e.to_string())?;
            c.set_value(v);
            Ok(c)
        });
        match made {
            Ok(constraint) => {
                self.prompt = Some(DimensionPrompt {
                    text: prompt_text(&constraint),
                    constraint,
                    existing: None,
                    at,
                    error: None,
                });
            }
            Err(e) => self.message = Some(e),
        }
    }

    /// Opens the prompt to change an existing dimension's value.
    pub fn edit_dimension(&mut self, model: &Model, id: SketchId, at: Pos2) {
        let Some(c) = model
            .sketch_def(self.feature)
            .and_then(|d| d.sketch.constraint(id))
        else {
            return;
        };
        if c.constraint.value().is_none() {
            return;
        }
        self.prompt = Some(DimensionPrompt {
            text: prompt_text(&c.constraint),
            constraint: c.constraint.clone(),
            existing: Some(id),
            at,
            error: None,
        });
    }

    /// Changes how the prompt's new point-to-point distance is measured,
    /// re-measuring its value.
    pub fn set_prompt_kind(&mut self, model: &Model, kind: DistanceKind) {
        let Some(d) = self.display(model) else {
            return;
        };
        let Some(p) = self.prompt.as_mut() else {
            return;
        };
        if let Constraint::Distance { a, b, .. } = p.constraint {
            let mut c = Constraint::Distance {
                a,
                b,
                kind,
                value: 0.0,
            };
            match measure(&d.sketch, &c) {
                Ok(v) if v > 0.0 => {
                    c.set_value(v);
                    p.text = prompt_text(&c);
                    p.constraint = c;
                    p.error = None;
                }
                Ok(_) => p.error = Some("the points are level in that direction".into()),
                Err(e) => p.error = Some(e.to_string()),
            }
        }
    }

    /// Applies the prompt's value: adds the new dimension or changes the
    /// existing one, in one undo step.
    pub fn commit_prompt(&mut self, model: &mut Model) {
        let Some(mut p) = self.prompt.take() else {
            return;
        };
        let value = match p.value() {
            Ok(v) => v,
            Err(e) => {
                p.error = Some(e);
                self.prompt = Some(p);
                return;
            }
        };
        let edits = match p.existing {
            Some(id) => vec![SketchEdit::SetValue {
                constraint: id,
                value,
            }],
            None => {
                let mut c = p.constraint.clone();
                c.set_value(value);
                let mut b = Batch::default();
                b.constraint(c);
                b.edits
            }
        };
        if self.commit(model, edits).is_some() {
            self.selection.clear();
        }
    }
}

/// `d` rounded to a step of about 1.5 screen points at `ppp` plane units
/// per point: 1, 2 or 5 times a power of ten, and at least one step.
fn round_distance(d: f64, ppp: f64) -> f64 {
    let target = 1.5 * ppp;
    if !(target.is_finite() && target > 0.0) {
        return d;
    }
    let p = 10f64.powf(target.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|k| k * p)
        .find(|s| *s >= target)
        .unwrap_or(10.0 * p);
    ((d / step).round() * step).max(step)
}

// ---------------------------------------------------------------------------
// Region fill
// ---------------------------------------------------------------------------

/// Horizontal spans, `step` points tall, filling the area enclosed by
/// closed loops of screen segments (even-odd rule, so holes stay empty).
pub fn fill_spans(segments: &[[Pos2; 2]], step: f32) -> Vec<Rect> {
    let mut out = Vec::new();
    if segments.is_empty() || step <= 0.0 {
        return out;
    }
    let ymin = segments
        .iter()
        .flatten()
        .map(|p| p.y)
        .fold(f32::INFINITY, f32::min);
    let ymax = segments
        .iter()
        .flatten()
        .map(|p| p.y)
        .fold(f32::NEG_INFINITY, f32::max);
    if !(ymin.is_finite() && ymax.is_finite()) {
        return out;
    }
    let mut y = (ymin / step).floor() * step + step * 0.5;
    let mut xs = Vec::new();
    while y < ymax {
        xs.clear();
        for [a, b] in segments {
            let (a, b) = if a.y <= b.y { (a, b) } else { (b, a) };
            // Half-open in y, so a vertex shared by two segments counts once.
            if y >= a.y && y < b.y {
                xs.push(a.x + (y - a.y) / (b.y - a.y) * (b.x - a.x));
            }
        }
        xs.sort_by(f32::total_cmp);
        for p in xs.chunks_exact(2) {
            out.push(Rect::from_min_max(
                Pos2::new(p[0], y - step * 0.5),
                Pos2::new(p[1], y + step * 0.5),
            ));
        }
        y += step;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circumcircle_of_a_right_triangle() {
        let (c, r) = circumcircle(
            DVec2::new(0.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(0.0, 2.0),
        )
        .unwrap();
        assert!((c - DVec2::new(1.0, 1.0)).length() < 1e-12);
        assert!((r - 2f64.sqrt()).abs() < 1e-12);
        assert!(circumcircle(DVec2::ZERO, DVec2::X, DVec2::X * 2.0).is_none());
    }

    #[test]
    fn fill_spans_cover_a_square_with_a_hole() {
        let sq = |x0: f32, y0: f32, s: f32| {
            let p = [
                Pos2::new(x0, y0),
                Pos2::new(x0 + s, y0),
                Pos2::new(x0 + s, y0 + s),
                Pos2::new(x0, y0 + s),
            ];
            (0..4).map(move |i| [p[i], p[(i + 1) % 4]])
        };
        let segs: Vec<[Pos2; 2]> = sq(0.0, 0.0, 100.0).chain(sq(25.0, 25.0, 50.0)).collect();
        let area: f32 = fill_spans(&segs, 1.0).iter().map(|r| r.area()).sum();
        assert!((area - (10000.0 - 2500.0)).abs() < 1.0, "{area}");
    }

    #[test]
    fn normal_view_faces_each_standard_plane() {
        for frame in [Frame::top(), Frame::front(), Frame::right()] {
            let (yaw, pitch) = normal_view(&frame);
            let cam = Camera {
                yaw,
                pitch,
                ..Camera::isometric()
            };
            assert!((cam.back_dir() - frame.z).length() < 1e-9, "{frame:?}");
            assert!((cam.right_dir() - frame.x).length() < 1e-9, "{frame:?}");
        }
    }

    /// Distance from `p` to the nearest of `lines`.
    fn nearest(lines: &[[DVec2; 2]], p: DVec2) -> f64 {
        lines
            .iter()
            .map(|[a, b]| (p - closest_on_segment(p, *a, *b)).length())
            .fold(f64::INFINITY, f64::min)
    }

    fn every_measure() -> Vec<Measured> {
        let (a, b) = (DVec2::new(-0.02, 0.0), DVec2::new(0.02, 0.01));
        vec![
            Measured::Aligned(a, b),
            Measured::Horizontal(a, b),
            Measured::Vertical(a, b),
            Measured::Foot {
                p: DVec2::new(0.0, 0.02),
                foot: DVec2::ZERO,
                near: DVec2::new(0.01, 0.0),
            },
            Measured::Radial {
                centre: DVec2::new(0.01, 0.0),
                radius: 0.01,
                dir: DVec2::X,
                diameter: false,
            },
            Measured::Radial {
                centre: DVec2::new(0.01, 0.0),
                radius: 0.01,
                dir: DVec2::X,
                diameter: true,
            },
            Measured::Angle {
                vertex: DVec2::ZERO,
                start: 0.0,
                sweep: 0.6,
                a: (DVec2::ZERO, DVec2::new(0.03, 0.0)),
                b: (DVec2::ZERO, DVec2::from_angle(0.6) * 0.03),
            },
            Measured::Across {
                from: DVec2::ZERO,
                to: DVec2::new(0.0, 0.005),
            },
            Measured::Turn {
                centre: DVec2::ZERO,
                source: DVec2::new(0.02, 0.0),
                sweep: -1.2,
            },
        ]
    }

    #[test]
    fn a_placed_label_sits_where_it_was_put_and_its_lines_reach_it() {
        let off = 0.008;
        for m in every_measure() {
            for placed in [
                DVec2::new(0.0, 0.04),
                DVec2::new(0.05, -0.01),
                DVec2::new(-0.03, -0.025),
                m.anchor(),
            ] {
                let (lines, at) = dimension_geometry(m, Some(placed), off, DVec2::ZERO);
                assert_eq!(at, placed, "{m:?}");
                assert!(lines.iter().flatten().all(|p| p.is_finite()), "{m:?}");
                // The label covers part of a line: it is attached (arcs
                // are polylines, within their sag).
                let sag = match m {
                    Measured::Angle { .. } | Measured::Turn { .. } => {
                        (placed - m.anchor()).length() * (1.0 - (TAU / 192.0).cos())
                    }
                    _ => 0.0,
                };
                assert!(nearest(&lines, placed) <= sag + 1e-12, "{m:?} at {placed}");
            }
        }
    }

    #[test]
    fn placed_dimension_lines_follow_the_label() {
        let off = 0.008;
        let (a, b) = (DVec2::new(-0.02, 0.0), DVec2::new(0.02, 0.01));
        // Aligned: parallel to the segment, through the label, reaching
        // past the end the label lies beyond; extension lines start at the
        // measured points.
        let l = DVec2::new(0.06, 0.04);
        let (lines, _) = dimension_geometry(Measured::Aligned(a, b), Some(l), off, DVec2::ZERO);
        let u = (b - a).normalize();
        let [p, q] = lines[2];
        assert!((q - p).normalize().perp_dot(u).abs() < 1e-12);
        assert!((l - p).perp_dot(u).abs() < 1e-12);
        assert!((q - l).length() < 1e-12, "the line ends at the label");
        assert_eq!((lines[0][0], lines[1][0]), (a, b));
        // Horizontal: at the label's height, spanning both points and it.
        let l = DVec2::new(-0.03, -0.02);
        let (lines, _) = dimension_geometry(Measured::Horizontal(a, b), Some(l), off, DVec2::ZERO);
        let [p, q] = lines[2];
        assert_eq!((p.y, q.y), (l.y, l.y));
        assert_eq!((p.x, q.x), (l.x, b.x));
        // The extension lines run down from the points, past the line.
        assert!(lines[0][1].y < l.y && lines[1][1].y < l.y);
        // Vertical: at the label's x.
        let (lines, _) = dimension_geometry(Measured::Vertical(a, b), Some(l), off, DVec2::ZERO);
        assert_eq!((lines[2][0].x, lines[2][1].x), (l.x, l.x));
        assert_eq!((lines[2][0].y, lines[2][1].y), (l.y, b.y));
        // Radius and diameter: along the direction of the label.
        let c = DVec2::new(0.01, 0.0);
        for diameter in [false, true] {
            let m = Measured::Radial {
                centre: c,
                radius: 0.01,
                dir: DVec2::X,
                diameter,
            };
            let l = c + DVec2::from_angle(2.0) * 0.03;
            let (lines, _) = dimension_geometry(m, Some(l), off, DVec2::ZERO);
            let [p, q] = lines[0];
            assert!((q - l).length() < 1e-12);
            let from = if diameter {
                c - DVec2::from_angle(2.0) * 0.01
            } else {
                c
            };
            assert!((p - from).length() < 1e-12, "{diameter}");
            // Inside the circle the line still reaches the circle.
            let (lines, _) = dimension_geometry(m, Some(c + DVec2::Y * 0.004), off, DVec2::ZERO);
            assert!((lines[0][1] - (c + DVec2::Y * 0.01)).length() < 1e-12);
        }
    }

    #[test]
    fn a_placed_angle_label_chooses_the_side_and_widens_the_arc() {
        let off = 0.008;
        let m = Measured::Angle {
            vertex: DVec2::ZERO,
            start: 0.0,
            sweep: 0.6,
            a: (DVec2::new(0.01, 0.0), DVec2::new(0.03, 0.0)),
            b: (DVec2::from_angle(0.6) * 0.01, DVec2::from_angle(0.6) * 0.03),
        };
        let arc_points = |lines: &[[DVec2; 2]], r: f64| -> Vec<DVec2> {
            lines
                .iter()
                .flatten()
                .copied()
                .filter(|p| (p.length() - r).abs() < 1e-12)
                .collect()
        };
        // Inside the angle: the arc spans it, at the label's radius, and
        // extension lines run from the lines' near ends to inside it.
        let l = DVec2::from_angle(0.3) * 0.005;
        let (lines, _) = dimension_geometry(m, Some(l), off, DVec2::ZERO);
        let r = l.length();
        let on_arc = arc_points(&lines, r);
        assert!(on_arc
            .iter()
            .all(|p| (-1e-9..=0.6 + 1e-9).contains(&p.to_angle())));
        // Beyond the lines' ends: extension lines from the far ends.
        let l = DVec2::from_angle(0.3) * 0.05;
        let (lines, _) = dimension_geometry(m, Some(l), off, DVec2::ZERO);
        let from_end = |e: DVec2| lines.iter().any(|[p, _]| (*p - e).length() < 1e-12);
        assert!(from_end(DVec2::new(0.03, 0.0)) && from_end(DVec2::from_angle(0.6) * 0.03));
        // Across the vertex: the opposite (equal) angle.
        let l = DVec2::from_angle(0.3 + PI) * 0.02;
        let (lines, _) = dimension_geometry(m, Some(l), off, DVec2::ZERO);
        let on_arc = arc_points(&lines, 0.02);
        assert!(!on_arc.is_empty());
        assert!(on_arc.iter().all(|p| p.dot(DVec2::from_angle(0.3)) < 0.0));
        // Outside both: the arc widens to reach the label's direction.
        let l = DVec2::from_angle(1.2) * 0.02;
        let (lines, _) = dimension_geometry(m, Some(l), off, DVec2::ZERO);
        let on_arc = arc_points(&lines, 0.02);
        let top = on_arc.iter().map(|p| p.to_angle()).fold(f64::MIN, f64::max);
        assert!((top - 1.2).abs() < 1e-9, "{top}");
    }

    #[test]
    fn automatic_labels_are_unchanged_by_the_placement_code() {
        // The automatic positions the shell tests rely on: outside, `off`
        // from what is measured.
        let off = 0.008;
        let (a, b) = (DVec2::new(-0.02, 0.01), DVec2::new(0.02, 0.01));
        let (_, at) = dimension_geometry(Measured::Aligned(a, b), None, off, DVec2::ZERO);
        assert!((at - DVec2::new(0.0, 0.01 + off)).length() < 1e-15);
        let (_, at) = dimension_geometry(Measured::Horizontal(a, b), None, off, DVec2::ZERO);
        assert!((at - DVec2::new(0.0, 0.01 + off)).length() < 1e-15);
        let (_, at) = dimension_geometry(
            Measured::Vertical(DVec2::new(-0.02, -0.01), DVec2::new(-0.02, 0.01)),
            None,
            off,
            DVec2::ZERO,
        );
        assert!((at - DVec2::new(-0.02 - off, 0.0)).length() < 1e-15);
        let m = Measured::Radial {
            centre: DVec2::ZERO,
            radius: 0.01,
            dir: DVec2::Y,
            diameter: false,
        };
        let (lines, at) = dimension_geometry(m, None, off, DVec2::ZERO);
        assert_eq!(lines, vec![[DVec2::ZERO, DVec2::Y * (0.01 + off)]]);
        assert!((at - DVec2::Y * (0.01 + off * 1.7)).length() < 1e-15);
    }

    proptest::proptest! {
        #[test]
        fn placed_labels_are_attached_and_finite(
            ax in -0.1f64..0.1, ay in -0.1f64..0.1,
            bx in -0.1f64..0.1, by in -0.1f64..0.1,
            lx in -0.2f64..0.2, ly in -0.2f64..0.2,
            sweep in -3.0f64..3.0, start in -4.0f64..4.0,
            r in 1e-4f64..0.05,
        ) {
            let (a, b, l) = (DVec2::new(ax, ay), DVec2::new(bx, by), DVec2::new(lx, ly));
            let off = 0.008;
            let ms = [
                Measured::Aligned(a, b),
                Measured::Horizontal(a, b),
                Measured::Vertical(a, b),
                Measured::Foot { p: a, foot: b, near: b + DVec2::X * 0.01 },
                Measured::Radial { centre: a, radius: r, dir: DVec2::X, diameter: true },
                Measured::Angle {
                    vertex: a,
                    start,
                    sweep,
                    a: (a, a + DVec2::from_angle(start) * r),
                    b: (b, b + DVec2::from_angle(start + sweep) * r),
                },
                Measured::Across { from: a, to: b },
                Measured::Turn { centre: a, source: b, sweep },
            ];
            for m in ms {
                let (lines, at) = dimension_geometry(m, Some(l), off, DVec2::ZERO);
                proptest::prop_assert_eq!(at, l);
                proptest::prop_assert!(lines.iter().flatten().all(|p| p.is_finite()));
                // Arcs are polylines: the label lies on a chord within its
                // sag of the arc.
                let sag = match m {
                    Measured::Angle { .. } | Measured::Turn { .. } => {
                        (l - a).length() * (1.0 - (TAU / 192.0).cos()) + 1e-12
                    }
                    _ => 1e-12,
                };
                proptest::prop_assert!(nearest(&lines, l) <= sag, "{:?}: {}", m, nearest(&lines, l));
            }
        }
    }
}
