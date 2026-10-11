//! Painting sketches, sketch-mode feedback and picked regions over the 3D
//! view with egui's painter.

use egui::{Align2, Color32, FontId, Painter, Pos2, Rect, Shape, Stroke, Vec2};
use openshape_core_model::{DVec2, SketchId};
use openshape_features::EntityState;
use openshape_sketch_solver::SolveStatus;

use crate::sketch_mode::{
    annotations, Annotation, Pick, ScreenMap, SketchDisplay, SketchEditor, SnapTarget,
};

/// Fully constrained geometry.
pub const CONSTRAINED: Color32 = Color32::from_rgb(18, 20, 26);
/// Geometry that can still move.
pub const FREE: Color32 = Color32::from_rgb(25, 95, 225);
/// Geometry or constraints in conflict.
pub const CONFLICT: Color32 = Color32::from_rgb(215, 35, 35);
/// Redundant constraints.
pub const REDUNDANT: Color32 = Color32::from_rgb(220, 120, 0);
/// Selected items.
pub const SELECTED: Color32 = Color32::from_rgb(255, 135, 0);
/// The item under the pointer.
pub const HOVERED: Color32 = Color32::from_rgb(255, 185, 70);
/// Sketches that are not being edited.
pub const INACTIVE: Color32 = Color32::from_rgb(45, 65, 120);
/// Picked profile regions.
pub const REGION_FILL: Color32 = Color32::from_rgba_premultiplied(30, 80, 170, 90);
/// The region under the pointer while picking.
pub const REGION_HOVER: Color32 = Color32::from_rgba_premultiplied(60, 110, 190, 45);

const LABEL_FONT: f32 = 13.0;

fn dimmed(c: Color32) -> Color32 {
    c.gamma_multiply(0.55)
}

fn stroke_polyline(
    painter: &Painter,
    map: &ScreenMap,
    pts: &[DVec2],
    stroke: Stroke,
    dashed: bool,
) {
    for run in map.polyline(pts) {
        if dashed {
            painter.extend(Shape::dashed_line(&run, stroke, 6.0, 4.0));
        } else {
            painter.add(Shape::line(run, stroke));
        }
    }
}

/// Draws a sketch that is not being edited.
pub fn paint_inactive(painter: &Painter, map: &ScreenMap, d: &SketchDisplay) {
    for c in &d.curves {
        let color = if c.construction {
            dimmed(INACTIVE)
        } else {
            INACTIVE
        };
        stroke_polyline(
            painter,
            map,
            &c.polyline,
            Stroke::new(1.25, color),
            c.construction,
        );
    }
}

/// Draws the sketch being edited with its diagnostics, constraint glyphs,
/// dimensions and the editor's feedback.
pub fn paint_editor(painter: &Painter, map: &ScreenMap, d: &SketchDisplay, editor: &SketchEditor) {
    let selected = |id: SketchId| editor.selection.contains(&id);
    let hovered = |id: SketchId| editor.hover.map(Pick::id) == Some(id);
    let color = |id: SketchId| {
        if selected(id) {
            SELECTED
        } else if hovered(id) {
            HOVERED
        } else if d.conflicted.contains(&id) {
            CONFLICT
        } else if d.states.get(&id) == Some(&EntityState::Constrained) {
            CONSTRAINED
        } else {
            FREE
        }
    };
    // Plane origin and axes.
    if let (Some(o), Some(x), Some(y)) = (
        map.to_screen(DVec2::ZERO),
        map.to_screen(DVec2::X),
        map.to_screen(DVec2::Y),
    ) {
        let axis = |to: Pos2| (to - o).normalized() * 28.0;
        painter.line_segment(
            [o, o + axis(x)],
            Stroke::new(1.5, Color32::from_rgb(200, 50, 50)),
        );
        painter.line_segment(
            [o, o + axis(y)],
            Stroke::new(1.5, Color32::from_rgb(40, 150, 60)),
        );
        painter.circle_filled(o, 3.0, CONSTRAINED);
    }
    for c in &d.curves {
        let col = color(c.id);
        let col = if c.construction { dimmed(col) } else { col };
        let width = if selected(c.id) || hovered(c.id) {
            2.5
        } else {
            1.75
        };
        stroke_polyline(
            painter,
            map,
            &c.polyline,
            Stroke::new(width, col),
            c.construction,
        );
    }
    for (id, at) in d.points() {
        if let Some(s) = map.to_screen(at) {
            let size = if selected(id) || hovered(id) {
                7.0
            } else {
                5.0
            };
            painter.rect_filled(Rect::from_center_size(s, Vec2::splat(size)), 0.0, color(id));
        }
    }
    for a in annotations(d, map) {
        paint_annotation(painter, &a, d, editor);
    }
    // The tool in progress.
    let band = Stroke::new(1.5, FREE);
    for pts in editor.rubber_band() {
        stroke_polyline(painter, map, &pts, band, true);
    }
    for s in &editor.clicks {
        if let Some(p) = map.to_screen(s.at) {
            painter.rect_filled(Rect::from_center_size(p, Vec2::splat(5.0)), 0.0, FREE);
        }
    }
    if let Some(cursor) = editor.cursor {
        if let Some(p) = map.to_screen(cursor.at) {
            let ring = Stroke::new(1.5, SELECTED);
            match cursor.target {
                SnapTarget::Point(_) | SnapTarget::Origin => {
                    painter.circle_stroke(p, 7.0, ring);
                }
                SnapTarget::Curve(_) => {
                    let k = 5.0;
                    painter.line_segment([p - Vec2::splat(k), p + Vec2::splat(k)], ring);
                    painter.line_segment([p + Vec2::new(-k, k), p + Vec2::new(k, -k)], ring);
                }
                SnapTarget::Free => {
                    painter.circle_filled(p, 2.0, FREE);
                }
            }
            if let Some(h) = editor.inferred {
                painter.text(
                    p + Vec2::new(12.0, 12.0),
                    Align2::LEFT_TOP,
                    if h { "H" } else { "V" },
                    FontId::proportional(LABEL_FONT),
                    FREE,
                );
            }
        }
    }
    if let Some(r) = editor.selection_box() {
        painter.rect_filled(r, 0.0, Color32::from_rgba_unmultiplied(255, 135, 0, 30));
        painter.rect_stroke(r, 0.0, Stroke::new(1.0, SELECTED), egui::StrokeKind::Inside);
    }
}

fn paint_annotation(painter: &Painter, a: &Annotation, d: &SketchDisplay, editor: &SketchEditor) {
    let id = a.constraint;
    let color = if editor.selection.contains(&id) {
        SELECTED
    } else if editor.hover == Some(Pick::Constraint(id)) {
        HOVERED
    } else if d.conflicting.contains(&id) {
        CONFLICT
    } else if d.redundant.contains(&id) {
        REDUNDANT
    } else if a.dimension {
        Color32::from_rgb(20, 40, 110)
    } else {
        Color32::from_rgb(60, 60, 70)
    };
    for l in &a.lines {
        painter.line_segment(*l, Stroke::new(1.0, color));
    }
    painter.rect_filled(
        a.rect,
        3.0,
        Color32::from_rgba_unmultiplied(255, 255, 255, 215),
    );
    painter.rect_stroke(
        a.rect,
        3.0,
        Stroke::new(1.0, color),
        egui::StrokeKind::Inside,
    );
    painter.text(
        a.rect.center(),
        Align2::CENTER_CENTER,
        &a.text,
        FontId::proportional(LABEL_FONT),
        color,
    );
}

/// Fills regions given as boundary polylines (plane coordinates).
pub fn paint_region(painter: &Painter, map: &ScreenMap, outline: &[Vec<DVec2>], fill: Color32) {
    let mut segments = Vec::new();
    for pts in outline {
        let screen: Vec<Option<Pos2>> = pts.iter().map(|p| map.to_screen(*p)).collect();
        for w in screen.windows(2) {
            if let (Some(a), Some(b)) = (w[0], w[1]) {
                segments.push([a, b]);
            }
        }
    }
    let mut mesh = egui::Mesh::default();
    for r in crate::sketch_mode::fill_spans(&segments, 2.0) {
        mesh.add_colored_rect(r, fill);
    }
    painter.add(Shape::mesh(mesh));
    let edge = Stroke::new(1.5, fill.to_opaque());
    for s in segments {
        painter.line_segment(s, edge);
    }
}

/// A one-line summary of a sketch's solve state, with its colour.
pub fn status_color(d: &SketchDisplay) -> Color32 {
    match d.status {
        SolveStatus::Failed(_) => CONFLICT,
        SolveStatus::Redundant => REDUNDANT,
        SolveStatus::Solved if d.dof == 0 => CONSTRAINED,
        SolveStatus::Solved => FREE,
    }
}
