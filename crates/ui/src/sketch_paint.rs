//! Painting sketches, sketch-mode feedback and picked regions over the 3D
//! view with egui's painter. Every colour comes from the theme
//! ([`Theme::sketch`], [`Theme::selection`], [`Theme::preselect`]);
//! selected and preselected geometry is drawn in those colours over a
//! translucent halo.

use egui::{Align2, Color32, FontId, Painter, Pos2, Rect, Shape, Stroke, Vec2};
use openshape_core_model::sketch::Geometry;
use openshape_core_model::{DVec2, SketchId};
use openshape_features::EntityState;
use openshape_sketch_solver::SolveStatus;
use openshape_ui_kit::Theme;

use crate::sketch_mode::{
    annotations, sample_resolved, Annotation, CurveEdit, Pick, ScreenMap, SketchDisplay,
    SketchEditor, SnapTarget, Tool,
};
use crate::sketch_ops::trim::Preview;

const LABEL_FONT: f32 = 13.0;

/// Width of the halo under selected and preselected curves, in points.
const HALO: f32 = 7.0;

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
pub fn paint_inactive(painter: &Painter, map: &ScreenMap, d: &SketchDisplay, theme: &Theme) {
    let inactive = theme.sketch.inactive;
    for c in &d.curves {
        let color = if c.construction {
            dimmed(inactive)
        } else {
            inactive
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
pub fn paint_editor(
    painter: &Painter,
    map: &ScreenMap,
    d: &SketchDisplay,
    editor: &SketchEditor,
    theme: &Theme,
) {
    let k = &theme.sketch;
    let selected = |id: SketchId| editor.selection.contains(&id);
    let hovered = |id: SketchId| {
        editor.hover.map(Pick::id) == Some(id)
            || editor.curve_edit.as_ref().map(|e| e.target) == Some(id)
    };
    let color = |id: SketchId| {
        if selected(id) {
            theme.selection
        } else if hovered(id) {
            theme.preselect
        } else if d.conflicted.contains(&id) {
            k.conflict
        } else if d.states.get(&id) == Some(&EntityState::Constrained) {
            k.constrained
        } else {
            k.free
        }
    };
    // Halos first, so every curve and point is drawn over them.
    for c in d.curves.iter().filter(|c| selected(c.id) || hovered(c.id)) {
        stroke_polyline(painter, map, &c.polyline, Stroke::new(HALO, k.halo), false);
    }
    for (_, at) in d.points().filter(|(id, _)| selected(*id) || hovered(*id)) {
        if let Some(s) = map.to_screen(at) {
            painter.circle_filled(s, HALO, k.halo);
        }
    }
    // Plane origin and axes.
    if let (Some(o), Some(x), Some(y)) = (
        map.to_screen(DVec2::ZERO),
        map.to_screen(DVec2::X),
        map.to_screen(DVec2::Y),
    ) {
        let axis = |to: Pos2| (to - o).normalized() * 28.0;
        painter.line_segment([o, o + axis(x)], Stroke::new(1.5, k.axis_x));
        painter.line_segment([o, o + axis(y)], Stroke::new(1.5, k.axis_y));
        painter.circle_filled(o, 3.0, k.constrained);
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
    // Control-point splines show their control polygons.
    for e in &d.sketch.entities {
        if let Geometry::ControlSpline { control, .. } = &e.geometry {
            let pts: Vec<DVec2> = control.iter().filter_map(|p| d.sketch.point(*p)).collect();
            stroke_polyline(
                painter,
                map,
                &pts,
                Stroke::new(1.0, dimmed(color(e.id))),
                true,
            );
        }
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
        paint_annotation(painter, &a, d, editor, theme);
    }
    // What trim, extend or split would do.
    if let Some(CurveEdit { plan: Ok(plan), .. }) = &editor.curve_edit {
        match &plan.preview {
            Preview::Piece(c) => {
                let pts = sample_resolved(c);
                let (color, dashed) = match editor.tool {
                    Tool::Trim => (theme.danger, false),
                    _ => (theme.preselect, true),
                };
                stroke_polyline(painter, map, &pts, Stroke::new(HALO, k.halo), false);
                stroke_polyline(painter, map, &pts, Stroke::new(2.5, color), dashed);
            }
            Preview::Point(q) => {
                if let Some(p) = map.to_screen(*q) {
                    painter.circle_filled(p, HALO, k.halo);
                    painter.circle_stroke(p, 5.0, Stroke::new(2.0, theme.selection));
                }
            }
            Preview::Shapes { curves, points } => {
                for c in curves {
                    let pts = sample_resolved(c);
                    stroke_polyline(painter, map, &pts, Stroke::new(HALO, k.halo), false);
                    stroke_polyline(painter, map, &pts, Stroke::new(2.0, theme.preselect), true);
                }
                for q in points {
                    if let Some(p) = map.to_screen(*q) {
                        painter.rect_filled(
                            Rect::from_center_size(p, Vec2::splat(6.0)),
                            0.0,
                            theme.preselect,
                        );
                    }
                }
            }
        }
    }
    // The tool in progress.
    let band = Stroke::new(1.5, k.free);
    for pts in editor.rubber_band() {
        stroke_polyline(painter, map, &pts, band, true);
    }
    for s in &editor.clicks {
        if let Some(p) = map.to_screen(s.at) {
            painter.rect_filled(Rect::from_center_size(p, Vec2::splat(5.0)), 0.0, k.free);
        }
    }
    if let Some(cursor) = editor.cursor {
        if let Some(p) = map.to_screen(cursor.at) {
            let ring = Stroke::new(1.5, theme.selection);
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
                    painter.circle_filled(p, 2.0, k.free);
                }
            }
            if let Some(h) = editor.inferred {
                painter.text(
                    p + Vec2::new(12.0, 12.0),
                    Align2::LEFT_TOP,
                    h.glyph(),
                    FontId::proportional(LABEL_FONT),
                    k.free,
                );
            }
        }
    }
    if let Some(r) = editor.selection_box() {
        painter.rect_filled(r, 0.0, k.halo.gamma_multiply(0.5));
        painter.rect_stroke(
            r,
            0.0,
            Stroke::new(1.0, theme.selection),
            egui::StrokeKind::Inside,
        );
    }
}

fn paint_annotation(
    painter: &Painter,
    a: &Annotation,
    d: &SketchDisplay,
    editor: &SketchEditor,
    theme: &Theme,
) {
    let k = &theme.sketch;
    let id = a.constraint;
    let color = if editor.selection.contains(&id) {
        theme.selection
    } else if editor.hover == Some(Pick::Constraint(id)) {
        theme.preselect
    } else if d.conflicting.contains(&id) {
        k.conflict
    } else if d.redundant.contains(&id) {
        k.redundant
    } else if a.dimension {
        k.dimension
    } else {
        k.glyph
    };
    for l in &a.lines {
        painter.line_segment(*l, Stroke::new(1.0, color));
    }
    painter.rect_filled(a.rect, theme.rounding, k.label_fill);
    painter.rect_stroke(
        a.rect,
        theme.rounding,
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

/// The colour of a sketch's solve state.
pub fn status_color(d: &SketchDisplay, theme: &Theme) -> Color32 {
    let k = &theme.sketch;
    match d.status {
        SolveStatus::Failed(_) => k.conflict,
        SolveStatus::Redundant => k.redundant,
        SolveStatus::Solved if d.dof == 0 => k.constrained,
        SolveStatus::Solved => k.free,
    }
}
