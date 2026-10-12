//! Sketch text (ADR 0015): placing a line of text in a box of four
//! construction lines, and changing a text's characters or font.
//!
//! The layout comes from `opendatum-text`; the document stores its outline
//! in box units, so these functions only lay out when someone types and
//! never when a sketch regenerates.

use std::collections::BTreeSet;

use glam::DVec2;
use opendatum_core_model::sketch::{Constraint, Geometry, TextPiece};
use opendatum_core_model::{SketchData, SketchEdit, SketchId};
use opendatum_sketch_solver::text_curves;
use opendatum_text::{lay_out, Font, Layout};

use super::{Batch, SAME_POINT};
use crate::sketch_mode::{sample_resolved, Snap};

/// The text the tool starts with.
pub const DEFAULT_TEXT: &str = "Text";

/// Why a constraint or dimension cannot take a text itself.
pub const CONSTRAIN_THE_BOX: &str =
    "Text is constrained through its box: select its corners or box lines instead";

/// Why trim, extend, split and fillet leave text alone.
pub const NOT_EDITABLE: &str =
    "Text cannot be trimmed, extended, split or filleted; change it in the dialog or delete it";

/// The outline of a layout as the document stores it.
pub fn outline(layout: &Layout) -> Vec<Vec<TextPiece>> {
    layout
        .contours
        .iter()
        .map(|c| {
            c.iter()
                .map(|p| TextPiece {
                    degree: p.degree,
                    points: p.points.clone(),
                })
                .collect()
        })
        .collect()
}

/// Lays out `content` in `font`, or says why it cannot be.
pub fn layout(content: &str, font: Font) -> Result<Layout, String> {
    lay_out(content, font).map_err(|e| e.to_string())
}

/// The box of a text of width over height `aspect` with one corner at
/// `first`, reaching towards `pointer`: its height is how far the pointer
/// is above or below `first`, and it reads left to right with a horizontal
/// bottom. Returns the corners `[bottom-left, bottom-right, top-right,
/// top-left]` and the index of the one at `first`, or `None` while the
/// pointer is level with `first`.
pub fn text_box(first: DVec2, pointer: DVec2, aspect: f64) -> Option<([DVec2; 4], usize)> {
    let dy = pointer.y - first.y;
    if dy.abs() <= SAME_POINT {
        return None;
    }
    let (h, w) = (dy.abs(), dy.abs() * aspect);
    let left = pointer.x < first.x;
    let below = dy < 0.0;
    let bl = first - DVec2::new(if left { w } else { 0.0 }, if below { h } else { 0.0 });
    let corners = [
        bl,
        bl + DVec2::new(w, 0.0),
        bl + DVec2::new(w, h),
        bl + DVec2::new(0.0, h),
    ];
    let at_first = match (left, below) {
        (false, false) => 0,
        (true, false) => 1,
        (true, true) => 2,
        (false, true) => 3,
    };
    Some((corners, at_first))
}

/// The text tool's preview: the box and each contour of the letters, as
/// plane polylines.
pub fn preview(first: DVec2, pointer: DVec2, layout: &Layout) -> Vec<Vec<DVec2>> {
    let Some((c, _)) = text_box(first, pointer, layout.aspect) else {
        return vec![vec![first, pointer]];
    };
    let mut out = vec![vec![c[0], c[1], c[2], c[3], c[0]]];
    let Some(curves) = text_curves(c, layout.aspect, &outline(layout)) else {
        return out;
    };
    // One polyline per contour, so the dashes run on round each letter.
    let mut pieces = curves.iter();
    for contour in &layout.contours {
        let mut line: Vec<DVec2> = Vec::new();
        for curve in pieces.by_ref().take(contour.len()) {
            let pts = sample_resolved(curve);
            let skip = usize::from(!line.is_empty());
            line.extend(pts.into_iter().skip(skip));
        }
        out.push(line);
    }
    out
}

/// Adds a text: its four corners, the text on them, the box as four
/// construction lines (bottom, right, top, left) with a horizontal bottom,
/// and the corner at the first click attached to what it snapped to.
/// Three freedoms remain: where it is and how tall.
pub fn place(
    b: &mut Batch,
    first: &Snap,
    pointer: DVec2,
    content: &str,
    font: Font,
    construction: bool,
) -> Result<(), String> {
    let layout = layout(content, font).map_err(|e| format!("Text not placed: {e}"))?;
    let (at, corner) = text_box(first.at, pointer, layout.aspect)
        .ok_or("Move the pointer above or below the first click to set the text's height")?;
    let corners = at.map(|p| b.point(p));
    b.entity(
        Geometry::Text {
            corners,
            content: content.to_string(),
            font: font.name().to_string(),
            aspect: layout.aspect,
            mirrored: false,
            outline: outline(&layout),
        },
        construction,
    );
    let lines: Vec<SketchId> = (0..4)
        .map(|i| {
            b.entity(
                Geometry::Line {
                    start: corners[i],
                    end: corners[(i + 1) % 4],
                },
                true,
            )
        })
        .collect();
    b.constraint(Constraint::Horizontal { line: lines[0] });
    b.attach(corners[corner], first);
    Ok(())
}

/// A text's characters and font name, if `id` is a text.
pub fn text_of(sketch: &SketchData, id: SketchId) -> Option<(String, String)> {
    match &sketch.entity(id)?.geometry {
        Geometry::Text { content, font, .. } => Some((content.clone(), font.clone())),
        _ => None,
    }
}

/// The lines joining consecutive corners of `text`'s box (the box the
/// tool draws), if `text` is a text.
pub fn box_lines(sketch: &SketchData, text: SketchId) -> Vec<SketchId> {
    let Some(Geometry::Text { corners: c, .. }) = sketch.entity(text).map(|e| &e.geometry) else {
        return Vec::new();
    };
    let sides: Vec<(SketchId, SketchId)> = (0..4).map(|i| (c[i], c[(i + 1) % 4])).collect();
    sketch
        .entities
        .iter()
        .filter(|e| {
            matches!(e.geometry, Geometry::Line { start, end }
                if sides.contains(&(start, end)) || sides.contains(&(end, start)))
        })
        .map(|e| e.id)
        .collect()
}

/// Every text of the sketch.
pub fn texts(sketch: &SketchData) -> BTreeSet<SketchId> {
    sketch
        .entities
        .iter()
        .filter(|e| matches!(e.geometry, Geometry::Text { .. }))
        .map(|e| e.id)
        .collect()
}

/// True when `id` is a text.
pub fn is_text(sketch: &SketchData, id: SketchId) -> bool {
    text_of(sketch, id).is_some()
}

/// `text` and the texts tied to it as mirror images or pattern copies (of
/// it, or of each other), which must keep its characters and font.
pub fn family(sketch: &SketchData, text: SketchId) -> BTreeSet<SketchId> {
    let mut pairs: Vec<(SketchId, SketchId)> = Vec::new();
    for c in &sketch.constraints {
        match &c.constraint {
            Constraint::Mirror { pairs: p, .. } => pairs.extend(p.iter().copied()),
            Constraint::LinearPattern { copies, .. }
            | Constraint::CircularPattern { copies, .. } => {
                pairs.extend(copies.iter().map(|(s, i, _)| (*s, *i)));
            }
            _ => {}
        }
    }
    let mut out = BTreeSet::from([text]);
    loop {
        let before = out.len();
        for (a, b) in &pairs {
            if out.contains(a) || out.contains(b) {
                out.insert(*a);
                out.insert(*b);
            }
        }
        if out.len() == before {
            return out;
        }
    }
}

/// The edits that give `text` (and its [`family`]) new characters or a new
/// font: the new outline, and the box's right-hand corners moved so that
/// the box keeps its bottom-left corner, height and direction and takes
/// the new width.
pub fn set_text(
    sketch: &SketchData,
    text: SketchId,
    content: &str,
    font: Font,
) -> Result<Vec<SketchEdit>, String> {
    if !is_text(sketch, text) {
        return Err("Select a text to change".into());
    }
    let layout = layout(content, font)?;
    let shape = outline(&layout);
    let mut edits = Vec::new();
    let mut moves = Vec::new();
    for id in family(sketch, text) {
        let Some(Geometry::Text {
            corners, mirrored, ..
        }) = sketch.entity(id).map(|e| &e.geometry)
        else {
            continue;
        };
        edits.push(SketchEdit::SetText {
            entity: id,
            content: content.to_string(),
            font: font.name().to_string(),
            aspect: layout.aspect,
            outline: shape.clone(),
        });
        let (Some(bl), Some(tl)) = (sketch.point(corners[0]), sketch.point(corners[3])) else {
            continue;
        };
        // br − bl = k·R(tl − bl), R a quarter turn clockwise, k negative
        // when mirrored (ADR 0015).
        let up = tl - bl;
        let k = if *mirrored {
            -layout.aspect
        } else {
            layout.aspect
        };
        let across = DVec2::new(up.y, -up.x) * k;
        moves.push((corners[1], bl + across));
        moves.push((corners[2], tl + across));
    }
    edits.push(SketchEdit::MovePoints { points: moves });
    Ok(edits)
}

/// The fonts offered, by name.
pub fn font_names() -> impl Iterator<Item = (Font, &'static str)> {
    Font::ALL.into_iter().map(|f| (f, f.name()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_box_reaches_towards_the_pointer_and_reads_left_to_right() {
        let f = DVec2::new(1.0, 2.0);
        for (p, at) in [
            (DVec2::new(5.0, 3.0), 0),
            (DVec2::new(-5.0, 3.0), 1),
            (DVec2::new(-5.0, 1.0), 2),
            (DVec2::new(5.0, 1.0), 3),
        ] {
            let (c, i) = text_box(f, p, 2.5).unwrap();
            assert_eq!(i, at);
            assert_eq!(c[i], f);
            assert!((c[1] - c[0] - DVec2::new(2.5, 0.0)).length() < 1e-15);
            assert!((c[3] - c[0] - DVec2::new(0.0, 1.0)).length() < 1e-15);
            assert_eq!(c[2] - c[3], c[1] - c[0]);
        }
        assert!(text_box(f, DVec2::new(9.0, 2.0), 2.5).is_none());
    }
}
