//! Sketch text (ADR 0015), driven through the editor: placing a line of
//! text in a box that grows towards the pointer, three freedoms left, a
//! corner held where the first click snapped; extruding letters; editing
//! the characters and font of a text and of its mirror image and pattern
//! copies; deleting a text with its box; and what other tools do with
//! letters (no snapping onto them, no constraints, cuts or fillets on
//! them).

mod common;

use common::*;
use openshape_app::sketch_mode::{
    constraints_for, dimension_for, ConstraintTool, CurveEdit, SnapTarget, Tool,
};
use openshape_app::sketch_ops::text::{self, CONSTRAIN_THE_BOX, NOT_EDITABLE};
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{DVec2, FeatureId, SketchId};
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};
use openshape_text::{lay_out, Font};

/// The texts of a sketch: (id, corners, content, font, aspect, mirrored).
#[allow(clippy::type_complexity)]
fn texts(m: &Model, f: FeatureId) -> Vec<(SketchId, [SketchId; 4], String, String, f64, bool)> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter_map(|e| match &e.geometry {
            Geometry::Text {
                corners,
                content,
                font,
                aspect,
                mirrored,
                ..
            } => Some((
                e.id,
                *corners,
                content.clone(),
                font.clone(),
                *aspect,
                *mirrored,
            )),
            _ => None,
        })
        .collect()
}

/// How close clicked geometry lands (millimetres): screen positions are
/// `f32`.
const NEAR: f64 = 1e-4;

/// Solved position of a point.
fn pos(m: &Model, f: FeatureId, id: SketchId) -> DVec2 {
    m.sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .point(id)
        .expect("a point")
}

/// Places `content` in `font` with the Text tool: first click at `a`, the
/// height from the pointer at `b` (millimetres). Returns the text.
fn place(
    m: &mut Model,
    ed: &mut openshape_app::sketch_mode::SketchEditor,
    map: &openshape_app::sketch_mode::ScreenMap,
    content: &str,
    font: Font,
    a: (f64, f64),
    b: (f64, f64),
) -> SketchId {
    let f = ed.feature;
    let before = texts(m, f).len();
    ed.set_tool(Tool::Text);
    ed.text_content = content.into();
    ed.text_font = font;
    click(m, ed, map, a.0, a.1);
    assert_eq!(ed.clicks.len(), 1, "{:?}", ed.message);
    hover(m, ed, map, b.0, b.1);
    // The preview shows the box and each contour of the letters, closed.
    let shown = ed.rubber_band();
    let contours = lay_out(content, font).unwrap().contours.len();
    assert_eq!(shown.len(), 1 + contours);
    assert!(shown.iter().all(|l| l.first() == l.last()));
    click(m, ed, map, b.0, b.1);
    assert!(ed.clicks.is_empty(), "{:?}", ed.message);
    let all = texts(m, f);
    assert_eq!(all.len(), before + 1);
    ed.set_tool(Tool::Select);
    all.last().unwrap().0
}

fn corners_mm(m: &Model, f: FeatureId, c: [SketchId; 4]) -> [DVec2; 4] {
    c.map(|k| pos(m, f, k) / MM)
}

/// The box-unit area of a layout whose pieces are all lines.
fn polygon_area(content: &str, font: Font) -> f64 {
    let l = lay_out(content, font).unwrap();
    l.contours
        .iter()
        .map(|c| {
            c.iter()
                .map(|p| {
                    assert_eq!(p.degree, 1, "a straight-sided letter");
                    p.points[0].perp_dot(p.points[1])
                })
                .sum::<f64>()
                * 0.5
        })
        .sum::<f64>()
        .abs()
}

#[test]
fn text_is_placed_in_a_box_with_three_freedoms() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Text);
    assert_eq!(ed.hint(), Tool::Text.hint(0));
    let t = place(
        &mut m,
        &mut ed,
        &map,
        "Hi",
        Font::Sans,
        (10.0, 10.0),
        (30.0, 18.0),
    );
    let (_, c, content, font, aspect, mirrored) = texts(&m, f).remove(0);
    assert_eq!(
        (content.as_str(), font.as_str(), mirrored),
        ("Hi", "Liberation Sans", false)
    );
    assert_eq!(aspect, lay_out("Hi", Font::Sans).unwrap().aspect);
    // The first click is the bottom-left corner; the pointer set the height.
    let k = corners_mm(&m, f, c);
    let w = 8.0 * aspect;
    for (got, want) in k.iter().zip([
        (10.0, 10.0),
        (10.0 + w, 10.0),
        (10.0 + w, 18.0),
        (10.0, 18.0),
    ]) {
        assert!((*got - DVec2::new(want.0, want.1)).length() < NEAR, "{k:?}");
    }
    // Four construction lines on the corners, the bottom one horizontal.
    assert_eq!(construction(&m, f), 4);
    assert_eq!(
        text::box_lines(&m.sketch_def(f).unwrap().sketch, t).len(),
        4
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Horizontal { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 3);
    // Its letters are curves of the text: "H" one region, "i" two (stem
    // and dot), each a region of the sketch.
    let letters = curves(&m, f, false);
    assert!(letters.iter().all(|(id, _)| *id == t));
    assert_eq!(region_areas(&m, f).len(), 3);
}

#[test]
fn the_box_grows_towards_the_pointer() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    place(
        &mut m,
        &mut ed,
        &map,
        "Ab",
        Font::Serif,
        (10.0, 10.0),
        (-5.0, 4.0),
    );
    let (_, c, _, _, aspect, _) = texts(&m, f).remove(0);
    let k = corners_mm(&m, f, c);
    // Below and to the left: the first click is the top-right corner.
    assert!((k[2] - DVec2::new(10.0, 10.0)).length() < NEAR, "{k:?}");
    assert!(
        (k[0] - DVec2::new(10.0 - 6.0 * aspect, 4.0)).length() < NEAR,
        "{k:?}"
    );
}

#[test]
fn a_first_click_on_a_point_holds_that_corner() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-20.0, 0.0), (0.0, 0.0)],
    );
    let before = dof(&m, f);
    place(
        &mut m,
        &mut ed,
        &map,
        "X",
        Font::Mono,
        (0.0, 0.0),
        (8.0, 6.0),
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Coincident { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before + 1);
}

#[test]
fn letters_extrude_to_their_exact_shapes() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    place(
        &mut m,
        &mut ed,
        &map,
        "I",
        Font::SansBold,
        (0.0, 0.0),
        (5.0, 20.0),
    );
    let (_, c, ..) = texts(&m, f).remove(0);
    // The box's height as solved (metres); the letter scales with it.
    let h = (pos(&m, f, c[3]) - pos(&m, f, c[0])).length();
    let area = polygon_area("I", Font::SansBold) * h * h;
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1);
    assert!(
        (areas[0] - area).abs() <= 1e-9 * area,
        "{areas:?} vs {area}"
    );
    // The stem's middle, in box units (the middle of its bounding box).
    let l = lay_out("I", Font::SansBold).unwrap();
    let pts: Vec<DVec2> = l.contours[0].iter().map(|p| p.points[0]).collect();
    let mid = pts.iter().fold(DVec2::ZERO, |a, p| a + *p) / pts.len() as f64;
    let inside = mid * 20.0;
    let x = extrude(&mut m, &map, &[(inside.x, inside.y)], "3 mm");
    let v = volume(&m, x);
    assert!((v - area * 3.0 * MM).abs() <= 1e-9 * v, "{v}");
}

#[test]
fn editing_a_text_keeps_its_corner_height_and_direction() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let t = place(
        &mut m,
        &mut ed,
        &map,
        "A",
        Font::Sans,
        (10.0, 10.0),
        (20.0, 22.0),
    );
    ed.selection = vec![t];
    ed.sync_text_edit(&m);
    let edit = ed
        .text_edit
        .as_mut()
        .expect("the selected text is editable");
    assert_eq!((edit.content.as_str(), edit.font), ("A", Font::Sans));
    edit.content = "AB 12".into();
    edit.font = Font::SerifBold;
    ed.apply_text(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Text changed"));
    let (_, c, content, font, aspect, _) = texts(&m, f).remove(0);
    assert_eq!(
        (content.as_str(), font.as_str()),
        ("AB 12", "Liberation Serif Bold")
    );
    assert_eq!(aspect, lay_out("AB 12", Font::SerifBold).unwrap().aspect);
    let k = corners_mm(&m, f, c);
    assert!((k[0] - DVec2::new(10.0, 10.0)).length() < NEAR, "{k:?}");
    assert!((k[3] - DVec2::new(10.0, 22.0)).length() < NEAR, "{k:?}");
    assert!(
        (k[1] - DVec2::new(10.0 + 12.0 * aspect, 10.0)).length() < NEAR,
        "{k:?}"
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 3);
    // Undo restores it, and the dialog shows what the document has.
    assert!(m.undo());
    ed.sync_text_edit(&m);
    assert_eq!(texts(&m, f)[0].2, "A");
    assert_eq!(ed.text_edit.as_ref().unwrap().content, "A");
    // A character the font lacks is refused by name; nothing changes.
    ed.text_edit.as_mut().unwrap().content = "A\u{65E5}".into();
    ed.apply_text(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some("Text not changed: Liberation Sans has no character '\u{65E5}' (U+65E5)")
    );
    assert_eq!(texts(&m, f)[0].2, "A");
    // Selecting something else closes the edit.
    ed.selection.clear();
    ed.sync_text_edit(&m);
    assert!(ed.text_edit.is_none());
}

#[test]
fn a_mirror_image_and_pattern_copies_take_the_new_text() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(0.0, -40.0), (0.0, 40.0)],
    );
    let t = place(
        &mut m,
        &mut ed,
        &map,
        "R",
        Font::Sans,
        (10.0, 0.0),
        (20.0, 10.0),
    );
    // Mirror the text (its box comes with it) about the line.
    ed.selection = vec![t];
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 0.0, 30.0);
    assert_eq!(status(&m, f), SolveStatus::Solved, "{:?}", ed.message);
    ed.escape();
    ed.set_tool(Tool::Select);
    let all = texts(&m, f);
    assert_eq!(all.len(), 2);
    let image = all.iter().find(|x| x.0 != t).unwrap().clone();
    assert!(image.5, "a mirror image of text is mirrored text");
    assert_eq!(construction(&m, f), 8, "the box is mirrored with the text");
    // Two copies of the text along x.
    ed.selection = vec![t];
    ed.set_tool(Tool::LinearPattern);
    click(&mut m, &mut ed, &map, 10.0, 0.0);
    click(&mut m, &mut ed, &map, 10.0, 15.0);
    ed.escape();
    ed.set_tool(Tool::Select);
    assert_eq!(texts(&m, f).len(), 4);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let before = dof(&m, f);
    // Changing the source changes every copy, mirrored ones still mirrored.
    ed.selection = vec![t];
    ed.sync_text_edit(&m);
    ed.text_edit.as_mut().unwrap().content = "Rx".into();
    ed.apply_text(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Text changed"));
    let all = texts(&m, f);
    assert!(all.iter().all(|x| x.2 == "Rx"));
    assert_eq!(all.iter().filter(|x| x.5).count(), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before);
    // The image's box is the source's box mirrored.
    let src = corners_mm(&m, f, all.iter().find(|x| x.0 == t).unwrap().1);
    let img = corners_mm(&m, f, all.iter().find(|x| x.0 == image.0).unwrap().1);
    for i in 0..4 {
        assert!(
            (img[i] - DVec2::new(-src[i].x, src[i].y)).length() < NEAR,
            "{i}: {img:?} {src:?}"
        );
    }
}

#[test]
fn deleting_a_text_takes_its_box_but_not_the_other_way_round() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let t = place(
        &mut m,
        &mut ed,
        &map,
        "Hi",
        Font::Sans,
        (0.0, 0.0),
        (10.0, 10.0),
    );
    // A box line alone goes alone.
    let line = text::box_lines(&m.sketch_def(f).unwrap().sketch, t)[1];
    ed.selection = vec![line];
    ed.delete_selection(&mut m);
    assert_eq!(texts(&m, f).len(), 1);
    assert_eq!(construction(&m, f), 3);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // The text takes its corners and the rest of its box.
    ed.selection = vec![t];
    ed.delete_selection(&mut m);
    assert!(m.sketch_def(f).unwrap().sketch.entities.is_empty());
    assert!(m.sketch_def(f).unwrap().sketch.constraints.is_empty());
}

/// The middle of a vertical straight edge of a letter of `t`, at least
/// `clear` mm from its box's lines (millimetres).
fn on_letter(m: &Model, f: FeatureId, t: SketchId, c: [DVec2; 4], clear: f64) -> DVec2 {
    let (x0, x1, y0, y1) = (c[0].x, c[1].x, c[0].y, c[3].y);
    curves(m, f, false)
        .into_iter()
        .filter(|(id, _)| *id == t)
        .filter_map(|(_, curve)| match curve {
            ResolvedCurve::Line { start, end } if (start.x - end.x).abs() < 1e-12 => {
                Some((start + end) * 0.5 / MM)
            }
            _ => None,
        })
        .find(|p| p.x - x0 > clear && x1 - p.x > clear && p.y - y0 > clear && y1 - p.y > clear)
        .expect("a letter point clear of the box")
}

#[test]
fn letters_take_no_snaps_constraints_cuts_or_fillets() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let t = place(
        &mut m,
        &mut ed,
        &map,
        "H",
        Font::Sans,
        (20.0, 20.0),
        (35.0, 50.0),
    );
    let (_, c, ..) = texts(&m, f).remove(0);
    let k = corners_mm(&m, f, c);
    let p = on_letter(&m, f, t, k, 3.0);
    // No snapping onto a letter.
    ed.set_tool(Tool::Line);
    hover(&m, &mut ed, &map, p.x, p.y);
    assert_eq!(ed.cursor.map(|s| s.target), Some(SnapTarget::Free));
    // No constraint or dimension on the text itself.
    let sketch = m.sketch_def(f).unwrap().sketch.clone();
    assert_eq!(
        constraints_for(ConstraintTool::Coincident, &[c[0], t], &sketch),
        Err(CONSTRAIN_THE_BOX.to_string())
    );
    assert_eq!(
        dimension_for(&[t], &sketch),
        Err(CONSTRAIN_THE_BOX.to_string())
    );
    // Trim, extend and split leave letters alone.
    for tool in [Tool::Trim, Tool::Extend, Tool::Split] {
        ed.set_tool(tool);
        hover(&m, &mut ed, &map, p.x, p.y);
        assert!(
            matches!(&ed.curve_edit, Some(CurveEdit { target, plan: Err(e) })
                if *target == t && e == NOT_EDITABLE),
            "{tool:?}: {:?}",
            ed.curve_edit
        );
    }
    // A fillet at a corner of its box is refused.
    ed.set_tool(Tool::SketchFillet);
    hover(&m, &mut ed, &map, k[2].x, k[2].y);
    assert!(
        matches!(&ed.curve_edit, Some(CurveEdit { plan: Err(e), .. }) if e.contains("square corners")),
        "{:?}",
        ed.curve_edit
    );
    ed.set_tool(Tool::Select);
    // A line across the letter is cut by it, but its new end is free.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(10.0, p.y), (60.0, p.y)],
    );
    let before = curves(&m, f, false).len();
    ed.set_tool(Tool::Trim);
    click(&mut m, &mut ed, &map, 12.0, p.y);
    assert_eq!(ed.message.as_deref(), Some("Trimmed"));
    assert_eq!(curves(&m, f, false).len(), before);
    assert_eq!(
        constraints(
            &m,
            f,
            |c| matches!(c, Constraint::PointOnCurve { curve, .. } if *curve == t)
        ),
        0
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn text_a_font_cannot_draw_is_refused_before_placing() {
    let (mut m, mut ed, map) = setup();
    ed.set_tool(Tool::Text);
    ed.text_content = "Gear \u{2603}".into();
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert!(ed.clicks.is_empty());
    assert_eq!(
        ed.message.as_deref(),
        Some("Text not placed: Liberation Sans has no character '\u{2603}' (U+2603)")
    );
    ed.text_content = String::new();
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert!(ed.clicks.is_empty());
    // A level pointer sets no height.
    ed.text_content = "ok".into();
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    click(&mut m, &mut ed, &map, 30.0, 0.0);
    assert_eq!(ed.clicks.len(), 1);
    assert_eq!(
        ed.message.as_deref(),
        Some("Move the pointer above or below the first click to set the text's height")
    );
}
