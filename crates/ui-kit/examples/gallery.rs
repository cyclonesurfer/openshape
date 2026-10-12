//! Renders the UI kit offscreen and writes four review images:
//!
//! - `icons.png`: every icon at 16, 20, 24 and 32 px on the panel colour,
//!   with its name (and `icons@2x.png`, the same at twice the scale);
//! - `view_cube.png`: the view cube at several orientations, one with a
//!   hovered cell;
//! - `widgets.png`: themed buttons, icon buttons in every state, a list with
//!   a selected and a hovered row, a text field, a tooltip-like panel and
//!   the colour tokens;
//! - `palette.png`: the top bar's command search field and, under it, the
//!   command palette open with the query "ex" (one match disabled).
//!
//! ```sh
//! cargo run -p openshape-ui-kit --example gallery -- <out-dir>
//! ```

use egui::emath::GuiRounding as _;
use egui::{Color32, Pos2, Rect, RichText, Sense, Vec2};
use openshape_render::{write_png, Camera};
use openshape_ui_kit::view_cube::{self, CubeEdge, CubeTarget, CENTER_HALF};
use openshape_ui_kit::{
    icon_button, paint_icon, search_field, CommandPalette, Icon, PaletteEntry, PaletteOutcome,
    Theme, UiSnapshot,
};
use std::path::{Path, PathBuf};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out: PathBuf = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("ui-kit-gallery"));
    std::fs::create_dir_all(&out)?;
    let theme = Theme::light();
    icons(&out, &theme, 1.0, "icons.png")?;
    icons(&out, &theme, 2.0, "icons@2x.png")?;
    cubes(&out, &theme)?;
    widgets(&out, &theme)?;
    palette(&out, &theme)?;
    println!("wrote gallery images to {}", out.display());
    Ok(())
}

/// Runs two passes and saves the second: the first lays the UI out (egui
/// hit-tests against the previous pass's widget rectangles), then the
/// pointer moves to `pointer()` (evaluated after the first pass, so it can
/// use rectangles recorded there) for the second.
fn shoot(
    snap: &mut UiSnapshot,
    path: &Path,
    pointer: impl Fn() -> Option<Pos2>,
    mut run: impl FnMut(&mut egui::Ui),
) -> Result<(), Box<dyn std::error::Error>> {
    snap.frame_ui(egui::RawInput::default(), &mut run)?;
    let mut input = egui::RawInput::default();
    if let Some(p) = pointer() {
        input.events.push(egui::Event::PointerMoved(p));
    }
    let image = snap.frame_ui(input, &mut run)?;
    write_png(path, &image)?;
    Ok(())
}

const SIZES: [f32; 4] = [16.0, 20.0, 24.0, 32.0];

fn icons(
    out: &Path,
    theme: &Theme,
    ppp: f32,
    file: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let columns = 3;
    let rows = Icon::ALL.len().div_ceil(columns);
    let (cell_w, cell_h, top) = (300.0, 40.0, 30.0);
    let size = [
        (columns as f32 * cell_w * ppp) as u32,
        ((top + rows as f32 * cell_h + 10.0) * ppp) as u32,
    ];
    let mut snap = UiSnapshot::new(size, ppp)?;
    theme.apply(snap.context());
    shoot(
        &mut snap,
        &out.join(file),
        || None,
        |ui| {
            let painter = ui.painter();
            painter.rect_filled(ui.max_rect(), 0.0, theme.panel);
            painter.text(
                Pos2::new(10.0, 8.0),
                egui::Align2::LEFT_TOP,
                "OpenShape icons at 16, 20, 24 and 32 px",
                egui::FontId::proportional(theme.heading_size),
                theme.text,
            );
            for (i, &icon) in Icon::ALL.iter().enumerate() {
                let (col, row) = (i / rows, i % rows);
                let origin = Pos2::new(col as f32 * cell_w + 10.0, top + row as f32 * cell_h);
                if row % 2 == 1 {
                    painter.rect_filled(
                        Rect::from_min_size(
                            origin - Vec2::new(10.0, 0.0),
                            Vec2::new(cell_w, cell_h),
                        ),
                        0.0,
                        theme.panel_header,
                    );
                }
                painter.text(
                    origin + Vec2::new(0.0, cell_h * 0.5),
                    egui::Align2::LEFT_CENTER,
                    icon.name(),
                    egui::FontId::proportional(theme.small_size),
                    theme.text_muted,
                );
                let mut x = origin.x + 130.0;
                for s in SIZES {
                    let rect = Rect::from_min_size(
                        Pos2::new(x, origin.y + (cell_h - s) * 0.5),
                        Vec2::splat(s),
                    );
                    // Snap to whole physical pixels like a laid-out widget.
                    let rect = rect.round_to_pixels(ppp);
                    paint_icon(painter, rect, icon, theme.text);
                    x += s + 10.0;
                }
            }
        },
    )
}

fn gradient(painter: &egui::Painter, rect: Rect, top: Color32, bottom: Color32) {
    let mut mesh = egui::Mesh::default();
    mesh.colored_vertex(rect.left_top(), top);
    mesh.colored_vertex(rect.right_top(), top);
    mesh.colored_vertex(rect.right_bottom(), bottom);
    mesh.colored_vertex(rect.left_bottom(), bottom);
    mesh.add_triangle(0, 1, 2);
    mesh.add_triangle(0, 2, 3);
    painter.add(egui::Shape::mesh(mesh));
}

fn cubes(out: &Path, theme: &Theme) -> Result<(), Box<dyn std::error::Error>> {
    let at = |yaw: f64, pitch: f64| Camera {
        yaw: yaw.to_radians(),
        pitch: pitch.to_radians(),
        ..Camera::default()
    };
    let (oblique, below) = (at(30.0, 20.0), at(-125.0, -25.0));
    let views: [(&str, Camera); 5] = [
        ("Isometric (hover: top-right edge)", Camera::isometric()),
        ("Front", Camera::front()),
        ("Top", Camera::top()),
        ("Yaw 30, pitch 20", oblique),
        ("Yaw -125, pitch -25", below),
    ];
    let (cell, side, ppp) = (170.0, 120.0, 2.0);
    let size = [
        (cell * views.len() as f32 * ppp) as u32,
        (200.0 * ppp) as u32,
    ];
    let rect_of = |i: usize| {
        Rect::from_min_size(
            Pos2::new(i as f32 * cell + (cell - side) * 0.5, 20.0),
            Vec2::splat(side),
        )
    };
    // Hover the edge cell between top and right on the first cube.
    let hover_target = CubeTarget::Edge(CubeEdge::TopRight);
    let [x, y, z] = hover_target.direction().map(|c| c as f64);
    // Middle of that edge's cell on the top face (edge bands span
    // CENTER_HALF..1 in face coordinates).
    let band_mid = (1.0 + CENTER_HALF) * 0.5;
    let hover = view_cube::project(&views[0].1, rect_of(0), [x * band_mid, y * band_mid, z]);

    let mut snap = UiSnapshot::new(size, ppp)?;
    theme.apply(snap.context());
    shoot(
        &mut snap,
        &out.join("view_cube.png"),
        || Some(hover),
        |ui| {
            let full = ui.max_rect();
            gradient(
                ui.painter(),
                full,
                theme.viewport_top,
                theme.viewport_bottom,
            );
            for (i, (caption, camera)) in views.iter().enumerate() {
                let rect = rect_of(i);
                // Nothing clicks in a snapshot; the returned target is unused.
                let _ = view_cube::show(ui, rect, camera, theme);
                ui.painter().text(
                    Pos2::new(rect.center().x, rect.bottom() + 22.0),
                    egui::Align2::CENTER_TOP,
                    *caption,
                    egui::FontId::proportional(theme.small_size),
                    theme.text_muted,
                );
            }
        },
    )?;
    // Sanity: the pointer really is over the target we meant to hover.
    let picked = view_cube::pick(&views[0].1, rect_of(0), hover);
    assert_eq!(picked, Some(hover_target), "hover point picks {picked:?}");
    Ok(())
}

/// A list row: icon and label, with selected and hovered backgrounds.
fn row(ui: &mut egui::Ui, theme: &Theme, icon: Icon, label: &str, selected: bool) {
    let (rect, response) = ui.allocate_exact_size(
        Vec2::new(ui.available_width(), theme.interact_height),
        Sense::click(),
    );
    let fill = if selected {
        Some(theme.accent_soft)
    } else if response.hovered() {
        Some(theme.hover_soft)
    } else {
        None
    };
    if let Some(fill) = fill {
        ui.painter().rect_filled(rect, 0.0, fill);
    }
    let ink = if selected { theme.accent } else { theme.text };
    let icon_rect = Rect::from_min_size(
        Pos2::new(rect.left() + 6.0, rect.center().y - 8.0),
        Vec2::splat(16.0),
    );
    paint_icon(ui.painter(), icon_rect, icon, ink);
    ui.painter().text(
        Pos2::new(rect.left() + 28.0, rect.center().y),
        egui::Align2::LEFT_CENTER,
        label,
        egui::FontId::proportional(theme.body_size),
        theme.text,
    );
}

fn widgets(out: &Path, theme: &Theme) -> Result<(), Box<dyn std::error::Error>> {
    let ppp = 2.0;
    let size = [(720.0 * ppp) as u32, (470.0 * ppp) as u32];
    // The pointer rests on the second icon button (its rectangle is
    // recorded in the first pass).
    let hovered_button = std::cell::Cell::new(None::<Rect>);
    let mut snap = UiSnapshot::new(size, ppp)?;
    theme.apply(snap.context());
    let mut text = String::from("25 mm");
    let pointer = || hovered_button.get().map(|r| r.center());
    shoot(&mut snap, &out.join("widgets.png"), pointer, |ui| {
        let full = ui.max_rect();
        ui.painter().rect_filled(full, 0.0, theme.chrome);
        let panel = Rect::from_min_size(Pos2::new(8.0, 8.0), Vec2::new(340.0, 454.0));
        ui.painter().rect(
            panel,
            theme.rounding,
            theme.panel,
            egui::Stroke::new(1.0, theme.border),
            egui::StrokeKind::Inside,
        );
        let mut left =
            ui.new_child(egui::UiBuilder::new().max_rect(panel.shrink(theme.panel_margin)));
        left.label(RichText::new("Icon buttons").size(theme.heading_size));
        left.horizontal(|ui| {
            icon_button(ui, theme, Icon::Sketch, "Sketch (at rest)", false, true);
            let r = icon_button(ui, theme, Icon::Extrude, "Extrude (hovered)", false, true);
            hovered_button.set(Some(r.rect));
            icon_button(ui, theme, Icon::Line, "Line (selected)", true, true);
            icon_button(ui, theme, Icon::Revolve, "Revolve (disabled)", false, false);
            icon_button(
                ui,
                theme,
                Icon::Circle,
                "Circle (selected, disabled)",
                true,
                false,
            );
        });
        left.label(
            RichText::new("rest, hovered, selected, disabled, selected + disabled")
                .size(theme.small_size)
                .color(theme.text_muted),
        );
        left.add_space(8.0);
        left.label(RichText::new("Buttons").size(theme.heading_size));
        left.horizontal(|ui| {
            let _ = ui.button("Cancel");
            let _ = ui.add(
                egui::Button::new(RichText::new("Accept").color(theme.on_accent))
                    .fill(theme.accent),
            );
            ui.add_enabled(false, egui::Button::new("Disabled"));
            let mut on = true;
            ui.checkbox(&mut on, "Merge");
        });
        left.add_space(8.0);
        left.label(RichText::new("Text field").size(theme.heading_size));
        left.horizontal(|ui| {
            ui.label("Depth");
            ui.add(egui::TextEdit::singleline(&mut text).desired_width(120.0));
        });
        left.add_space(8.0);
        left.label(RichText::new("Feature list").size(theme.heading_size));
        let list = [
            (Icon::Origin, "Origin", false),
            (Icon::Plane, "Top plane", false),
            (Icon::Sketch, "Sketch 1", true),
            (Icon::Extrude, "Extrude 1", false),
            (Icon::Revolve, "Revolve 1", false),
        ];
        left.spacing_mut().item_spacing.y = 0.0;
        for (icon, label, selected) in list {
            row(&mut left, theme, icon, label, selected);
        }

        // Right column: a hovered list row (the pointer can only be in one
        // place, so this row is shown hovered by drawing its state), a
        // tooltip-like panel and the colour tokens.
        let right = Rect::from_min_size(Pos2::new(360.0, 8.0), Vec2::new(352.0, 454.0));
        ui.painter().rect(
            right,
            theme.rounding,
            theme.panel,
            egui::Stroke::new(1.0, theme.border),
            egui::StrokeKind::Inside,
        );
        let mut r = ui.new_child(egui::UiBuilder::new().max_rect(right.shrink(theme.panel_margin)));
        r.label(RichText::new("Rows: selected and hovered").size(theme.heading_size));
        let (rect, _) = r.allocate_exact_size(
            Vec2::new(r.available_width(), theme.interact_height),
            Sense::hover(),
        );
        r.painter().rect_filled(rect, 0.0, theme.accent_soft);
        paint_icon(
            r.painter(),
            Rect::from_min_size(
                Pos2::new(rect.left() + 6.0, rect.center().y - 8.0),
                Vec2::splat(16.0),
            ),
            Icon::Part,
            theme.accent,
        );
        r.painter().text(
            Pos2::new(rect.left() + 28.0, rect.center().y),
            egui::Align2::LEFT_CENTER,
            "Part 1 (selected)",
            egui::FontId::proportional(theme.body_size),
            theme.text,
        );
        let (rect, _) = r.allocate_exact_size(
            Vec2::new(r.available_width(), theme.interact_height),
            Sense::hover(),
        );
        r.painter().rect_filled(rect, 0.0, theme.hover_soft);
        paint_icon(
            r.painter(),
            Rect::from_min_size(
                Pos2::new(rect.left() + 6.0, rect.center().y - 8.0),
                Vec2::splat(16.0),
            ),
            Icon::Part,
            theme.text,
        );
        r.painter().text(
            Pos2::new(rect.left() + 28.0, rect.center().y),
            egui::Align2::LEFT_CENTER,
            "Part 2 (hovered)",
            egui::FontId::proportional(theme.body_size),
            theme.text,
        );
        r.add_space(10.0);
        r.label(RichText::new("Tooltip-like panel").size(theme.heading_size));
        egui::Frame::popup(r.style()).show(&mut r, |ui| {
            ui.label(RichText::new("Extrude").strong());
            ui.label(
                RichText::new("Adds or removes material along a direction. Shortcut: E")
                    .size(theme.small_size)
                    .color(theme.text_muted),
            );
        });
        r.add_space(10.0);
        r.label(RichText::new("Tokens").size(theme.heading_size));
        let tokens = [
            ("chrome", theme.chrome),
            ("panel", theme.panel),
            ("panel_header", theme.panel_header),
            ("border", theme.border),
            ("border_strong", theme.border_strong),
            ("text", theme.text),
            ("text_muted", theme.text_muted),
            ("text_disabled", theme.text_disabled),
            ("accent", theme.accent),
            ("accent_hover", theme.accent_hover),
            ("accent_soft", theme.accent_soft),
            ("hover_soft", theme.hover_soft),
            ("selection", theme.selection),
            ("preselect", theme.preselect),
            ("success", theme.success),
            ("warning", theme.warning),
            ("danger", theme.danger),
            ("viewport_bottom", theme.viewport_bottom),
        ];
        egui::Grid::new("tokens")
            .num_columns(4)
            .spacing([8.0, 2.0])
            .show(&mut r, |ui| {
                for (i, (name, color)) in tokens.iter().enumerate() {
                    let (rect, _) = ui.allocate_exact_size(Vec2::new(22.0, 14.0), Sense::hover());
                    ui.painter().rect(
                        rect,
                        2.0,
                        *color,
                        egui::Stroke::new(1.0, theme.border),
                        egui::StrokeKind::Inside,
                    );
                    ui.label(RichText::new(*name).size(theme.small_size));
                    if i % 2 == 1 {
                        ui.end_row();
                    }
                }
            });
    })
}

/// The commands of a part studio, as the app would list them.
fn commands() -> Vec<PaletteEntry<'static>> {
    let entry = |label, detail, icon, shortcut, keywords, enabled| PaletteEntry {
        label,
        detail,
        icon: Some(icon),
        shortcut,
        keywords,
        enabled,
    };
    vec![
        entry(
            "Sketch",
            "Feature",
            Icon::Sketch,
            Some("Shift+S"),
            &["draw", "profile"],
            true,
        ),
        entry(
            "Extrude",
            "Feature",
            Icon::Extrude,
            Some("Shift+E"),
            &["pad", "boss", "pocket", "cut"],
            true,
        ),
        entry(
            "Revolve",
            "Feature",
            Icon::Revolve,
            Some("Shift+W"),
            &["lathe", "turn"],
            true,
        ),
        entry("Box", "Feature", Icon::Box, None, &["cube", "block"], true),
        entry("Cylinder", "Feature", Icon::Cylinder, None, &["rod"], true),
        entry(
            "Export STL",
            "File",
            Icon::Export,
            None,
            &["save", "mesh", "print"],
            true,
        ),
        entry("Undo", "Edit", Icon::Undo, Some("Ctrl+Z"), &[], true),
        entry("Redo", "Edit", Icon::Redo, Some("Ctrl+Y"), &[], true),
        entry(
            "Fit view",
            "View",
            Icon::ZoomFit,
            Some("F"),
            &["zoom", "extents"],
            true,
        ),
        // Only available while a sketch is being edited.
        entry(
            "Exit sketch",
            "Sketch",
            Icon::FinishSketch,
            None,
            &["finish", "done"],
            false,
        ),
    ]
}

fn palette(out: &Path, theme: &Theme) -> Result<(), Box<dyn std::error::Error>> {
    let (ppp, size) = (2.0, Vec2::new(720.0, 340.0));
    let mut snap = UiSnapshot::new([(size.x * ppp) as u32, (size.y * ppp) as u32], ppp)?;
    theme.apply(snap.context());
    let entries = commands();
    let mut palette = CommandPalette::new();
    palette.query = "ex".to_owned();
    let bar = Rect::from_min_size(Pos2::ZERO, Vec2::new(size.x, 40.0));
    let field = Rect::from_center_size(bar.center(), Vec2::new(320.0, theme.interact_height));
    let anchor = field.center_bottom() + Vec2::new(0.0, theme.item_spacing.y + 2.0);
    let mut run = |ui: &mut egui::Ui| {
        let full = ui.max_rect();
        gradient(
            ui.painter(),
            full,
            theme.viewport_top,
            theme.viewport_bottom,
        );
        ui.painter().rect_filled(bar, 0.0, theme.panel);
        ui.painter().rect_filled(
            Rect::from_min_size(bar.left_bottom(), Vec2::new(bar.width(), 1.0)),
            0.0,
            theme.border,
        );
        ui.painter().text(
            Pos2::new(bar.left() + 12.0, bar.center().y),
            egui::Align2::LEFT_CENTER,
            "Part Studio 1",
            egui::FontId::proportional(theme.body_size),
            theme.text,
        );
        let mut slot = ui.new_child(egui::UiBuilder::new().max_rect(field));
        let _ = search_field(&mut slot, theme, "Search commands", "Ctrl+K", field.width());
        let outcome = palette.show(ui.ctx(), theme, anchor, 480.0, &entries);
        assert_eq!(outcome, PaletteOutcome::Open);
    };
    // The palette's area is sized in its first frame and drawn from the
    // second; a few more let the text cursor and fades settle.
    for _ in 0..5 {
        snap.frame_ui(egui::RawInput::default(), &mut run)?;
    }
    let image = snap.frame_ui(egui::RawInput::default(), &mut run)?;
    write_png(out.join("palette.png"), &image)?;
    Ok(())
}
