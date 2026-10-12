//! The icon set: every icon draws inside its square and nowhere else, no two
//! icons look the same, and the toolbar button shows its states. Pixel tests
//! run on the GPU (lavapipe in CI).

mod common;

use common::{max_channel_diff, maybe_dump, rgba, snapshot};
use egui::{Color32, LayerId, Pos2, Rect, Vec2};
use openshape_render::RgbaImage;
use openshape_ui_kit::{icon_button, paint_icon, Icon, Theme};
use std::collections::HashSet;

const COLUMNS: usize = 10;

#[test]
fn the_set_is_complete_and_named() {
    // 19 app, 7 view, 11 tree and tab, 39 sketch icons.
    assert_eq!(Icon::ALL.len(), 77);
    let unique: HashSet<Icon> = Icon::ALL.iter().copied().collect();
    assert_eq!(unique.len(), Icon::ALL.len(), "duplicates in Icon::ALL");
    let names: HashSet<&str> = Icon::ALL.iter().map(|i| i.name()).collect();
    assert_eq!(names.len(), Icon::ALL.len(), "duplicate names");
    assert!(Icon::ALL.iter().all(|i| !i.name().trim().is_empty()));
    assert_eq!(Icon::CenterArc.name(), "Center arc");
}

/// Renders every icon at `size` px in its own `size + 2 * margin` cell
/// (pixel-aligned) on white, in dark ink.
fn sheet(size: u32, margin: u32) -> (RgbaImage, impl Fn(usize) -> (u32, u32)) {
    let cell = size + 2 * margin;
    let rows = Icon::ALL.len().div_ceil(COLUMNS) as u32;
    let mut snap = snapshot([cell * COLUMNS as u32, cell * rows], 1.0);
    let origin = move |i: usize| {
        (
            (i % COLUMNS) as u32 * cell + margin,
            (i / COLUMNS) as u32 * cell + margin,
        )
    };
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            let painter = ctx.layer_painter(LayerId::background());
            painter.rect_filled(ctx.content_rect(), 0.0, Color32::WHITE);
            for (i, &icon) in Icon::ALL.iter().enumerate() {
                let (x, y) = origin(i);
                let rect =
                    Rect::from_min_size(Pos2::new(x as f32, y as f32), Vec2::splat(size as f32));
                paint_icon(&painter, rect, icon, Color32::from_rgb(20, 20, 20));
            }
        })
        .unwrap();
    (image, origin)
}

/// Any visible trace of ink, antialiasing included.
fn is_ink(p: [u8; 4]) -> bool {
    max_channel_diff(p, [255, 255, 255, 255]) > 24
}

/// A clearly inked pixel (at least about a third covered).
fn is_solid_ink(p: [u8; 4]) -> bool {
    max_channel_diff(p, [255, 255, 255, 255]) > 80
}

#[test]
fn every_icon_draws_inside_its_square_and_nowhere_else() {
    for size in [16, 20, 24, 32] {
        let margin = 6;
        let (image, origin) = sheet(size, margin);
        maybe_dump(&format!("icons_{size}px"), &image);
        let cell = size + 2 * margin;
        for (i, icon) in Icon::ALL.iter().enumerate() {
            let (x0, y0) = origin(i);
            let (cx, cy) = (x0 - margin, y0 - margin);
            let mut ink_inside = 0;
            let mut outside = Vec::new();
            for y in cy..cy + cell {
                for x in cx..cx + cell {
                    let p = image.pixel(x, y);
                    let inside = (x0..x0 + size).contains(&x) && (y0..y0 + size).contains(&y);
                    if inside {
                        ink_inside += is_ink(p) as u32;
                    } else if p != [255, 255, 255, 255] {
                        outside.push((x - cx, y - cy, p));
                    }
                }
            }
            assert!(
                ink_inside >= size,
                "{icon:?} at {size}px has only {ink_inside} ink pixels"
            );
            assert!(
                outside.is_empty(),
                "{icon:?} at {size}px paints outside its square: {:?}",
                &outside[..outside.len().min(5)]
            );
        }
    }
}

#[test]
fn no_two_icons_look_alike_at_24px() {
    let (image, origin) = sheet(24, 2);
    let blocks: Vec<Vec<bool>> = (0..Icon::ALL.len())
        .map(|i| {
            let (x0, y0) = origin(i);
            (0..24 * 24)
                .map(|k| is_solid_ink(image.pixel(x0 + k % 24, y0 + k / 24)))
                .collect()
        })
        .collect();
    let mut closest = (usize::MAX, Icon::New, Icon::New);
    for i in 0..blocks.len() {
        for j in i + 1..blocks.len() {
            let differing = blocks[i]
                .iter()
                .zip(&blocks[j])
                .filter(|(a, b)| a != b)
                .count();
            if differing < closest.0 {
                closest = (differing, Icon::ALL[i], Icon::ALL[j]);
            }
        }
    }
    // Identical images would differ in 0 pixels; demand a visible margin.
    assert!(
        closest.0 >= 20,
        "{:?} and {:?} differ in only {} pixels",
        closest.1,
        closest.2,
        closest.0
    );
}

#[test]
fn icons_centre_in_non_square_rectangles_and_take_the_colour() {
    let mut snap = snapshot([48, 24], 1.0);
    let accent = Theme::light().accent;
    let image = snap
        .frame(egui::RawInput::default(), |ctx| {
            let painter = ctx.layer_painter(LayerId::background());
            painter.rect_filled(ctx.content_rect(), 0.0, Color32::WHITE);
            paint_icon(
                &painter,
                Rect::from_min_size(Pos2::ZERO, Vec2::new(48.0, 24.0)),
                Icon::Point,
                accent,
            );
        })
        .unwrap();
    for y in 0..24 {
        for x in (0..12).chain(36..48) {
            assert_eq!(image.pixel(x, y), [255, 255, 255, 255], "ink at ({x}, {y})");
        }
    }
    // The point's dot is solid at the centre of the 24-unit square.
    assert!(max_channel_diff(image.pixel(24, 12), rgba(accent)) <= 1);
}

#[test]
fn icon_buttons_show_rest_selected_disabled_and_hovered_states() {
    let theme = Theme::light();
    let mut snap = snapshot([200, 40], 1.0);
    theme.apply(snap.context());
    let rects = std::cell::RefCell::new(Vec::new());
    let mut run = |ui: &mut egui::Ui| {
        ui.painter().rect_filled(ui.max_rect(), 0.0, theme.panel);
        ui.horizontal(|ui| {
            let mut r = rects.borrow_mut();
            r.clear();
            // Point icons have a solid dot in the middle that shows the ink.
            r.push(icon_button(ui, &theme, Icon::Point, "rest", false, true).rect);
            r.push(icon_button(ui, &theme, Icon::Point, "selected", true, true).rect);
            r.push(icon_button(ui, &theme, Icon::Point, "disabled", false, false).rect);
            r.push(icon_button(ui, &theme, Icon::Point, "hovered", false, true).rect);
        });
    };
    snap.frame_ui(egui::RawInput::default(), &mut run).unwrap();
    let hover = rects.borrow()[3].center() + Vec2::new(6.0, 6.0);
    let mut input = egui::RawInput::default();
    input.events.push(egui::Event::PointerMoved(hover));
    let image = snap.frame_ui(input, &mut run).unwrap();
    maybe_dump("icon_button_states", &image);

    let rects = rects.borrow();
    for r in rects.iter() {
        assert_eq!(r.size(), Vec2::splat(theme.toolbar_button));
    }
    let background = |i: usize| {
        let p = rects[i].left_top() + Vec2::splat(2.0);
        image.pixel(p.x as u32, p.y as u32)
    };
    let ink = |i: usize| {
        let p = rects[i].center();
        image.pixel(p.x as u32, p.y as u32)
    };
    let expect = |got: [u8; 4], want: Color32, what: &str| {
        assert!(
            max_channel_diff(got, rgba(want)) <= 1,
            "{what}: {got:?} vs {want:?}"
        );
    };
    expect(background(0), theme.panel, "rest background");
    expect(background(1), theme.accent_soft, "selected background");
    expect(background(2), theme.panel, "disabled background");
    expect(background(3), theme.hover_soft, "hovered background");
    expect(ink(0), theme.text, "rest ink");
    expect(ink(1), theme.accent, "selected ink");
    expect(ink(2), theme.text_disabled, "disabled ink");
    expect(ink(3), theme.text, "hovered ink");
}

/// Clicks the first (and only) button of `ui` built by `add`, over three
/// passes; returns whether the button reported a click.
fn click(enabled: bool) -> bool {
    let ctx = egui::Context::default();
    let theme = Theme::light();
    let mut clicked = false;
    let mut rect = Rect::NOTHING;
    let pos = |r: Rect| r.center();
    let events = |r: Rect| {
        vec![
            vec![],
            vec![
                egui::Event::PointerMoved(pos(r)),
                egui::Event::PointerButton {
                    pos: pos(r),
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: Default::default(),
                },
            ],
            vec![egui::Event::PointerButton {
                pos: pos(r),
                button: egui::PointerButton::Primary,
                pressed: false,
                modifiers: Default::default(),
            }],
        ]
    };
    for pass in 0..3 {
        let input = egui::RawInput {
            events: events(rect)[pass].clone(),
            ..Default::default()
        };
        let mut out = ctx.run_ui(input, |ui| {
            let r = icon_button(ui, &theme, Icon::Check, "Accept", false, enabled);
            rect = r.rect;
            clicked |= r.clicked();
        });
        out.textures_delta.clear();
    }
    clicked
}

#[test]
fn enabled_buttons_click_and_disabled_ones_do_not() {
    assert!(click(true));
    assert!(!click(false));
}
