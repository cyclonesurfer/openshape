//! The theme installs light visuals with OpenDatum's tokens, its spacing and
//! text sizes, and fonts that cover the constraint glyphs. No GPU needed.

use egui::{Color32, FontFamily, TextStyle};
use opendatum_ui_kit::{font_definitions, Theme};

/// Runs one pass so pending fonts and styles take effect, dropping its
/// texture uploads (no renderer here).
fn settle(ctx: &egui::Context) {
    let mut out = ctx.run_ui(egui::RawInput::default(), |_| {});
    out.textures_delta.clear();
}

#[test]
fn light_is_the_default_and_uses_the_specified_tokens() {
    let t = Theme::default();
    assert_eq!(t, Theme::light());
    let hex = |c: Color32| format!("#{:02X}{:02X}{:02X}", c.r(), c.g(), c.b());
    let tokens = [
        (t.chrome, "#EEF0F3"),
        (t.panel, "#FFFFFF"),
        (t.panel_header, "#F6F7F9"),
        (t.border, "#D5D9DF"),
        (t.border_strong, "#AEB4BD"),
        (t.text, "#1E2329"),
        (t.text_muted, "#5B6470"),
        (t.text_disabled, "#A3A9B2"),
        (t.accent, "#2B6CD4"),
        (t.accent_hover, "#1F5BBF"),
        (t.accent_soft, "#E1ECFB"),
        (t.on_accent, "#FFFFFF"),
        (t.hover_soft, "#EEF3FA"),
        (t.selection, "#2B6CD4"),
        (t.preselect, "#6FA3EE"),
        (t.success, "#2E9A57"),
        (t.warning, "#C98112"),
        (t.danger, "#CC3D33"),
        (t.viewport_top, "#FFFFFF"),
        (t.viewport_bottom, "#D9E0E8"),
    ];
    for (color, expected) in tokens {
        assert_eq!(hex(color), expected);
        assert!(color.is_opaque());
    }
    assert_eq!(t.rounding, 3.0);
    assert_eq!(t.item_spacing, egui::vec2(6.0, 4.0));
    assert_eq!(t.button_padding, egui::vec2(6.0, 3.0));
    assert_eq!(t.interact_height, 24.0);
    assert_eq!(t.panel_margin, 6.0);
    assert_eq!(
        (t.body_size, t.small_size, t.heading_size),
        (13.0, 11.5, 14.0)
    );
    assert_eq!((t.toolbar_button, t.toolbar_icon), (28.0, 20.0));
}

#[test]
fn apply_installs_light_visuals_with_the_token_colours() {
    let theme = Theme::light();
    let ctx = egui::Context::default();
    // Even if the system asks for dark mode, the app stays light.
    ctx.set_theme(egui::ThemePreference::Dark);
    theme.apply(&ctx);
    settle(&ctx);

    let style = ctx.global_style();
    let v = &style.visuals;
    assert!(!v.dark_mode);
    assert_eq!(v.panel_fill, theme.panel);
    assert_eq!(v.window_fill, theme.panel);
    assert_eq!(v.window_stroke, egui::Stroke::new(1.0, theme.border));
    assert_eq!(v.extreme_bg_color, theme.panel);
    assert_eq!(v.faint_bg_color, theme.panel_header);
    assert_eq!(v.hyperlink_color, theme.accent);
    assert_eq!(v.warn_fg_color, theme.warning);
    assert_eq!(v.error_fg_color, theme.danger);
    assert_eq!(v.selection.bg_fill, theme.accent_soft);
    assert_eq!(v.selection.stroke.color, theme.accent);
    assert_eq!(v.weak_text_color, Some(theme.text_muted));
    assert_eq!(v.text_color(), theme.text);
    assert_eq!(v.widgets.noninteractive.bg_stroke.color, theme.border);
    assert_eq!(v.widgets.inactive.weak_bg_fill, theme.panel);
    assert_eq!(v.widgets.inactive.bg_stroke.color, theme.border);
    assert_eq!(v.widgets.hovered.weak_bg_fill, theme.hover_soft);
    assert_eq!(v.widgets.hovered.bg_stroke.color, theme.border_strong);
    assert_eq!(v.widgets.active.weak_bg_fill, theme.accent_soft);
    for w in [
        &v.widgets.noninteractive,
        &v.widgets.inactive,
        &v.widgets.hovered,
        &v.widgets.active,
        &v.widgets.open,
    ] {
        assert_eq!(w.corner_radius, egui::CornerRadius::same(3));
        assert_eq!(w.fg_stroke.color, theme.text);
        assert!(w.bg_stroke.width <= 1.0, "borders are thin");
        assert_eq!(w.expansion, 0.0);
    }
    assert_eq!(v.window_corner_radius, egui::CornerRadius::same(3));

    let s = &style.spacing;
    assert_eq!(s.item_spacing, theme.item_spacing);
    assert_eq!(s.button_padding, theme.button_padding);
    assert_eq!(s.interact_size.y, theme.interact_height);
    assert_eq!(s.window_margin, egui::Margin::same(6));

    let size = |ts: TextStyle| style.text_styles[&ts].size;
    assert_eq!(size(TextStyle::Body), 13.0);
    assert_eq!(size(TextStyle::Button), 13.0);
    assert_eq!(size(TextStyle::Small), 11.5);
    assert_eq!(size(TextStyle::Heading), 14.0);
    assert_eq!(
        style.text_styles[&TextStyle::Body].family,
        FontFamily::Proportional
    );
    assert_eq!(
        style.text_styles[&TextStyle::Monospace].family,
        FontFamily::Monospace
    );

    // The panel frame helper uses the same tokens.
    let frame = theme.panel_frame();
    assert_eq!(frame.fill, theme.panel);
    assert_eq!(frame.inner_margin, egui::Margin::same(6));
}

#[test]
fn font_definitions_append_hack_once_as_the_last_proportional_fallback() {
    let defaults = egui::FontDefinitions::default();
    let fonts = font_definitions();
    let proportional = &fonts.families[&FontFamily::Proportional];
    assert_eq!(proportional.last().map(String::as_str), Some("Hack"));
    assert_eq!(proportional.iter().filter(|n| *n == "Hack").count(), 1);
    // Everything else is egui's defaults, in order.
    assert_eq!(
        &proportional[..proportional.len() - 1],
        &defaults.families[&FontFamily::Proportional][..]
    );
    assert_eq!(
        fonts.families[&FontFamily::Monospace],
        defaults.families[&FontFamily::Monospace]
    );
}

#[test]
fn proportional_text_has_the_constraint_glyphs() {
    let ctx = egui::Context::default();
    Theme::light().apply(&ctx);
    settle(&ctx);
    // `has_glyph` misreports characters in the same face as the replacement
    // glyph, so read the family's character map instead.
    let chars = ctx.fonts_mut(|f| f.fonts.font(&FontFamily::Proportional).characters().clone());
    for c in ['\u{2016}', '\u{27C2}', '\u{00B0}', 'A'] {
        assert!(
            chars.contains_key(&c),
            "proportional fonts lack U+{:04X}",
            c as u32
        );
    }
}

/// WCAG 2.x relative luminance of an sRGB colour.
fn luminance(c: Color32) -> f64 {
    let lin = |v: u8| {
        let v = f64::from(v) / 255.0;
        if v <= 0.04045 {
            v / 12.92
        } else {
            ((v + 0.055) / 1.055).powf(2.4)
        }
    };
    0.2126 * lin(c.r()) + 0.7152 * lin(c.g()) + 0.0722 * lin(c.b())
}

/// WCAG 2.x contrast ratio.
fn contrast(a: Color32, b: Color32) -> f64 {
    let (la, lb) = (luminance(a), luminance(b));
    (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
}

#[test]
fn the_palette_keeps_its_contrast_ratios() {
    // The criteria ADR 0012 chose the palette by.
    let t = Theme::light();
    let at_least = |name: &str, a: Color32, b: Color32, min: f64| {
        let r = contrast(a, b);
        assert!(r >= min, "{name}: {r:.2} < {min}");
    };
    at_least("text on panel", t.text, t.panel, 7.0);
    at_least("text on chrome", t.text, t.chrome, 7.0);
    at_least("text on selected row", t.text, t.accent_soft, 7.0);
    at_least("muted text on panel", t.text_muted, t.panel, 4.5);
    at_least("muted text on chrome", t.text_muted, t.chrome, 4.5);
    at_least("text on accent", t.on_accent, t.accent, 4.5);
    at_least("accent text on panel", t.accent, t.panel, 4.5);
    at_least("error text on panel", t.danger, t.panel, 4.5);
    at_least("warning icon on panel", t.warning, t.panel, 3.0);
    at_least("success icon on panel", t.success, t.panel, 3.0);
    at_least("accent icon on selected row", t.accent, t.accent_soft, 3.0);
    let k = &t.sketch;
    for bg in [t.viewport_top, t.viewport_bottom] {
        at_least("constrained sketch", k.constrained, bg, 7.0);
        at_least("free sketch", k.free, bg, 3.0);
        at_least("conflicting sketch", k.conflict, bg, 3.0);
        at_least("inactive sketch", k.inactive, bg, 3.0);
    }
    at_least("dimension on its label", k.dimension, t.panel, 7.0);
}
