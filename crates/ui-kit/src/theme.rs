//! OpenShape's theme: colour, spacing and type tokens, and their mapping
//! onto egui's [`egui::Style`].
//!
//! The look is compact, light and crisp: white panels on a cool light-grey
//! chrome, thin 1 px borders, small corner radii, subtle shadows and a
//! single blue accent for selection and focus. The palette is OpenShape's
//! own; every token is a public field of [`Theme`] so custom widgets draw
//! with the same values that [`Theme::apply`] gives egui's built-in ones.

use egui::{Color32, CornerRadius, FontFamily, FontId, Margin, Shadow, Stroke, TextStyle, Vec2};

/// `0xRRGGBB` as an opaque sRGB colour.
const fn hex(rgb: u32) -> Color32 {
    Color32::from_rgb((rgb >> 16) as u8, (rgb >> 8) as u8, rgb as u8)
}

/// Colours of sketch geometry, annotations and other graphics-area overlays.
/// Selection and preselection use [`Theme::selection`] and
/// [`Theme::preselect`], drawn with a [`SketchColors::halo`] under them.
#[derive(Clone, Debug, PartialEq)]
pub struct SketchColors {
    /// Fully constrained geometry.
    pub constrained: Color32,
    /// Geometry that can still move (under-constrained).
    pub free: Color32,
    /// Geometry and constraints in conflict.
    pub conflict: Color32,
    /// Redundant constraints.
    pub redundant: Color32,
    /// Sketches that are not being edited.
    pub inactive: Color32,
    /// Dimension text and lines.
    pub dimension: Color32,
    /// Constraint glyphs.
    pub glyph: Color32,
    /// Background of dimension and glyph labels.
    pub label_fill: Color32,
    /// Translucent band drawn under selected or preselected geometry.
    pub halo: Color32,
    /// Fill of picked profile regions.
    pub region: Color32,
    /// Fill of the region under the pointer while picking.
    pub region_hover: Color32,
    /// The sketch plane's x axis.
    pub axis_x: Color32,
    /// The sketch plane's y axis.
    pub axis_y: Color32,
    /// Fill of the default planes in the 3D view.
    pub plane: Color32,
    /// Outline of the default planes.
    pub plane_edge: Color32,
    /// Fill of a selected or hovered default plane.
    pub plane_active: Color32,
    /// Outline of a selected or hovered default plane.
    pub plane_active_edge: Color32,
    /// Plane labels.
    pub plane_label: Color32,
    /// Teeth of curvature combs.
    pub comb: Color32,
    /// The line joining a curvature comb's tooth tips.
    pub comb_envelope: Color32,
}

/// The app's visual tokens. [`Theme::light`] (also [`Default`]) is the only
/// palette today.
#[derive(Clone, Debug, PartialEq)]
pub struct Theme {
    /// App background behind panels (title strip, gaps between panels).
    pub chrome: Color32,
    /// Panel, window and menu background.
    pub panel: Color32,
    /// Panel header strips (section titles, table headers).
    pub panel_header: Color32,
    /// Thin borders and separators.
    pub border: Color32,
    /// Emphasised borders (hovered widgets, the view cube's edges).
    pub border_strong: Color32,
    /// Normal text and icons.
    pub text: Color32,
    /// Secondary text (hints, units, captions).
    pub text_muted: Color32,
    /// Text and icons of disabled controls.
    pub text_disabled: Color32,
    /// Accent: primary buttons, focus rings, selected icons.
    pub accent: Color32,
    /// Accent of a hovered primary button.
    pub accent_hover: Color32,
    /// Background of selected rows and toggled tool buttons.
    pub accent_soft: Color32,
    /// Text and icons drawn on an [`Theme::accent`] background.
    pub on_accent: Color32,
    /// Background of hovered rows and tool buttons.
    pub hover_soft: Color32,
    /// Selected geometry in the 3D view and in sketches.
    pub selection: Color32,
    /// Geometry (or view-cube cell) under the pointer, before a click.
    pub preselect: Color32,
    /// Success state (regenerated feature, satisfied sketch).
    pub success: Color32,
    /// Warning state.
    pub warning: Color32,
    /// Error state (failed feature, over-constrained sketch).
    pub danger: Color32,
    /// Top colour of the graphics-area background gradient.
    pub viewport_top: Color32,
    /// Bottom colour of the graphics-area background gradient.
    pub viewport_bottom: Color32,
    /// Sketch and 3D-view colours.
    pub sketch: SketchColors,

    /// Corner radius of widgets, windows and menus, in points.
    pub rounding: f32,
    /// Gap between neighbouring widgets (x, y), in points.
    pub item_spacing: Vec2,
    /// Padding inside text buttons (x, y), in points.
    pub button_padding: Vec2,
    /// Minimum height of interactive rows and buttons, in points.
    pub interact_height: f32,
    /// Inner margin of panels and windows, in points.
    pub panel_margin: f32,
    /// Body text size, in points.
    pub body_size: f32,
    /// Small text size (captions, status bar), in points.
    pub small_size: f32,
    /// Heading text size, in points.
    pub heading_size: f32,
    /// Side of a square toolbar button, in points.
    pub toolbar_button: f32,
    /// Side of the icon inside a toolbar button, in points.
    pub toolbar_icon: f32,
}

impl Default for Theme {
    fn default() -> Self {
        Theme::light()
    }
}

impl Theme {
    /// OpenShape's light theme.
    pub fn light() -> Theme {
        Theme {
            chrome: hex(0xEEF0F3),
            panel: hex(0xFFFFFF),
            panel_header: hex(0xF6F7F9),
            border: hex(0xD5D9DF),
            border_strong: hex(0xAEB4BD),
            text: hex(0x1E2329),
            text_muted: hex(0x5B6470),
            text_disabled: hex(0xA3A9B2),
            accent: hex(0x2B6CD4),
            accent_hover: hex(0x1F5BBF),
            accent_soft: hex(0xE1ECFB),
            on_accent: hex(0xFFFFFF),
            hover_soft: hex(0xEEF3FA),
            selection: hex(0x2B6CD4),
            preselect: hex(0x6FA3EE),
            success: hex(0x2E9A57),
            warning: hex(0xC98112),
            danger: hex(0xCC3D33),
            viewport_top: hex(0xFFFFFF),
            viewport_bottom: hex(0xD9E0E8),
            sketch: SketchColors {
                constrained: hex(0x1E2329),
                free: hex(0x2F62C4),
                conflict: hex(0xCC3D33),
                redundant: hex(0xC98112),
                inactive: hex(0x55627A),
                dimension: hex(0x24407E),
                glyph: hex(0x4A525E),
                label_fill: Color32::from_rgba_unmultiplied(0xFF, 0xFF, 0xFF, 0xE0),
                halo: Color32::from_rgba_unmultiplied(0x2B, 0x6C, 0xD4, 0x40),
                region: Color32::from_rgba_unmultiplied(0x2B, 0x6C, 0xD4, 0x55),
                region_hover: Color32::from_rgba_unmultiplied(0x2B, 0x6C, 0xD4, 0x28),
                axis_x: hex(0xC8423A),
                axis_y: hex(0x2E9A57),
                plane: Color32::from_rgba_unmultiplied(0x7E, 0x9C, 0xC8, 0x18),
                plane_edge: Color32::from_rgba_unmultiplied(0x6A, 0x86, 0xB4, 0x90),
                plane_active: Color32::from_rgba_unmultiplied(0x2B, 0x6C, 0xD4, 0x30),
                plane_active_edge: hex(0x2B6CD4),
                plane_label: hex(0x5B6470),
                comb: Color32::from_rgba_unmultiplied(0xB0, 0x4A, 0x9A, 0x90),
                comb_envelope: hex(0x9A3486),
            },

            rounding: 3.0,
            item_spacing: Vec2::new(6.0, 4.0),
            button_padding: Vec2::new(6.0, 3.0),
            interact_height: 24.0,
            panel_margin: 6.0,
            body_size: 13.0,
            small_size: 11.5,
            heading_size: 14.0,
            toolbar_button: 28.0,
            toolbar_icon: 20.0,
        }
    }

    /// Installs the theme on `ctx`: selects egui's light theme, replaces its
    /// light [`egui::Style`] with [`Theme::style`] and sets the fonts from
    /// [`font_definitions`]. Takes effect from the next pass.
    pub fn apply(&self, ctx: &egui::Context) {
        ctx.set_fonts(font_definitions());
        ctx.set_theme(egui::ThemePreference::Light);
        ctx.set_style_of(egui::Theme::Light, self.style());
    }

    /// The egui style for this theme (light visuals, spacing, text sizes).
    pub fn style(&self) -> egui::Style {
        let mut style = egui::Style {
            visuals: self.visuals(),
            ..egui::Style::default()
        };

        style.text_styles = [
            (
                TextStyle::Small,
                FontId::new(self.small_size, FontFamily::Proportional),
            ),
            (
                TextStyle::Body,
                FontId::new(self.body_size, FontFamily::Proportional),
            ),
            (
                TextStyle::Button,
                FontId::new(self.body_size, FontFamily::Proportional),
            ),
            (
                TextStyle::Heading,
                FontId::new(self.heading_size, FontFamily::Proportional),
            ),
            (
                TextStyle::Monospace,
                FontId::new(self.body_size - 0.5, FontFamily::Monospace),
            ),
        ]
        .into();

        let s = &mut style.spacing;
        s.item_spacing = self.item_spacing;
        s.button_padding = self.button_padding;
        s.interact_size = Vec2::new(36.0, self.interact_height);
        s.window_margin = Margin::same(round_i8(self.panel_margin));
        s.menu_margin = Margin::same(round_i8(self.panel_margin * 0.5));
        s.indent = 16.0;
        s.icon_width = 13.0;
        s.icon_width_inner = 7.0;
        s.icon_spacing = 5.0;
        s.combo_height = 240.0;
        style
    }

    /// The egui visuals for this theme.
    pub fn visuals(&self) -> egui::Visuals {
        let radius = CornerRadius::same(round_u8(self.rounding));
        let thin = |color| Stroke::new(1.0, color);
        let mut v = egui::Visuals::light();

        v.dark_mode = false;
        v.override_text_color = None;
        v.weak_text_color = Some(self.text_muted);
        v.hyperlink_color = self.accent;
        v.faint_bg_color = self.panel_header;
        v.extreme_bg_color = self.panel;
        v.text_edit_bg_color = Some(self.panel);
        v.code_bg_color = self.chrome;
        v.warn_fg_color = self.warning;
        v.error_fg_color = self.danger;

        v.panel_fill = self.panel;
        v.window_fill = self.panel;
        v.window_stroke = thin(self.border);
        v.window_corner_radius = radius;
        v.menu_corner_radius = radius;
        v.window_shadow = Shadow {
            offset: [0, 2],
            blur: 10,
            spread: 0,
            color: Color32::from_black_alpha(30),
        };
        v.popup_shadow = Shadow {
            offset: [0, 2],
            blur: 6,
            spread: 0,
            color: Color32::from_black_alpha(24),
        };
        v.text_cursor.stroke = Stroke::new(1.5, self.accent);

        v.selection.bg_fill = self.accent_soft;
        v.selection.stroke = thin(self.accent);

        let w = &mut v.widgets;
        // Labels, separators, frames.
        w.noninteractive.bg_fill = self.panel;
        w.noninteractive.weak_bg_fill = self.panel;
        w.noninteractive.bg_stroke = thin(self.border);
        w.noninteractive.fg_stroke = thin(self.text);
        w.noninteractive.corner_radius = radius;
        w.noninteractive.expansion = 0.0;
        // Buttons, check boxes, sliders at rest: white with a thin border.
        w.inactive.bg_fill = self.panel;
        w.inactive.weak_bg_fill = self.panel;
        w.inactive.bg_stroke = thin(self.border);
        w.inactive.fg_stroke = thin(self.text);
        w.inactive.corner_radius = radius;
        w.inactive.expansion = 0.0;
        // Under the pointer.
        w.hovered.bg_fill = self.hover_soft;
        w.hovered.weak_bg_fill = self.hover_soft;
        w.hovered.bg_stroke = thin(self.border_strong);
        w.hovered.fg_stroke = thin(self.text);
        w.hovered.corner_radius = radius;
        w.hovered.expansion = 0.0;
        // Being pressed or dragged.
        w.active.bg_fill = self.accent_soft;
        w.active.weak_bg_fill = self.accent_soft;
        w.active.bg_stroke = thin(self.accent);
        w.active.fg_stroke = thin(self.text);
        w.active.corner_radius = radius;
        w.active.expansion = 0.0;
        // Open combo boxes and menus.
        w.open.bg_fill = self.hover_soft;
        w.open.weak_bg_fill = self.hover_soft;
        w.open.bg_stroke = thin(self.border_strong);
        w.open.fg_stroke = thin(self.text);
        w.open.corner_radius = radius;
        w.open.expansion = 0.0;

        v.button_frame = true;
        v.collapsing_header_frame = false;
        v.indent_has_left_vline = true;
        v.striped = false;
        v.slider_trailing_fill = true;
        v
    }

    /// Frame for side, top and bottom panels: panel fill and the theme's
    /// panel margin (egui's default panel frames use their own margins).
    pub fn panel_frame(&self) -> egui::Frame {
        egui::Frame::new()
            .fill(self.panel)
            .inner_margin(Margin::same(round_i8(self.panel_margin)))
    }
}

fn round_i8(v: f32) -> i8 {
    v.round().clamp(i8::MIN as f32, i8::MAX as f32) as i8
}

fn round_u8(v: f32) -> u8 {
    v.round().clamp(0.0, u8::MAX as f32) as u8
}

/// The app's fonts: egui's defaults with the monospace face `"Hack"`
/// appended as the last fallback of the proportional family (when that font
/// is built in), so proportional text can show the constraint glyphs U+2016
/// (double vertical line, "parallel") and U+27C2 (perpendicular), which the
/// default proportional faces lack.
pub fn font_definitions() -> egui::FontDefinitions {
    let mut fonts = egui::FontDefinitions::default();
    if fonts.font_data.contains_key("Hack") {
        let proportional = fonts.families.entry(FontFamily::Proportional).or_default();
        if !proportional.iter().any(|name| name == "Hack") {
            proportional.push("Hack".to_owned());
        }
    }
    fonts
}
