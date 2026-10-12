//! The command palette: a floating search box that lists the app's
//! commands, narrows them down as the user types and runs the chosen one,
//! and the top-bar [`search_field`] that opens it.
//!
//! The application owns the commands. It describes them as a slice of
//! [`PaletteEntry`] each frame, keeps a [`CommandPalette`] while the palette
//! is open and calls [`CommandPalette::show`], which reports a
//! [`PaletteOutcome`]: still open, run this entry, or closed.
//!
//! # Look
//!
//! A white panel ([`Theme::panel`]) with a thin [`Theme::border`], the
//! theme's corner radius and its window shadow, in three bands:
//!
//! - the query field: the Search icon in [`Theme::text_muted`] and the text
//!   at [`Theme::heading_size`];
//! - the list: one row per command, [`Theme::toolbar_button`] tall, with
//!   its icon, its label (the characters the query matched in
//!   [`Theme::accent`], the rest in [`Theme::text`]), its detail in
//!   [`Theme::text_muted`] at [`Theme::small_size`], and its shortcut
//!   right-aligned in a small key cap. The highlighted row has an
//!   [`Theme::accent_soft`] background and an accent icon. Disabled
//!   commands are listed after the enabled ones in
//!   [`Theme::text_disabled`]. At most [`VISIBLE_ROWS`] rows show at once;
//!   longer lists scroll;
//! - a footer on [`Theme::panel_header`] with the keys: arrows to move,
//!   Enter to run, Esc to close.
//!
//! Every colour, radius, spacing and text size comes from [`Theme`].
//!
//! # Matching
//!
//! [`match_score`] is a fuzzy subsequence matcher. The query's whitespace is
//! dropped, and its remaining characters must appear in the text in order,
//! not necessarily next to each other ("exs" matches "Exit sketch").
//! Characters compare case-insensitively (each character is lowered on its
//! own; characters whose lower case is several characters compare as
//! themselves).
//!
//! Among all the ways the query can be laid onto the text, the one with the
//! highest score counts. A layout scores, per matched character:
//!
//! | rule | points |
//! |---|---|
//! | every matched character | +16 |
//! | it is the first character of the text | +32 |
//! | it starts a word (follows a character that is not a letter or digit: a space, `-`, `_`, …) | +20 |
//! | it is an upper-case letter after a lower-case one (`camelCase` hump) | +16 |
//! | it directly follows the previous matched character | +24 |
//! | each text character skipped between two matched characters | −2 |
//! | each text character before the first match (at most 8 in all) | −1 |
//!
//! So a prefix ("ext" in "Extrude": 48 + 40 + 40 = 128) beats the same
//! letters inside a word ("ext" in "Next feature": 15 + 40 + 40 = 95), and a
//! contiguous run beats letters scattered over the text. An empty query
//! scores 0 against every text.
//!
//! The best layout is found by dynamic programming over (query character,
//! text position), keeping for each cell the best score of a layout that
//! ends there; this is the classic alignment recurrence of Smith & Waterman
//! (*Identification of common molecular subsequences*, J. Mol. Biol. 147,
//! 1981) restricted to matches only, with a linear gap penalty. It takes
//! O(query × text) time. Ties between layouts prefer extending a run, then
//! the earliest position, so scores and highlighted positions are
//! deterministic.
//!
//! [`CommandPalette::ranked`] scores each entry by its label and its
//! keywords (the best of those) and sorts stably, so equal scores keep the
//! order the application gave.
//!
//! # Interaction
//!
//! While open, the query field holds the keyboard focus (it is asked for
//! again in any frame it lacks it). Typing filters the list and moves the
//! highlight back to the first row; Up and Down (and Page Up and Page Down,
//! by [`VISIBLE_ROWS`]) move the highlight; moving the pointer over a row
//! highlights it. Enter or a click runs the highlighted or clicked row if
//! it is enabled. Escape, or a primary press outside the panel, closes the
//! palette; a press in the very frame the palette opens is ignored, so the
//! click that opened it cannot close it again.
//!
//! The panel is an [`egui::Area`] on the [`egui::Order::Foreground`] layer.
//! Like every area, it is laid out invisibly in the first frame it is ever
//! shown (egui's sizing pass) and drawn from the next.

use crate::icons::{paint_icon, Icon};
use crate::theme::Theme;
use egui::emath::GuiRounding as _;
use egui::style::{ScrollAnimation, ScrollStyle};
use egui::text::{LayoutJob, TextFormat, TextWrapping};
use egui::text_edit::TextEditState;
use egui::{
    pos2, vec2, Align, Align2, Area, CornerRadius, CursorIcon, Event, EventFilter, FontId, Frame,
    Galley, Id, Key, Layout, Margin, Modifiers, Order, Painter, PointerButton, Pos2, Rect,
    Response, ScrollArea, Sense, Shape, Stroke, StrokeKind, TextEdit, Ui, UiBuilder, Vec2,
    WidgetInfo, WidgetType,
};
use std::sync::Arc;

/// Rows the list shows at once; longer lists scroll.
pub const VISIBLE_ROWS: usize = 10;

/// Placeholder of the palette's query field.
pub const QUERY_HINT: &str = "Search commands";

/// What the list shows when no entry matches the query.
pub const NO_MATCHES: &str = "No matching commands";

/// The narrowest palette, in toolbar buttons ([`Theme::toolbar_button`]).
const MIN_WIDTH_BUTTONS: f32 = 8.0;

/// The palette's area id (one palette is open at a time).
const AREA_ID: &str = "opendatum_command_palette";

// Scoring weights; the table in the module documentation explains them.
const SCORE_MATCH: i32 = 16;
const BONUS_TEXT_START: i32 = 32;
const BONUS_WORD_START: i32 = 20;
const BONUS_CASE_CHANGE: i32 = 16;
const BONUS_CONSECUTIVE: i32 = 24;
const PENALTY_GAP: i32 = 2;
const PENALTY_LEADING: i32 = 1;
const MAX_LEADING_PENALTY: i32 = 8;

/// One command the palette offers.
#[derive(Clone, Copy, Debug)]
pub struct PaletteEntry<'a> {
    /// Name shown and matched, e.g. "Extrude".
    pub label: &'a str,
    /// Muted text after the name, e.g. the command's group "Feature"; may be empty.
    pub detail: &'a str,
    /// Icon drawn before the label.
    pub icon: Option<Icon>,
    /// Shortcut shown right-aligned, e.g. "Shift+E" or "Ctrl+Z".
    pub shortcut: Option<&'a str>,
    /// Extra words that also match (e.g. "pad", "boss" for Extrude); never shown.
    pub keywords: &'a [&'a str],
    /// Disabled entries are listed (greyed, after enabled ones) but cannot run.
    pub enabled: bool,
}

/// What happened in one frame of the palette.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PaletteOutcome {
    /// Still open; nothing chosen yet.
    Open,
    /// The user chose this entry (an index into the `entries` slice); the palette should be closed by the caller.
    Run(usize),
    /// The user dismissed it (Escape, or a click outside).
    Closed,
}

/// Match quality of `query` against `text`: `None` if the query's characters
/// (ignoring case and whitespace) are not a subsequence of `text`, otherwise a
/// score where higher is better. Deterministic. Empty query: `Some(0)`.
///
/// The scoring rules are in the [module documentation](self#matching).
pub fn match_score(query: &str, text: &str) -> Option<i32> {
    best_match(query, text).map(|(score, _)| score)
}

/// The char indices of `text` that `query` matched for the best score
/// (for highlighting), or `None` when it does not match.
///
/// The indices are strictly increasing, one per non-whitespace character of
/// the query; an empty query matches with no indices.
pub fn match_positions(query: &str, text: &str) -> Option<Vec<usize>> {
    best_match(query, text).map(|(_, positions)| positions)
}

/// Lower case of one character, as one character: the comparison key of
/// the matcher. Characters whose lower case is several characters (such as
/// U+0130) compare as themselves, so text indices stay one-to-one.
fn fold(c: char) -> char {
    if c.is_ascii() {
        return c.to_ascii_lowercase();
    }
    let mut lower = c.to_lowercase();
    match (lower.next(), lower.next()) {
        (Some(l), None) => l,
        _ => c,
    }
}

/// The bonus for matching the character at `text[j]`, from its place in
/// the text.
fn position_bonus(text: &[char], j: usize) -> i32 {
    if j == 0 {
        return BONUS_TEXT_START;
    }
    let (prev, cur) = (text[j - 1], text[j]);
    if !prev.is_alphanumeric() && cur.is_alphanumeric() {
        BONUS_WORD_START
    } else if prev.is_lowercase() && cur.is_uppercase() {
        BONUS_CASE_CHANGE
    } else {
        0
    }
}

/// The best score of `query` laid onto `text` and the text positions of
/// that layout (see the module documentation).
fn best_match(query: &str, text: &str) -> Option<(i32, Vec<usize>)> {
    let q: Vec<char> = query
        .chars()
        .filter(|c| !c.is_whitespace())
        .map(fold)
        .collect();
    if q.is_empty() {
        return Some((0, Vec::new()));
    }
    let t: Vec<char> = text.chars().collect();
    let tf: Vec<char> = t.iter().map(|&c| fold(c)).collect();
    let (m, n) = (q.len(), t.len());
    {
        // Cheap rejection before the quadratic part.
        let mut rest = tf.iter();
        if !q.iter().all(|c| rest.any(|d| d == c)) {
            return None;
        }
    }
    let bonus: Vec<i32> = (0..n).map(|j| position_bonus(&t, j)).collect();

    // score[i * n + j]: best score of a layout of q[..=i] whose q[i] sits on
    // text[j]; from[i * n + j]: where q[i - 1] sits in that layout.
    let mut score: Vec<Option<i32>> = vec![None; m * n];
    let mut from: Vec<usize> = vec![0; m * n];
    for (j, (&c, &b)) in tf.iter().zip(&bonus).enumerate() {
        if c == q[0] {
            let leading = (j as i32).saturating_mul(PENALTY_LEADING);
            score[j] = Some(SCORE_MATCH + b - leading.min(MAX_LEADING_PENALTY));
        }
    }
    for i in 1..m {
        let (done, rest) = score.split_at_mut(i * n);
        let prev = &done[(i - 1) * n..];
        let cur = &mut rest[..n];
        // Best predecessor k <= j - 2 (at least one skipped character),
        // its gap penalty already applied for the current j.
        let mut gapped: Option<(i32, usize)> = None;
        for j in 0..n {
            if j >= i && tf[j] == q[i] {
                let run = if j >= 1 {
                    prev[j - 1].map(|s| (s + BONUS_CONSECUTIVE, j - 1))
                } else {
                    None
                };
                let best = match (run, gapped) {
                    (Some(r), Some(g)) => Some(if g.0 > r.0 { g } else { r }),
                    (r, g) => r.or(g),
                };
                if let Some((s, k)) = best {
                    cur[j] = Some(SCORE_MATCH + bonus[j] + s);
                    from[i * n + j] = k;
                }
            }
            // Step to j + 1: every gap grows by one character, and k = j - 1
            // becomes a gapped predecessor (skipping text[j]).
            gapped = gapped.map(|(s, k)| (s - PENALTY_GAP, k));
            if j >= 1 {
                if let Some(s) = prev[j - 1] {
                    let s = s - PENALTY_GAP;
                    if gapped.is_none_or(|(g, _)| s > g) {
                        gapped = Some((s, j - 1));
                    }
                }
            }
        }
    }

    let mut best: Option<(usize, i32)> = None;
    for (j, s) in score[(m - 1) * n..].iter().enumerate() {
        if let Some(s) = *s {
            if best.is_none_or(|(_, b)| s > b) {
                best = Some((j, s));
            }
        }
    }
    let (mut j, best) = best?;
    let mut positions = vec![0; m];
    for i in (0..m).rev() {
        positions[i] = j;
        if i > 0 {
            j = from[i * n + j];
        }
    }
    Some((best, positions))
}

/// An entry's score: the best of its label's and its keywords' scores.
fn entry_score(query: &str, entry: &PaletteEntry<'_>) -> Option<i32> {
    std::iter::once(entry.label)
        .chain(entry.keywords.iter().copied())
        .filter_map(|text| match_score(query, text))
        .max()
}

/// An open command palette's state (create one when the palette opens; drop it when it closes).
#[derive(Clone, Debug, Default)]
pub struct CommandPalette {
    /// The text typed so far (the caller may preset it).
    pub query: String,
    /// Highlighted row, as a position in [`CommandPalette::ranked`].
    highlight: usize,
    /// The query the highlight belongs to; `None` before the first frame.
    seen_query: Option<String>,
    /// Scroll the highlighted row into view in the next laid-out frame.
    scroll_to_highlight: bool,
    /// Frames shown so far.
    frames: u64,
    /// The panel's rectangle in the last frame.
    panel: Option<Rect>,
    /// Visible part of each listed row in the last frame, by entry index.
    rows: Vec<(usize, Rect)>,
}

/// Navigation keys taken from the input before the query field sees them.
#[derive(Clone, Copy, Debug, Default)]
struct Keys {
    escape: bool,
    enter: bool,
    down: usize,
    up: usize,
}

/// What [`CommandPalette::contents`] reports back to [`CommandPalette::show`].
struct Contents {
    ranked: Vec<usize>,
    clicked: Option<usize>,
}

impl CommandPalette {
    /// A palette with an empty query and the first row highlighted.
    pub fn new() -> Self {
        Self::default()
    }

    /// Indices of the entries to list for the current query, best first:
    /// matching enabled entries by descending score (ties keep the given
    /// order), then matching disabled ones the same way. An entry matches
    /// through its label or any keyword (its score is the best of those).
    /// With an empty query every entry is listed in the given order (enabled first).
    pub fn ranked(&self, entries: &[PaletteEntry<'_>]) -> Vec<usize> {
        let mut hits: Vec<(bool, i32, usize)> = entries
            .iter()
            .enumerate()
            .filter_map(|(i, e)| entry_score(&self.query, e).map(|s| (e.enabled, s, i)))
            .collect();
        hits.sort_by(|a, b| b.0.cmp(&a.0).then(b.1.cmp(&a.1)).then(a.2.cmp(&b.2)));
        hits.into_iter().map(|(_, _, i)| i).collect()
    }

    /// The highlighted row: a position in [`CommandPalette::ranked`] (0 is
    /// the first row).
    pub fn highlighted(&self) -> usize {
        self.highlight
    }

    /// The panel's rectangle (without its shadow) in the last frame shown,
    /// in points; `None` before the first frame.
    pub fn panel_rect(&self) -> Option<Rect> {
        self.panel
    }

    /// The visible part of the row listing `entry` (an index into the
    /// `entries` slice) in the last frame shown, in points; `None` if that
    /// entry was not listed or was scrolled out of view.
    pub fn row_rect(&self, entry: usize) -> Option<Rect> {
        self.rows.iter().find(|(e, _)| *e == entry).map(|(_, r)| *r)
    }

    /// Draws the palette as a floating panel (an egui Area on the Foreground order)
    /// of width `width` whose top edge is centred horizontally on `anchor`.
    ///
    /// The query field at the top has the keyboard focus while open (it is
    /// requested in every frame it lacks it); typing filters the list. Up and
    /// Down move the highlighted row one step (Page Up and Page Down by
    /// [`VISIBLE_ROWS`]), stopping at the ends; the highlight goes back to
    /// the first row when the query changes. Enter runs the highlighted row
    /// if it is enabled, and clicking an enabled row runs it
    /// ([`PaletteOutcome::Run`]). Escape, or a primary press outside the
    /// panel, returns [`PaletteOutcome::Closed`].
    ///
    /// Rows show the icon, the label with the characters the query matched
    /// in [`Theme::accent`] and the rest in [`Theme::text`], the detail and
    /// the right-aligned shortcut muted. At most [`VISIBLE_ROWS`] rows show;
    /// more scroll. With no match the list says [`NO_MATCHES`].
    ///
    /// The palette takes Up, Down, Page Up, Page Down, Enter and Escape out of
    /// the frame's input, so widgets drawn after it do not see them. The
    /// width is limited to the screen less a [`Theme::panel_margin`] on each
    /// side, but is at least eight [`Theme::toolbar_button`]s. See the
    /// [module documentation](self#interaction) for the details.
    pub fn show(
        &mut self,
        ctx: &egui::Context,
        theme: &Theme,
        anchor: Pos2,
        width: f32,
        entries: &[PaletteEntry<'_>],
    ) -> PaletteOutcome {
        let keys = ctx.input_mut(|i| {
            let page_down = i.count_and_consume_key(Modifiers::NONE, Key::PageDown);
            let page_up = i.count_and_consume_key(Modifiers::NONE, Key::PageUp);
            Keys {
                escape: i.consume_key(Modifiers::NONE, Key::Escape),
                enter: i.consume_key(Modifiers::NONE, Key::Enter),
                down: i.count_and_consume_key(Modifiers::NONE, Key::ArrowDown)
                    + page_down * VISIBLE_ROWS,
                up: i.count_and_consume_key(Modifiers::NONE, Key::ArrowUp) + page_up * VISIBLE_ROWS,
            }
        });
        if keys.escape {
            return PaletteOutcome::Closed;
        }
        let first_frame = self.frames == 0;
        self.frames += 1;

        let area_id = Id::new(AREA_ID);
        if first_frame {
            // Forget the cursor of an earlier palette: start at the end of
            // the (possibly preset) query.
            TextEditState::default().store(ctx, query_id(area_id));
        }
        let screen = ctx.content_rect();
        let width = width
            .min(screen.width() - 2.0 * theme.panel_margin)
            .max(MIN_WIDTH_BUTTONS * theme.toolbar_button);
        let area = Area::new(area_id)
            .order(Order::Foreground)
            .pivot(Align2::CENTER_TOP)
            .fixed_pos(anchor)
            .fade_in(false)
            .default_width(width)
            .show(ctx, |ui| self.contents(ui, theme, width, entries, keys));
        let Contents { ranked, clicked } = area.inner;
        let panel = area.response.rect;
        self.panel = Some(panel);

        if let Some(entry) = clicked {
            return PaletteOutcome::Run(entry);
        }
        if keys.enter {
            if let Some(&entry) = ranked.get(self.highlight) {
                if entries[entry].enabled {
                    return PaletteOutcome::Run(entry);
                }
            }
        }
        let pressed_outside = ctx.input(|i| {
            i.events.iter().any(|e| {
                matches!(
                    e,
                    Event::PointerButton {
                        pos,
                        button: PointerButton::Primary,
                        pressed: true,
                        ..
                    } if !panel.contains(*pos)
                )
            })
        });
        if pressed_outside && !first_frame {
            return PaletteOutcome::Closed;
        }
        PaletteOutcome::Open
    }

    /// Lays out and paints the panel inside the palette's area.
    fn contents(
        &mut self,
        ui: &mut Ui,
        theme: &Theme,
        width: f32,
        entries: &[PaletteEntry<'_>],
        keys: Keys,
    ) -> Contents {
        let m = Metrics::new(theme);
        ui.style_mut().visuals = theme.visuals();
        ui.spacing_mut().item_spacing = Vec2::ZERO;
        // A narrow scroll bar, inset like the rows' highlight, always shown
        // while the list overflows and wider under the pointer. Floating, so
        // the rows keep their width while it appears.
        ui.spacing_mut().scroll = ScrollStyle {
            floating: true,
            bar_width: 2.0 * m.inset,
            floating_width: m.inset,
            floating_allocated_width: 2.0 * m.inset,
            bar_inner_margin: 0.0,
            bar_outer_margin: m.inset,
            handle_min_length: m.row_height,
            foreground_color: false,
            dormant_background_opacity: 0.0,
            active_background_opacity: 0.0,
            interact_background_opacity: 0.0,
            dormant_handle_opacity: 1.0,
            active_handle_opacity: 1.0,
            interact_handle_opacity: 1.0,
            ..ScrollStyle::solid()
        };
        ui.set_width(width);
        // Shadow and fill go under everything; their size is known at the end.
        let background = ui.painter().add(Shape::Noop);

        // ---- Query field ----
        let (header, _) = ui.allocate_exact_size(vec2(width, m.header_height), Sense::hover());
        let icon = m.icon_rect(header.left() + m.content_left, header.center().y);
        paint_icon(ui.painter(), icon, Icon::Search, theme.text_muted);
        let field = Rect::from_min_max(
            pos2(icon.right() + m.gap, header.top()),
            pos2(header.right() - m.content_left, header.bottom()),
        );
        let id = query_id(Id::new(AREA_ID));
        if !ui.memory(|mem| mem.has_focus(id)) {
            ui.memory_mut(|mem| mem.request_focus(id));
        }
        let mut field_ui = ui.new_child(
            UiBuilder::new()
                .max_rect(field)
                .layout(Layout::left_to_right(Align::Center)),
        );
        field_ui.add(
            TextEdit::singleline(&mut self.query)
                .id(id)
                .frame(Frame::NONE)
                .margin(Margin::ZERO)
                .font(FontId::proportional(theme.heading_size))
                .text_color(theme.text)
                .hint_text(QUERY_HINT)
                .desired_width(field.width())
                .return_key(None)
                // Keep the focus: arrows, Tab and Escape belong to the palette.
                .event_filter(EventFilter {
                    tab: true,
                    horizontal_arrows: true,
                    vertical_arrows: true,
                    escape: true,
                }),
        );
        let (rule, _) = ui.allocate_exact_size(vec2(width, 1.0), Sense::hover());
        ui.painter().rect_filled(rule, 0.0, theme.border);
        // The list's scroll bar handle: the strong border colour, darker under
        // the pointer (egui would use widget fills, the panel's white here).
        {
            let v = ui.visuals_mut();
            v.widgets.inactive.bg_fill = theme.border_strong;
            v.widgets.hovered.bg_fill = theme.text_muted;
            v.widgets.active.bg_fill = theme.text_muted;
        }

        // ---- Highlight ----
        let ranked = self.ranked(entries);
        if self.seen_query.as_deref() != Some(self.query.as_str()) {
            self.seen_query = Some(self.query.clone());
            self.highlight = 0;
            self.scroll_to_highlight = true;
        }
        let last = ranked.len().saturating_sub(1);
        if keys.down > 0 || keys.up > 0 {
            self.highlight = self
                .highlight
                .saturating_add(keys.down)
                .min(last)
                .saturating_sub(keys.up);
            self.scroll_to_highlight = true;
        }
        self.highlight = self.highlight.min(last);

        // ---- List ----
        let mut clicked = None;
        self.rows.clear();
        if ranked.is_empty() {
            let (rect, _) =
                ui.allocate_exact_size(vec2(width, m.row_height + 2.0 * m.inset), Sense::hover());
            ui.painter().text(
                rect.center(),
                Align2::CENTER_CENTER,
                NO_MATCHES,
                FontId::proportional(theme.body_size),
                theme.text_muted,
            );
        } else {
            let pointer_moved = ui.input(|i| i.pointer.delta() != Vec2::ZERO);
            ScrollArea::vertical()
                .id_salt(AREA_ID)
                .max_height(VISIBLE_ROWS as f32 * m.row_height + 2.0 * m.inset)
                .auto_shrink([false, true])
                // Keyboard navigation jumps; it never animates.
                .animated(false)
                .show(ui, |ui| {
                    ui.add_space(m.inset);
                    let row_width = ui.available_width();
                    let rows: Vec<(Rect, Response)> = ranked
                        .iter()
                        .map(|_| {
                            ui.allocate_exact_size(vec2(row_width, m.row_height), Sense::CLICK)
                        })
                        .collect();
                    ui.add_space(m.inset);

                    if pointer_moved {
                        if let Some(row) = rows.iter().position(|(_, r)| r.hovered()) {
                            self.highlight = row;
                        }
                    }
                    if self.scroll_to_highlight && !ui.is_sizing_pass() {
                        self.scroll_to_highlight = false;
                        if let Some((rect, _)) = rows.get(self.highlight) {
                            ui.scroll_to_rect_animation(
                                rect.expand2(vec2(0.0, m.inset)),
                                None,
                                ScrollAnimation::none(),
                            );
                        }
                    }

                    let clip = ui.clip_rect();
                    let style = RowStyle {
                        theme,
                        metrics: &m,
                        query: &self.query,
                    };
                    for (row, ((rect, response), &entry)) in rows.iter().zip(&ranked).enumerate() {
                        let e = &entries[entry];
                        let highlighted = row == self.highlight;
                        response.widget_info(|| {
                            WidgetInfo::selected(
                                WidgetType::SelectableLabel,
                                e.enabled,
                                highlighted,
                                e.label,
                            )
                        });
                        if e.enabled && response.clicked() {
                            clicked = Some(entry);
                        }
                        let visible = rect.intersect(clip);
                        if visible.is_positive() {
                            self.rows.push((entry, visible));
                        }
                        if ui.is_rect_visible(*rect) {
                            paint_row(ui.painter(), &style, *rect, e, highlighted);
                        }
                    }
                });
        }

        // ---- Footer ----
        let (footer, _) = ui.allocate_exact_size(vec2(width, m.footer_height), Sense::hover());
        paint_footer(&ui.painter_at(footer), theme, &m, footer);

        // ---- Panel ----
        let panel = ui.min_rect();
        let shadow = theme.visuals().window_shadow.as_shape(panel, m.radius);
        ui.painter().set(
            background,
            Shape::Vec(vec![
                shadow.into(),
                Shape::rect_filled(panel, m.radius, theme.panel),
            ]),
        );
        ui.painter().rect_stroke(
            panel,
            m.radius,
            Stroke::new(1.0, theme.border),
            StrokeKind::Inside,
        );
        Contents { ranked, clicked }
    }
}

/// The id of the palette's query field.
fn query_id(area_id: Id) -> Id {
    area_id.with("query")
}

/// Sizes of the palette's parts, all derived from [`Theme`] tokens.
#[derive(Clone, Copy, Debug)]
struct Metrics {
    /// Height of a list row: a toolbar button's side.
    row_height: f32,
    /// Height of the query field band.
    header_height: f32,
    /// Height of the key-hint footer.
    footer_height: f32,
    /// Gap between the panel's edge and a row's highlight (also the list's
    /// top and bottom padding).
    inset: f32,
    /// Distance from the panel's left edge to the icons.
    content_left: f32,
    /// Side of the icons in rows and fields.
    icon: f32,
    /// Gap between an icon and its text, and between neighbouring parts.
    gap: f32,
    /// Corner radius of the panel, highlights and key caps.
    radius: CornerRadius,
}

impl Metrics {
    fn new(t: &Theme) -> Self {
        let inset = t.panel_margin * 0.5;
        Metrics {
            row_height: t.toolbar_button,
            header_height: t.toolbar_button + 2.0 * t.item_spacing.y,
            footer_height: t.interact_height,
            inset,
            content_left: inset + t.button_padding.x,
            icon: t.interact_height - 2.0 * t.item_spacing.y,
            gap: t.item_spacing.x * 1.5,
            radius: CornerRadius::same(t.rounding.round().clamp(0.0, 255.0) as u8),
        }
    }

    /// An icon's square with its left edge at `left`, centred on `center_y`.
    fn icon_rect(&self, left: f32, center_y: f32) -> Rect {
        Rect::from_min_size(
            pos2(left, center_y - self.icon * 0.5),
            Vec2::splat(self.icon),
        )
    }
}

/// What every row of one frame is painted with.
struct RowStyle<'a> {
    theme: &'a Theme,
    metrics: &'a Metrics,
    query: &'a str,
}

/// Paints one list row: highlight, icon, label (matched characters
/// accented), detail and shortcut key cap.
fn paint_row(
    painter: &Painter,
    style: &RowStyle<'_>,
    rect: Rect,
    entry: &PaletteEntry<'_>,
    highlighted: bool,
) {
    let (t, m) = (style.theme, style.metrics);
    if highlighted {
        painter.rect_filled(rect.shrink2(vec2(m.inset, 0.0)), m.radius, t.accent_soft);
    }
    let (ink, matched, muted, icon_ink) = if !entry.enabled {
        (
            t.text_disabled,
            t.text_muted,
            t.text_disabled,
            t.text_disabled,
        )
    } else if highlighted {
        (t.text, t.accent, t.text_muted, t.accent)
    } else {
        (t.text, t.accent, t.text_muted, t.text)
    };

    let icon = m.icon_rect(rect.left() + m.content_left, rect.center().y);
    if let Some(i) = entry.icon {
        paint_icon(painter, icon, i, icon_ink);
    }
    // The icon column is kept without an icon, so labels line up.
    let left = icon.right() + m.gap;
    let mut right = rect.right() - m.content_left;
    if let Some(shortcut) = entry.shortcut.filter(|s| !s.is_empty()) {
        let cap = Keycap::new(painter, t, shortcut, muted);
        right = cap
            .paint(painter, t, pos2(right - cap.size.x, rect.center().y))
            .left()
            - m.gap;
    }

    let mut job = LayoutJob {
        wrap: TextWrapping {
            max_width: (right - left).max(0.0),
            max_rows: 1,
            break_anywhere: true,
            overflow_character: Some('…'),
        },
        ..LayoutJob::default()
    };
    let font = FontId::proportional(t.body_size);
    let positions = match_positions(style.query, entry.label).unwrap_or_default();
    let mut run = String::new();
    let mut run_matched = false;
    for (i, c) in entry.label.chars().enumerate() {
        let is_matched = positions.binary_search(&i).is_ok();
        if is_matched != run_matched && !run.is_empty() {
            job.append(
                &run,
                0.0,
                text_format(&font, if run_matched { matched } else { ink }),
            );
            run.clear();
        }
        run_matched = is_matched;
        run.push(c);
    }
    if !run.is_empty() {
        job.append(
            &run,
            0.0,
            text_format(&font, if run_matched { matched } else { ink }),
        );
    }
    if !entry.detail.is_empty() {
        job.append(
            entry.detail,
            m.gap,
            text_format(&FontId::proportional(t.small_size), muted),
        );
    }
    let galley = painter.layout_job(job);
    let pos = pos2(left, rect.center().y - galley.size().y * 0.5);
    painter.galley(pos, galley, ink);
}

/// A text section in `font` and `color`, bottom-aligned on its line.
fn text_format(font: &FontId, color: egui::Color32) -> TextFormat {
    TextFormat {
        font_id: font.clone(),
        color,
        // Mixed sizes share the bottom of the line, close to a common baseline.
        valign: Align::BOTTOM,
        ..TextFormat::default()
    }
}

/// The footer: key caps with what they do.
fn paint_footer(painter: &Painter, theme: &Theme, m: &Metrics, rect: Rect) {
    let mut radius = m.radius;
    radius.nw = 0;
    radius.ne = 0;
    painter.rect_filled(rect, radius, theme.panel_header);
    painter.rect_filled(
        Rect::from_min_size(rect.min, vec2(rect.width(), 1.0)),
        0.0,
        theme.border,
    );
    let font = FontId::proportional(theme.small_size);
    let y = rect.center().y + 0.5;
    let mut x = rect.left() + m.content_left;
    let hints: [(&[&str], &str); 3] = [
        (&["\u{2191}", "\u{2193}"], "to move"),
        (&["Enter"], "to run"),
        (&["Esc"], "to close"),
    ];
    for (caps, what) in hints {
        for key in caps {
            let cap = Keycap::new(painter, theme, key, theme.text_muted);
            x = cap.paint(painter, theme, pos2(x, y)).right() + m.inset;
        }
        let galley = painter.layout_no_wrap(what.to_owned(), font.clone(), theme.text_muted);
        let size = galley.size();
        painter.galley(pos2(x, y - size.y * 0.5), galley, theme.text_muted);
        x += size.x + 2.0 * m.gap;
    }
}

/// A key name in a small rounded box (shortcuts and key hints).
struct Keycap {
    galley: Arc<Galley>,
    size: Vec2,
}

impl Keycap {
    fn new(painter: &Painter, theme: &Theme, text: &str, color: egui::Color32) -> Self {
        let galley = painter.layout_no_wrap(
            text.to_owned(),
            FontId::proportional(theme.small_size),
            color,
        );
        let padding = vec2(theme.button_padding.x * 0.75, theme.button_padding.y * 0.5);
        // Never narrower than tall, so single symbols get square caps.
        let mut size = (galley.size() + 2.0 * padding).round();
        size.x = size.x.max(size.y);
        Keycap { galley, size }
    }

    /// Paints the cap with its left edge at `left_center.x`, centred on
    /// `left_center.y`, snapped to whole pixels; returns its rectangle.
    fn paint(&self, painter: &Painter, theme: &Theme, left_center: Pos2) -> Rect {
        let ppp = painter.pixels_per_point();
        let rect = Rect::from_min_size(
            pos2(left_center.x, left_center.y - self.size.y * 0.5),
            self.size,
        )
        .round_to_pixels(ppp);
        painter.rect(
            rect,
            CornerRadius::same(theme.rounding.round().clamp(0.0, 255.0) as u8),
            theme.panel_header,
            Stroke::new(1.0, theme.border),
            StrokeKind::Inside,
        );
        let text = self.galley.size();
        let x = rect.center().x - text.x * 0.5;
        let y = rect.center().y - text.y * 0.5;
        painter.galley(pos2(x, y), Arc::clone(&self.galley), theme.text_muted);
        rect
    }
}

/// A search field look-alike for the top bar: a rounded box `width` wide and
/// the theme's interact height, with the Search icon, the muted `hint` text and,
/// right-aligned, the muted `shortcut` (e.g. "Ctrl+K"). It is not editable; it
/// senses clicks (the app opens the palette on click). Hover changes its border
/// to `theme.border_strong`.
///
/// It is focusable like a button (Space or Enter on it count as a click);
/// focused, its border is [`Theme::accent`]. An empty `shortcut` shows no key
/// cap. The hint is cut short with an ellipsis when the box is too narrow.
pub fn search_field(
    ui: &mut Ui,
    theme: &Theme,
    hint: &str,
    shortcut: &str,
    width: f32,
) -> Response {
    let m = Metrics::new(theme);
    let (rect, response) =
        ui.allocate_exact_size(vec2(width, theme.interact_height), Sense::click());
    let enabled = ui.is_enabled();
    response.widget_info(|| WidgetInfo::labeled(WidgetType::Button, enabled, hint));
    if ui.is_rect_visible(rect) {
        let painter = ui.painter();
        let border = if response.has_focus() {
            theme.accent
        } else if response.hovered() {
            theme.border_strong
        } else {
            theme.border
        };
        painter.rect(
            rect,
            m.radius,
            theme.panel,
            Stroke::new(1.0, border),
            StrokeKind::Inside,
        );
        let (ink, muted) = if enabled {
            (theme.text_muted, theme.text_muted)
        } else {
            (theme.text_disabled, theme.text_disabled)
        };
        let icon = m.icon_rect(rect.left() + theme.button_padding.x, rect.center().y);
        paint_icon(painter, icon, Icon::Search, ink);
        let left = icon.right() + m.gap;
        let mut right = rect.right() - theme.button_padding.x;
        if !shortcut.is_empty() {
            let cap = Keycap::new(painter, theme, shortcut, muted);
            // As far from the right edge as from the top and bottom.
            let margin = ((rect.height() - cap.size.y) * 0.5).max(0.0);
            let left_edge = rect.right() - margin - cap.size.x;
            right = cap
                .paint(painter, theme, pos2(left_edge, rect.center().y))
                .left()
                - m.gap;
        }
        let job = LayoutJob {
            wrap: TextWrapping {
                max_width: (right - left).max(0.0),
                max_rows: 1,
                break_anywhere: true,
                overflow_character: Some('…'),
            },
            ..LayoutJob::single_section(
                hint.to_owned(),
                text_format(&FontId::proportional(theme.body_size), muted),
            )
        };
        let galley = painter.layout_job(job);
        let pos = pos2(left, rect.center().y - galley.size().y * 0.5);
        painter.galley(pos, galley, muted);
    }
    if enabled {
        response.on_hover_cursor(CursorIcon::PointingHand)
    } else {
        response
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fold_is_one_to_one() {
        assert_eq!(fold('E'), 'e');
        assert_eq!(fold('É'), 'é');
        // U+0130 lowers to two characters, so it stays itself.
        assert_eq!(fold('\u{130}'), '\u{130}');
    }

    #[test]
    fn bonuses_follow_word_and_case_boundaries() {
        let t: Vec<char> = "Fit view_x-yZ".chars().collect();
        assert_eq!(position_bonus(&t, 0), BONUS_TEXT_START);
        assert_eq!(position_bonus(&t, 1), 0);
        assert_eq!(position_bonus(&t, 4), BONUS_WORD_START);
        assert_eq!(position_bonus(&t, 9), BONUS_WORD_START);
        assert_eq!(position_bonus(&t, 11), BONUS_WORD_START);
        assert_eq!(position_bonus(&t, 12), BONUS_CASE_CHANGE);
    }

    #[test]
    fn documented_examples_score_as_stated() {
        assert_eq!(match_score("ext", "Extrude"), Some(128));
        assert_eq!(match_score("ext", "Next feature"), Some(95));
    }
}
