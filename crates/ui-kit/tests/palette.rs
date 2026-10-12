//! The command palette: fuzzy matching (pure and property-tested), ranking,
//! keyboard and pointer interaction driven through an egui context, and a
//! render on the GPU (lavapipe in CI) for the palette and the search field.

mod common;

use common::{max_channel_diff, maybe_dump, rgba, snapshot};
use egui::{pos2, vec2, Event, Key, LayerId, Modifiers, PointerButton, Pos2, Rect, Vec2};
use opendatum_ui_kit::palette::{match_positions, match_score, VISIBLE_ROWS};
use opendatum_ui_kit::{search_field, CommandPalette, Icon, PaletteEntry, PaletteOutcome, Theme};
use proptest::prelude::*;

// ---------------------------------------------------------------------------
// Matching

#[test]
fn matching_ignores_case() {
    let lower = match_score("ext", "extrude");
    assert!(lower.is_some());
    assert_eq!(match_score("EXT", "extrude"), lower);
    assert_eq!(match_score("ext", "EXTRUDE"), lower);
    assert_eq!(match_score("eXt", "Extrude"), lower);
    assert_eq!(match_positions("eXt", "EXTRUDE"), Some(vec![0, 1, 2]));
    assert_eq!(match_positions("é", "É"), Some(vec![0]));
}

#[test]
fn a_word_start_beats_the_same_letters_mid_word() {
    let extrude = match_score("ext", "Extrude").unwrap();
    let next = match_score("ext", "Next feature").unwrap();
    assert!(extrude > next, "{extrude} vs {next}");
    // Inside a word the letters still match, as one run.
    assert_eq!(match_positions("ext", "Next feature"), Some(vec![1, 2, 3]));
    // After a space, '-' or '_', or at a lower-to-upper case change.
    for text in ["fit tool", "fit-tool", "fit_tool"] {
        assert_eq!(match_positions("t", text), Some(vec![4]), "{text}");
    }
    assert_eq!(match_positions("t", "fitTool"), Some(vec![3]));
    assert_eq!(match_positions("t", "fittool"), Some(vec![2]));
    assert_eq!(match_positions("se", "ShadedEdges"), Some(vec![0, 6]));
    assert_eq!(match_positions("fv", "Fit view"), Some(vec![0, 4]));
}

#[test]
fn an_exact_prefix_beats_a_scattered_subsequence() {
    let pairs = [
        ("rev", "Revolve", "Remove vertex"),
        ("fit", "Fit view", "Fillet tool"),
        ("box", "Box", "Bounding box"),
        ("ex", "Extrude", "Hexagon"),
        ("sketch", "Sketch", "Exit sketch"),
        ("und", "Undo", "Bounding sphere"),
    ];
    for (query, prefix, scattered) in pairs {
        let p = match_score(query, prefix).unwrap();
        let s = match_score(query, scattered).unwrap();
        assert!(p > s, "{query:?}: {prefix:?} {p} vs {scattered:?} {s}");
    }
}

#[test]
fn non_subsequences_do_not_match() {
    for (query, text) in [
        ("dx", "Extrude"),
        ("extrudes", "Extrude"),
        ("a", ""),
        ("zz", "Fit view"),
        ("vf", "Fit view"),
    ] {
        assert_eq!(match_score(query, text), None, "{query:?} in {text:?}");
        assert_eq!(match_positions(query, text), None, "{query:?} in {text:?}");
    }
}

#[test]
fn whitespace_in_the_query_is_ignored() {
    assert_eq!(
        match_score("fit view", "Fit view"),
        match_score("fitview", "Fit view")
    );
    assert_eq!(
        match_score(" e x ", "Extrude"),
        match_score("ex", "Extrude")
    );
    assert_eq!(
        match_positions("fit view", "Fit view"),
        Some(vec![0, 1, 2, 4, 5, 6, 7])
    );
    assert_eq!(match_score(" \t ", "Extrude"), Some(0));
    assert_eq!(match_positions(" \t ", "Extrude"), Some(vec![]));
}

#[test]
fn the_empty_query_matches_everything_with_score_zero() {
    for text in ["", "Extrude", "Fit view"] {
        assert_eq!(match_score("", text), Some(0));
        assert_eq!(match_positions("", text), Some(vec![]));
    }
}

#[test]
fn positions_are_char_indices_not_bytes() {
    // "É" is two bytes; positions count characters.
    assert_eq!(match_positions("ex", "Éxe x"), Some(vec![2, 4]));
    assert_eq!(match_positions("éx", "Été x"), Some(vec![0, 4]));
    // Accented and plain letters are different letters.
    assert_eq!(match_score("e", "É"), None);
}

#[test]
fn matching_is_deterministic() {
    for (query, text) in [("ex", "Export STL"), ("es", "Exit sketch"), ("e", "eee")] {
        assert_eq!(match_score(query, text), match_score(query, text));
        assert_eq!(match_positions(query, text), match_positions(query, text));
    }
    // Equal alternatives resolve to the earliest position: both x start a
    // word past the capped leading penalty.
    assert_eq!(
        match_score("x", "aaaaaaaaaa-x"),
        match_score("x", "aaaaaaaaaa-a-x")
    );
    assert_eq!(match_positions("x", "aaaaaaaaaa-x-x"), Some(vec![11]));
}

/// Case-insensitive equality of two characters, as Unicode lower-casing
/// gives it.
fn same_letter(a: char, b: char) -> bool {
    a.to_lowercase().eq(b.to_lowercase())
}

fn check_consistent(query: &str, text: &str) -> Result<(), TestCaseError> {
    let score = match_score(query, text);
    let positions = match_positions(query, text);
    prop_assert_eq!(score.is_some(), positions.is_some());
    // Deterministic.
    prop_assert_eq!(match_score(query, text), score);
    prop_assert_eq!(&match_positions(query, text), &positions);
    if let Some(positions) = positions {
        let wanted: Vec<char> = query.chars().filter(|c| !c.is_whitespace()).collect();
        let chars: Vec<char> = text.chars().collect();
        prop_assert_eq!(positions.len(), wanted.len());
        prop_assert!(positions.windows(2).all(|w| w[0] < w[1]), "{:?}", positions);
        for (&p, &q) in positions.iter().zip(&wanted) {
            prop_assert!(p < chars.len());
            prop_assert!(
                same_letter(chars[p], q),
                "{:?} at {} for {:?}",
                chars[p],
                p,
                q
            );
        }
    }
    Ok(())
}

/// Texts and queries from a mixed alphabet: both cases, digits, word
/// separators, whitespace and a few non-ASCII letters.
const ALPHABET: &str = "[a-dA-Dexyz0-9 _\\-\tÉéß]";

/// Strings of up to `max` characters from [`ALPHABET`].
fn strings(max: usize) -> impl Strategy<Value = String> {
    proptest::string::string_regex(&format!("{ALPHABET}{{0,{max}}}")).expect("valid regex")
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn score_and_positions_agree(
        query in strings(6),
        text in strings(24),
    ) {
        check_consistent(&query, &text)?;
    }

    #[test]
    fn a_subsequence_with_any_case_and_spaces_matches(
        text in strings(24),
        picks in proptest::collection::vec((any::<bool>(), any::<bool>(), any::<bool>()), 0..24),
    ) {
        // Keep some characters of the text, flip the case of some, add spaces.
        let mut query = String::new();
        for (c, (keep, upper, space)) in text.chars().zip(picks) {
            if !keep {
                continue;
            }
            let mut upper_case = c.to_uppercase();
            if let (true, Some(u), None) = (upper, upper_case.next(), upper_case.next()) {
                // Only one-to-one case pairs ("ß" upper-cases to "SS").
                query.push(u);
            } else {
                query.push(c);
            }
            if space {
                query.push(' ');
            }
        }
        prop_assert!(match_score(&query, &text).is_some(), "{:?} in {:?}", query, text);
        check_consistent(&query, &text)?;
    }

    #[test]
    fn every_text_matches_itself_and_the_empty_query(text in strings(32)) {
        prop_assert!(match_score(&text, &text).is_some());
        prop_assert_eq!(match_score("", &text), Some(0));
        prop_assert_eq!(match_positions("", &text), Some(vec![]));
        check_consistent(&text, &text)?;
    }
}

// ---------------------------------------------------------------------------
// Ranking

const fn entry(label: &'static str, enabled: bool) -> PaletteEntry<'static> {
    PaletteEntry {
        label,
        detail: "",
        icon: None,
        shortcut: None,
        keywords: &[],
        enabled,
    }
}

fn ranked(query: &str, entries: &[PaletteEntry<'_>]) -> Vec<usize> {
    let mut palette = CommandPalette::new();
    palette.query = query.to_owned();
    palette.ranked(entries)
}

#[test]
fn the_empty_query_lists_everything_enabled_first_in_the_given_order() {
    let entries = [
        entry("Undo", true),
        entry("Exit sketch", false),
        entry("Redo", true),
        entry("Finish sketch", false),
        entry("Box", true),
    ];
    assert_eq!(ranked("", &entries), vec![0, 2, 4, 1, 3]);
    assert_eq!(ranked("   ", &entries), vec![0, 2, 4, 1, 3]);
}

#[test]
fn enabled_matches_come_before_disabled_ones_whatever_their_score() {
    let entries = [
        entry("Exit sketch", false),
        entry("Next feature", true),
        entry("Extrude", true),
        entry("Ext", false),
    ];
    // The disabled "Ext" outscores the enabled "Next feature" but is listed
    // after every enabled match.
    assert!(match_score("ext", "Ext") > match_score("ext", "Next feature"));
    assert_eq!(ranked("ext", &entries), vec![2, 1, 3, 0]);
}

#[test]
fn entries_rank_by_score_and_non_matches_are_left_out() {
    let entries = [
        entry("Next feature", true),
        entry("Revolve", true),
        entry("Extrude", true),
    ];
    assert_eq!(ranked("ext", &entries), vec![2, 0]);
    assert_eq!(ranked("zzz", &entries), Vec::<usize>::new());
}

#[test]
fn keywords_match_and_the_best_of_label_and_keywords_counts() {
    let extrude = PaletteEntry {
        keywords: &["pad", "boss"],
        ..entry("Extrude", true)
    };
    let entries = [entry("Sketch", true), extrude, entry("Spline", true)];
    assert_eq!(ranked("boss", &entries), vec![1]);
    assert_eq!(ranked("pad", &entries), vec![1]);
    // "s" starts both labels but sits mid-word in the keyword "boss".
    assert_eq!(ranked("s", &entries), vec![0, 2, 1]);
    // A keyword prefix ("pa") beats a scattered label match ("S…p…a…").
    let entries = [entry("Spline array", true), extrude];
    assert_eq!(ranked("pa", &entries), vec![1, 0]);
}

#[test]
fn equal_scores_keep_the_given_order() {
    let entries = [
        entry("Extrude", true),
        entry("Export STL", true),
        entry("Exit sketch", true),
    ];
    assert_eq!(ranked("ex", &entries), vec![0, 1, 2]);
    let reversed = [entries[2], entries[1], entries[0]];
    assert_eq!(ranked("ex", &reversed), vec![0, 1, 2]);
    let dup = [entry("Box", true), entry("Box", true), entry("Box", true)];
    assert_eq!(ranked("bx", &dup), vec![0, 1, 2]);
}

// ---------------------------------------------------------------------------
// Interaction

const EXTRUDE: usize = 0;
const REVOLVE: usize = 1;
const EXIT_SKETCH: usize = 9;

/// Fifteen commands (more than fit, so the list scrolls); "Exit sketch" is
/// disabled.
fn commands() -> Vec<PaletteEntry<'static>> {
    let e = |label, detail, icon, shortcut, enabled| PaletteEntry {
        label,
        detail,
        icon: Some(icon),
        shortcut,
        keywords: &[],
        enabled,
    };
    vec![
        PaletteEntry {
            keywords: &["pad", "boss"],
            ..e("Extrude", "Feature", Icon::Extrude, Some("Shift+E"), true)
        },
        e("Revolve", "Feature", Icon::Revolve, Some("Shift+W"), true),
        e("Sketch", "Feature", Icon::Sketch, Some("Shift+S"), true),
        e("Box", "Feature", Icon::Box, None, true),
        e("Cylinder", "Feature", Icon::Cylinder, None, true),
        e("Export STL", "File", Icon::Export, None, true),
        e("Undo", "Edit", Icon::Undo, Some("Ctrl+Z"), true),
        e("Redo", "Edit", Icon::Redo, Some("Ctrl+Y"), true),
        e("Fit view", "View", Icon::ZoomFit, Some("F"), true),
        e("Exit sketch", "Sketch", Icon::FinishSketch, None, false),
        e("Shaded", "View", Icon::Shaded, None, true),
        e("Wireframe", "View", Icon::Wireframe, None, true),
        e("Home view", "View", Icon::Home, Some("H"), true),
        e("Perspective", "View", Icon::Perspective, None, true),
        e("Orthographic", "View", Icon::Orthographic, None, true),
    ]
}

const SCREEN: Vec2 = vec2(800.0, 600.0);
const ANCHOR: Pos2 = pos2(400.0, 40.0);
const WIDTH: f32 = 480.0;

/// A palette shown in a bare egui context (no GPU).
struct Harness {
    ctx: egui::Context,
    theme: Theme,
    palette: CommandPalette,
    entries: Vec<PaletteEntry<'static>>,
}

impl Harness {
    fn new() -> Self {
        let ctx = egui::Context::default();
        let theme = Theme::light();
        theme.apply(&ctx);
        Harness {
            ctx,
            theme,
            palette: CommandPalette::new(),
            entries: commands(),
        }
    }

    /// A harness after the frames egui needs to size and place the panel.
    fn settled() -> Self {
        let mut h = Harness::new();
        for _ in 0..2 {
            assert_eq!(h.frame(vec![]), PaletteOutcome::Open);
        }
        h
    }

    fn frame(&mut self, events: Vec<Event>) -> PaletteOutcome {
        let input = egui::RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, SCREEN)),
            events,
            ..Default::default()
        };
        let mut outcome = PaletteOutcome::Open;
        let (palette, theme, entries) = (&mut self.palette, &self.theme, &self.entries);
        let mut out = self.ctx.run_ui(input, |ui| {
            outcome = palette.show(ui.ctx(), theme, ANCHOR, WIDTH, entries);
        });
        out.textures_delta.clear();
        outcome
    }

    fn ranked(&self) -> Vec<usize> {
        self.palette.ranked(&self.entries)
    }

    /// Moves the pointer to `pos`, presses and releases the primary button
    /// over three frames; returns the outcomes.
    fn click(&mut self, pos: Pos2) -> [PaletteOutcome; 3] {
        let button = |pressed| Event::PointerButton {
            pos,
            button: PointerButton::Primary,
            pressed,
            modifiers: Modifiers::NONE,
        };
        [
            self.frame(vec![Event::PointerMoved(pos)]),
            self.frame(vec![button(true)]),
            self.frame(vec![button(false)]),
        ]
    }
}

fn key(key: Key) -> Event {
    Event::Key {
        key,
        physical_key: None,
        pressed: true,
        repeat: false,
        modifiers: Modifiers::NONE,
    }
}

#[test]
fn the_panel_hangs_centred_below_the_anchor() {
    let h = Harness::settled();
    let panel = h.palette.panel_rect().expect("shown");
    assert_eq!(panel.width(), WIDTH);
    assert!((panel.center().x - ANCHOR.x).abs() <= 0.5, "{panel:?}");
    assert_eq!(panel.top(), ANCHOR.y);
    // Every listed row lies inside it.
    for entry in h.ranked() {
        if let Some(r) = h.palette.row_rect(entry) {
            assert!(panel.contains_rect(r), "{entry}: {r:?} outside {panel:?}");
        }
    }
}

#[test]
fn at_most_ten_rows_show_and_the_rest_scroll() {
    let h = Harness::settled();
    let row_height = h.theme.toolbar_button;
    let ranked = h.ranked();
    assert!(ranked.len() > VISIBLE_ROWS);
    let full: Vec<usize> = ranked
        .iter()
        .copied()
        .filter(|&e| {
            h.palette
                .row_rect(e)
                .is_some_and(|r| r.height() == row_height)
        })
        .collect();
    assert_eq!(full, ranked[..VISIBLE_ROWS].to_vec());
    assert_eq!(h.palette.row_rect(*ranked.last().unwrap()), None);
}

#[test]
fn typing_filters_the_list() {
    let mut h = Harness::settled();
    assert_eq!(
        h.frame(vec![Event::Text("re".into())]),
        PaletteOutcome::Open
    );
    assert_eq!(h.frame(vec![Event::Text("v".into())]), PaletteOutcome::Open);
    assert_eq!(h.palette.query, "rev");
    let ranked = h.ranked();
    assert_eq!(ranked[0], REVOLVE);
    assert!(h.palette.row_rect(REVOLVE).is_some());
    assert_eq!(h.palette.row_rect(EXTRUDE), None);
    assert!(ranked.len() < h.entries.len());
}

#[test]
fn arrow_down_then_enter_runs_the_second_row() {
    let mut h = Harness::settled();
    let ranked = h.ranked();
    assert_eq!(h.palette.highlighted(), 0);
    assert_eq!(h.frame(vec![key(Key::ArrowDown)]), PaletteOutcome::Open);
    assert_eq!(h.palette.highlighted(), 1);
    assert_eq!(
        h.frame(vec![key(Key::Enter)]),
        PaletteOutcome::Run(ranked[1])
    );

    // The same within one frame, after filtering.
    let mut h = Harness::settled();
    h.frame(vec![Event::Text("e".into())]);
    let ranked = h.ranked();
    assert_eq!(
        h.frame(vec![key(Key::ArrowDown), key(Key::Enter)]),
        PaletteOutcome::Run(ranked[1])
    );
}

#[test]
fn up_and_down_stop_at_the_ends_and_typing_resets_the_highlight() {
    let mut h = Harness::settled();
    let last = h.ranked().len() - 1;
    h.frame(vec![key(Key::ArrowUp)]);
    assert_eq!(h.palette.highlighted(), 0);
    h.frame(vec![key(Key::ArrowDown); 3]);
    assert_eq!(h.palette.highlighted(), 3);
    h.frame(vec![key(Key::PageDown); 3]);
    assert_eq!(h.palette.highlighted(), last);
    h.frame(vec![key(Key::ArrowUp)]);
    assert_eq!(h.palette.highlighted(), last - 1);
    h.frame(vec![key(Key::PageUp)]);
    assert_eq!(h.palette.highlighted(), last - 1 - VISIBLE_ROWS);
    h.frame(vec![Event::Text("s".into())]);
    assert_eq!(h.palette.highlighted(), 0);
}

#[test]
fn the_highlighted_row_scrolls_into_view() {
    let mut h = Harness::settled();
    let ranked = h.ranked();
    let target = VISIBLE_ROWS + 2;
    h.frame(vec![key(Key::ArrowDown); target]);
    h.frame(vec![]);
    let rect = h
        .palette
        .row_rect(ranked[target])
        .expect("scrolled into view");
    assert_eq!(rect.height(), h.theme.toolbar_button);
    // And back up to the top.
    h.frame(vec![key(Key::PageUp); 2]);
    h.frame(vec![]);
    let rect = h.palette.row_rect(ranked[0]).expect("first row visible");
    assert_eq!(rect.height(), h.theme.toolbar_button);
}

#[test]
fn enter_on_a_disabled_row_does_not_run() {
    let mut h = Harness::settled();
    // The disabled entry is listed last.
    let ranked = h.ranked();
    assert_eq!(*ranked.last().unwrap(), EXIT_SKETCH);
    h.frame(vec![key(Key::ArrowDown); ranked.len()]);
    assert_eq!(h.palette.highlighted(), ranked.len() - 1);
    assert_eq!(h.frame(vec![key(Key::Enter)]), PaletteOutcome::Open);

    // Also when it is the only match.
    let mut h = Harness::settled();
    h.frame(vec![Event::Text("exit".into())]);
    assert_eq!(h.ranked(), vec![EXIT_SKETCH]);
    assert_eq!(h.frame(vec![key(Key::Enter)]), PaletteOutcome::Open);
}

#[test]
fn enter_with_no_match_does_nothing() {
    let mut h = Harness::settled();
    h.frame(vec![Event::Text("zzzz".into())]);
    assert!(h.ranked().is_empty());
    assert_eq!(h.frame(vec![key(Key::Enter)]), PaletteOutcome::Open);
    assert!(h
        .entries
        .iter()
        .enumerate()
        .all(|(i, _)| h.palette.row_rect(i).is_none()));
}

#[test]
fn escape_closes() {
    let mut h = Harness::settled();
    assert_eq!(h.frame(vec![key(Key::Escape)]), PaletteOutcome::Closed);
    // Even with text typed.
    let mut h = Harness::settled();
    h.frame(vec![Event::Text("ex".into())]);
    assert_eq!(h.frame(vec![key(Key::Escape)]), PaletteOutcome::Closed);
}

#[test]
fn clicking_an_enabled_row_runs_it() {
    let mut h = Harness::settled();
    let ranked = h.ranked();
    let entry = ranked[2];
    let rect = h.palette.row_rect(entry).expect("row shown");
    let outcomes = h.click(rect.center());
    assert_eq!(
        outcomes,
        [
            PaletteOutcome::Open,
            PaletteOutcome::Open,
            PaletteOutcome::Run(entry)
        ]
    );
}

#[test]
fn moving_the_pointer_over_a_row_highlights_it() {
    let mut h = Harness::settled();
    let ranked = h.ranked();
    let rect = h.palette.row_rect(ranked[4]).unwrap();
    h.frame(vec![Event::PointerMoved(rect.center())]);
    h.frame(vec![Event::PointerMoved(rect.center() + vec2(5.0, 0.0))]);
    assert_eq!(h.palette.highlighted(), 4);
}

#[test]
fn clicking_a_disabled_row_does_not_run_it() {
    let mut h = Harness::settled();
    h.frame(vec![Event::Text("sketch".into())]);
    h.frame(vec![]);
    let rect = h
        .palette
        .row_rect(EXIT_SKETCH)
        .expect("disabled row listed");
    assert_eq!(h.click(rect.center()), [PaletteOutcome::Open; 3]);
}

#[test]
fn a_click_outside_the_panel_closes_and_inside_does_not() {
    let mut h = Harness::settled();
    let panel = h.palette.panel_rect().unwrap();
    // On the query field: stays open.
    let field = pos2(panel.center().x, panel.top() + 10.0);
    assert_eq!(h.click(field), [PaletteOutcome::Open; 3]);
    // Outside: closes on the press.
    let outside = pos2(20.0, SCREEN.y - 20.0);
    assert!(!panel.contains(outside));
    let [moved, pressed, _] = h.click(outside);
    assert_eq!(moved, PaletteOutcome::Open);
    assert_eq!(pressed, PaletteOutcome::Closed);
}

#[test]
fn the_press_that_opened_the_palette_does_not_close_it() {
    let mut h = Harness::new();
    let pos = pos2(20.0, 20.0);
    let outcome = h.frame(vec![
        Event::PointerMoved(pos),
        Event::PointerButton {
            pos,
            button: PointerButton::Primary,
            pressed: true,
            modifiers: Modifiers::NONE,
        },
    ]);
    assert_eq!(outcome, PaletteOutcome::Open);
}

#[test]
fn the_query_field_keeps_the_keyboard_focus() {
    let mut h = Harness::settled();
    // Tab would normally move the focus on; the palette keeps it.
    h.frame(vec![key(Key::Tab)]);
    h.frame(vec![Event::Text("box".into())]);
    assert_eq!(h.palette.query, "box");
    // Someone else takes the focus; it comes back the next frame.
    h.ctx
        .memory_mut(|m| m.request_focus(egui::Id::new("elsewhere")));
    h.frame(vec![]);
    h.frame(vec![Event::Text("x".into())]);
    assert_eq!(h.palette.query, "boxx");
}

#[test]
fn the_search_field_senses_clicks() {
    let ctx = egui::Context::default();
    let theme = Theme::light();
    let rect = std::cell::Cell::new(Rect::NOTHING);
    let clicked = std::cell::Cell::new(false);
    let run = |events: Vec<Event>| {
        let input = egui::RawInput {
            events,
            ..Default::default()
        };
        let mut out = ctx.run_ui(input, |ui| {
            let r = search_field(ui, &theme, "Search commands", "Ctrl+K", 260.0);
            rect.set(r.rect);
            clicked.set(clicked.get() | r.clicked());
        });
        out.textures_delta.clear();
    };
    run(vec![]);
    let pos = rect.get().center();
    assert_eq!(rect.get().size(), vec2(260.0, theme.interact_height));
    run(vec![
        Event::PointerMoved(pos),
        Event::PointerButton {
            pos,
            button: PointerButton::Primary,
            pressed: true,
            modifiers: Modifiers::NONE,
        },
    ]);
    run(vec![Event::PointerButton {
        pos,
        button: PointerButton::Primary,
        pressed: false,
        modifiers: Modifiers::NONE,
    }]);
    assert!(clicked.get());
}

// ---------------------------------------------------------------------------
// Rendering (GPU)

#[test]
fn the_palette_draws_inside_its_panel_and_shadow_only() {
    let theme = Theme::light();
    let entries = commands();
    let mut snap = snapshot([640, 480], 1.0);
    theme.apply(snap.context());
    let mut palette = CommandPalette::new();
    let mut image = None;
    // Frame 1 sizes the area (egui draws areas from their second frame).
    for _ in 0..4 {
        image = Some(
            snap.frame(egui::RawInput::default(), |ctx| {
                ctx.layer_painter(LayerId::background()).rect_filled(
                    ctx.content_rect(),
                    0.0,
                    theme.chrome,
                );
                let outcome = palette.show(ctx, &theme, pos2(320.0, 30.0), 400.0, &entries);
                assert_eq!(outcome, PaletteOutcome::Open);
            })
            .unwrap(),
        );
    }
    let image = image.unwrap();
    maybe_dump("palette", &image);

    let panel = palette.panel_rect().unwrap();
    assert_eq!(panel.width(), 400.0);
    let shadow = panel + theme.visuals().window_shadow.margin();
    let reach = shadow.expand(1.0);
    let (mut ink, mut fill, mut area) = (0, 0, 0);
    for y in 0..image.height {
        for x in 0..image.width {
            let p = image.pixel(x, y);
            let centre = pos2(x as f32 + 0.5, y as f32 + 0.5);
            if panel.shrink(1.0).contains(centre) {
                area += 1;
                if p == rgba(theme.panel) {
                    fill += 1;
                } else if max_channel_diff(p, rgba(theme.panel)) > 64 {
                    ink += 1;
                }
            } else if !reach.contains(centre) {
                assert!(
                    max_channel_diff(p, rgba(theme.chrome)) <= 1,
                    "paint at ({x}, {y}) outside the panel and its shadow: {p:?}"
                );
            }
        }
    }
    // Mostly the panel's white, with text, icons and the highlight on it.
    assert!(
        fill * 2 > area,
        "panel fill covers only {fill} of {area} px"
    );
    assert!(ink > 2000, "only {ink} inked pixels in the panel");
    // The highlighted first row has the accent-soft background.
    let first = palette.row_rect(palette.ranked(&entries)[0]).unwrap();
    let probe = pos2(first.right() - 40.0, first.top() + 3.0);
    assert!(
        max_channel_diff(
            image.pixel(probe.x as u32, probe.y as u32),
            rgba(theme.accent_soft)
        ) <= 1
    );
}

#[test]
fn the_search_field_border_darkens_on_hover() {
    let theme = Theme::light();
    let mut snap = snapshot([320, 40], 1.0);
    theme.apply(snap.context());
    let rect = std::cell::Cell::new(Rect::NOTHING);
    let mut run = |ui: &mut egui::Ui| {
        ui.painter().rect_filled(ui.max_rect(), 0.0, theme.chrome);
        let r = search_field(ui, &theme, "Search commands", "Ctrl+K", 280.0);
        rect.set(r.rect);
    };
    let rest = snap.frame_ui(egui::RawInput::default(), &mut run).unwrap();
    let r = rect.get();
    let mut input = egui::RawInput::default();
    input.events.push(Event::PointerMoved(r.center()));
    snap.frame_ui(input, &mut run).unwrap();
    let hovered = snap.frame_ui(egui::RawInput::default(), &mut run).unwrap();
    maybe_dump("search_field_hovered", &hovered);

    // The middle of the top edge: the 1 px border.
    let (x, y) = (r.center().x as u32, r.top() as u32);
    assert!(max_channel_diff(rest.pixel(x, y), rgba(theme.border)) <= 1);
    assert!(max_channel_diff(hovered.pixel(x, y), rgba(theme.border_strong)) <= 1);
    // Inside: the panel fill, with the hint text and icon drawn on it.
    assert!(max_channel_diff(rest.pixel(x, y + 2), rgba(theme.panel)) <= 1);
    let inked = (r.left() as u32..r.right() as u32)
        .flat_map(|x| (r.top() as u32 + 2..r.bottom() as u32 - 2).map(move |y| (x, y)))
        .filter(|&(x, y)| max_channel_diff(rest.pixel(x, y), rgba(theme.panel)) > 64)
        .count();
    assert!(inked > 100, "only {inked} inked pixels in the field");
}
