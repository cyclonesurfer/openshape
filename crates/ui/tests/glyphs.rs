//! Every non-ASCII character in the workspace's source code outside comments
//! (labels, constraint glyphs, units, error messages the app shows) exists in
//! the fonts the app installs, so none of it renders as a missing-glyph box.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

/// Non-ASCII characters outside `//` comments, with one place each occurs.
fn source_characters() -> BTreeMap<char, String> {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut files = Vec::new();
    for entry in std::fs::read_dir(&crates).unwrap() {
        let src = entry.unwrap().path().join("src");
        if src.is_dir() {
            rust_files(&src, &mut files);
        }
    }
    assert!(files.len() > 20, "found only {} source files", files.len());
    let mut chars = BTreeMap::new();
    for file in files {
        let text = std::fs::read_to_string(&file).unwrap();
        for (n, line) in text.lines().enumerate() {
            // Dropping everything after `//` also drops string text after a
            // `//` inside a literal; that only makes the scan miss a
            // character, never report a false one.
            let code = line.split("//").next().unwrap_or("");
            for c in code.chars().filter(|c| !c.is_ascii()) {
                chars
                    .entry(c)
                    .or_insert_with(|| format!("{}:{}", file.display(), n + 1));
            }
        }
    }
    chars
}

#[test]
fn every_source_character_has_a_glyph() {
    let ctx = egui::Context::default();
    ctx.set_fonts(openshape_ui_kit::font_definitions());
    let mut out = ctx.run_ui(egui::RawInput::default(), |_| {});
    out.textures_delta.clear();
    // `Fonts::has_glyph` reports false for any character that lives in the
    // same face as egui's replacement glyph, so read the character map.
    let available = ctx.fonts_mut(|f| {
        f.fonts
            .font(&egui::FontFamily::Proportional)
            .characters()
            .clone()
    });
    let chars = source_characters();
    assert!(
        chars.contains_key(&'°'),
        "the scan should see the degree sign"
    );
    let missing: Vec<String> = chars
        .iter()
        .filter(|(c, _)| !available.contains_key(c))
        .map(|(c, at)| format!("{c:?} (U+{:04X}) at {at}", *c as u32))
        .collect();
    assert!(
        missing.is_empty(),
        "characters with no glyph in the app's proportional fonts:\n{}",
        missing.join("\n")
    );
}
