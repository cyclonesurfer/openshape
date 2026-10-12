//! The bundled font files are exactly the ones fonts/README.md records.

use openshape_text::Font;
use sha2::{Digest, Sha256};

#[test]
fn bundled_fonts_match_their_recorded_hashes() {
    let readme = include_str!("../fonts/README.md");
    for (font, file) in [
        (Font::Sans, "LiberationSans-Regular.ttf"),
        (Font::SansBold, "LiberationSans-Bold.ttf"),
        (Font::Serif, "LiberationSerif-Regular.ttf"),
        (Font::SerifBold, "LiberationSerif-Bold.ttf"),
        (Font::Mono, "LiberationMono-Regular.ttf"),
        (Font::MonoBold, "LiberationMono-Bold.ttf"),
    ] {
        let hash: String = Sha256::digest(font.data())
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        let row = format!("| {file} | {hash} |");
        assert!(readme.contains(&row), "{file}: {hash} is not recorded");
    }
}

#[test]
fn the_licence_travels_with_the_fonts() {
    let licence = include_str!("../fonts/OFL.txt");
    assert!(licence.contains("SIL OPEN FONT LICENSE Version 1.1"));
    assert!(licence.contains("Copyright (c) 2012 Red Hat, Inc."));
}
