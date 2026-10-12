//! Lays out a line of text in one of OpenDatum's bundled fonts as closed
//! outline contours, for sketch text (ADR 0015).
//!
//! The layout is in **box units**: the text's box is 1 high (the font's line
//! height, ascender to descender) and [`Layout::aspect`] wide (the shaped
//! advance width), with the baseline at the descender's depth above the
//! bottom. Sketches store the result and map it into their text box, so a
//! document's letters never depend on this crate's version.
//!
//! Shaping (kerning, ligatures, contextual forms) is harfrust, the HarfBuzz
//! project's port of its shaping algorithm to Rust; line metrics and glyph
//! outlines are read with skrifa (unhinted, in font units). Outlines come
//! from the font's
//! `glyf` table (quadratic, OpenType specification §"glyf"; Apple TrueType
//! Reference Manual ch. 1) or `CFF` table (cubic, Adobe Technical Note
//! #5177). Runs of quadratic segments joined at implied on-curve points (the
//! midpoint of two successive off-curve points) are exactly a clamped
//! quadratic B-spline with uniform interior knots whose control points are
//! the run's start, its off-curve points and its end (Piegl & Tiller, *The
//! NURBS Book*, §5.6 knot points), so each run becomes one piece.

use glam::DVec2;
use harfrust::{BufferFlags, ShapeOptions, ShaperFont};
use skrifa::instance::{LocationRef, Size};
use skrifa::outline::{DrawSettings, OutlinePen};
use skrifa::{FontRef, GlyphId, MetadataProvider};
use thiserror::Error;

/// The longest text accepted, in characters.
pub const MAX_CHARS: usize = 1000;

/// A bundled font face (ADR 0015): the Liberation fonts 2.1.5, under the
/// SIL Open Font Licence 1.1 (`fonts/OFL.txt`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Font {
    /// Liberation Sans.
    Sans,
    /// Liberation Sans Bold.
    SansBold,
    /// Liberation Serif.
    Serif,
    /// Liberation Serif Bold.
    SerifBold,
    /// Liberation Mono.
    Mono,
    /// Liberation Mono Bold.
    MonoBold,
}

impl Font {
    /// Every bundled face, in menu order.
    pub const ALL: [Font; 6] = [
        Font::Sans,
        Font::SansBold,
        Font::Serif,
        Font::SerifBold,
        Font::Mono,
        Font::MonoBold,
    ];

    /// The face new text uses.
    pub const DEFAULT: Font = Font::Sans;

    /// The face's name, as documents store it.
    pub fn name(self) -> &'static str {
        match self {
            Font::Sans => "Liberation Sans",
            Font::SansBold => "Liberation Sans Bold",
            Font::Serif => "Liberation Serif",
            Font::SerifBold => "Liberation Serif Bold",
            Font::Mono => "Liberation Mono",
            Font::MonoBold => "Liberation Mono Bold",
        }
    }

    /// The bundled face named `name`, if there is one.
    pub fn from_name(name: &str) -> Option<Font> {
        Font::ALL.into_iter().find(|f| f.name() == name)
    }

    /// The font file.
    pub fn data(self) -> &'static [u8] {
        match self {
            Font::Sans => include_bytes!("../fonts/LiberationSans-Regular.ttf"),
            Font::SansBold => include_bytes!("../fonts/LiberationSans-Bold.ttf"),
            Font::Serif => include_bytes!("../fonts/LiberationSerif-Regular.ttf"),
            Font::SerifBold => include_bytes!("../fonts/LiberationSerif-Bold.ttf"),
            Font::Mono => include_bytes!("../fonts/LiberationMono-Regular.ttf"),
            Font::MonoBold => include_bytes!("../fonts/LiberationMono-Bold.ttf"),
        }
    }
}

/// One piece of a contour, in box units: a line (degree 1, two points) or
/// a clamped B-spline of degree 2 or 3 with uniform interior knots on
/// [0, 1] (`points` are its control points, at least `degree + 1`).
#[derive(Clone, Debug, PartialEq)]
pub struct Piece {
    /// 1, 2 or 3.
    pub degree: u32,
    /// The line's ends or the spline's control points; the first and last
    /// lie on the outline.
    pub points: Vec<DVec2>,
}

impl Piece {
    /// Where the piece starts.
    pub fn start(&self) -> DVec2 {
        self.points[0]
    }

    /// Where the piece ends.
    pub fn end(&self) -> DVec2 {
        self.points[self.points.len() - 1]
    }
}

/// A laid-out line of text.
#[derive(Clone, Debug, PartialEq)]
pub struct Layout {
    /// Box width over box height (> 0).
    pub aspect: f64,
    /// Closed contours, glyph by glyph in visual order. Each piece starts
    /// exactly where the previous one ends, and the last ends exactly where
    /// the first starts. Outer contours run clockwise and holes
    /// counter-clockwise (the TrueType convention); contours of different
    /// glyphs may touch or overlap where the font's spacing makes them.
    pub contours: Vec<Vec<Piece>>,
}

/// Why a text cannot be laid out.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum LayoutError {
    /// Nothing to lay out.
    #[error("the text is empty")]
    Empty,
    /// More than [`MAX_CHARS`] characters.
    #[error("the text is longer than {MAX_CHARS} characters")]
    TooLong,
    /// Sketch text is a single line.
    #[error("sketch text is a single line: it cannot contain a line break")]
    LineBreak,
    /// A tab or other control character.
    #[error("the text cannot contain control characters such as tabs")]
    Control,
    /// The font has no glyph for a character.
    #[error("{font} has no character {ch:?} (U+{code:04X})")]
    Missing {
        /// The font's name.
        font: &'static str,
        /// The character.
        ch: char,
        /// Its code point.
        code: u32,
    },
    /// The bundled font file could not be read (never expected).
    #[error("{0} could not be read")]
    Font(&'static str),
}

/// Lays out `text` in `font` (see the crate docs for the units).
pub fn lay_out(text: &str, font: Font) -> Result<Layout, LayoutError> {
    if text.is_empty() {
        return Err(LayoutError::Empty);
    }
    if text.chars().count() > MAX_CHARS {
        return Err(LayoutError::TooLong);
    }
    for c in text.chars() {
        if matches!(
            c,
            '\n' | '\r' | '\u{0B}' | '\u{0C}' | '\u{85}' | '\u{2028}' | '\u{2029}'
        ) {
            return Err(LayoutError::LineBreak);
        }
        if c.is_control() {
            return Err(LayoutError::Control);
        }
    }
    let unreadable = || LayoutError::Font(font.name());
    let face = FontRef::new(font.data()).map_err(|_| unreadable())?;
    let metrics = face.metrics(Size::unscaled(), LocationRef::default());
    let ascender = f64::from(metrics.ascent);
    let descender = f64::from(metrics.descent);
    let height = ascender - descender;
    if !height.is_finite() || height <= 0.0 {
        return Err(unreadable());
    }
    let shaping = harfrust::Font::new(font.data(), 0).ok_or_else(unreadable)?;
    let mut buffer = harfrust::Buffer::new();
    buffer.push_str(text);
    buffer.guess_segment_properties();
    // The buffer is the whole text, so its ends are the text's ends.
    buffer.set_flags(BufferFlags::BEGINNING_OF_TEXT | BufferFlags::END_OF_TEXT);
    harfrust::shape(&ShaperFont::new(&shaping), &mut buffer, ShapeOptions::new())
        .map_err(|_| unreadable())?;
    let glyphs = face.outline_glyphs();
    let mut pen = 0i64;
    let mut contours = Vec::new();
    for (info, pos) in buffer.glyph_infos().iter().zip(buffer.glyph_positions()) {
        if info.glyph_id == 0 {
            let ch = text
                .get(info.cluster as usize..)
                .and_then(|s| s.chars().next())
                .unwrap_or('\u{FFFD}');
            return Err(LayoutError::Missing {
                font: font.name(),
                ch,
                code: u32::from(ch),
            });
        }
        let glyph = glyphs
            .get(GlyphId::new(info.glyph_id))
            .ok_or_else(unreadable)?;
        let mut outline = Outline {
            dx: (pen + i64::from(pos.x_offset)) as f64,
            dy: f64::from(pos.y_offset) - descender,
            height,
            contours: &mut contours,
            segs: Vec::new(),
            start: DVec2::ZERO,
            last: DVec2::ZERO,
        };
        let unhinted = DrawSettings::unhinted(Size::unscaled(), LocationRef::default());
        glyph
            .draw(unhinted, &mut outline)
            .map_err(|_| unreadable())?;
        outline.finish();
        pen += i64::from(pos.x_advance);
    }
    let aspect = pen as f64 / height;
    if !(aspect > 0.0 && aspect.is_finite()) {
        // Only zero-width marks: nothing a box can hold.
        return Err(LayoutError::Empty);
    }
    Ok(Layout { aspect, contours })
}

/// A segment in font units, before conversion.
#[derive(Clone, Copy, Debug)]
enum Seg {
    Line(DVec2, DVec2),
    Quad(DVec2, DVec2, DVec2),
    Cubic(DVec2, DVec2, DVec2, DVec2),
}

/// Collects one glyph's contours, offset to its pen position.
struct Outline<'a> {
    /// Horizontal offset (font units).
    dx: f64,
    /// Vertical offset from the baseline to the box bottom (font units).
    dy: f64,
    /// Box height (font units).
    height: f64,
    contours: &'a mut Vec<Vec<Piece>>,
    segs: Vec<Seg>,
    start: DVec2,
    last: DVec2,
}

impl Outline<'_> {
    fn p(&self, x: f32, y: f32) -> DVec2 {
        DVec2::new(f64::from(x), f64::from(y))
    }

    /// Box units of a font-unit point. Equal inputs give equal outputs, so
    /// pieces that share an end in font units share it exactly here.
    fn to_box(&self, p: DVec2) -> DVec2 {
        DVec2::new((p.x + self.dx) / self.height, (p.y + self.dy) / self.height)
    }

    fn finish(&mut self) {
        if !self.segs.is_empty() {
            self.close();
        }
    }
}

impl OutlinePen for Outline<'_> {
    fn move_to(&mut self, x: f32, y: f32) {
        self.finish();
        self.start = self.p(x, y);
        self.last = self.start;
    }

    fn line_to(&mut self, x: f32, y: f32) {
        let p = self.p(x, y);
        self.segs.push(Seg::Line(self.last, p));
        self.last = p;
    }

    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        let (c, p) = (self.p(x1, y1), self.p(x, y));
        self.segs.push(Seg::Quad(self.last, c, p));
        self.last = p;
    }

    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        let (c1, c2, p) = (self.p(x1, y1), self.p(x2, y2), self.p(x, y));
        self.segs.push(Seg::Cubic(self.last, c1, c2, p));
        self.last = p;
    }

    fn close(&mut self) {
        if self.segs.is_empty() {
            return;
        }
        if self.last != self.start {
            self.segs.push(Seg::Line(self.last, self.start));
        }
        let segs = std::mem::take(&mut self.segs);
        self.last = self.start;
        let pieces: Vec<Piece> = split_loops(pieces(&segs))
            .into_iter()
            .map(|(degree, pts)| Piece {
                degree,
                points: pts.into_iter().map(|p| self.to_box(p)).collect(),
            })
            .collect();
        if pieces.len() >= 2 && contour_area(&pieces).abs() > 1e-9 {
            self.contours.push(pieces);
        }
    }
}

/// Turns one closed contour's segments (font units) into pieces: drops
/// zero-length segments, turns segments whose control points sit on their
/// ends into lines, and joins runs of quadratics meeting at implied
/// on-curve points into one quadratic B-spline each.
fn pieces(segs: &[Seg]) -> Vec<(u32, Vec<DVec2>)> {
    let mut out: Vec<(u32, Vec<DVec2>)> = Vec::new();
    // The quadratic run being built: its start, then its off-curve points,
    // and where it currently ends.
    let mut run: Option<(Vec<DVec2>, DVec2)> = None;
    let flush = |run: &mut Option<(Vec<DVec2>, DVec2)>, out: &mut Vec<(u32, Vec<DVec2>)>| {
        if let Some((mut pts, end)) = run.take() {
            pts.push(end);
            out.push((2, pts));
        }
    };
    for seg in segs {
        match *seg {
            Seg::Line(a, b) => {
                if a == b {
                    continue;
                }
                flush(&mut run, &mut out);
                out.push((1, vec![a, b]));
            }
            Seg::Quad(a, c, b) => {
                if a == b {
                    // Out and back along the same path: no area, no gap.
                    continue;
                }
                if c == a || c == b {
                    flush(&mut run, &mut out);
                    out.push((1, vec![a, b]));
                    continue;
                }
                match &mut run {
                    Some((pts, end)) => {
                        let prev = pts[pts.len() - 1];
                        // Joined at the implied midpoint of two distinct
                        // off-curve points: the same B-spline continues.
                        if *end == a && prev != c && (prev + c) * 0.5 == a {
                            pts.push(c);
                            *end = b;
                        } else {
                            flush(&mut run, &mut out);
                            run = Some((vec![a, c], b));
                        }
                    }
                    None => run = Some((vec![a, c], b)),
                }
            }
            Seg::Cubic(a, c1, c2, b) => {
                flush(&mut run, &mut out);
                if a == b && c1 == a && c2 == a {
                    continue;
                }
                if c1 == a && c2 == b {
                    out.push((1, vec![a, b]));
                } else {
                    out.push((3, vec![a, c1, c2, b]));
                }
            }
        }
    }
    flush(&mut run, &mut out);
    out
}

/// Splits every piece that ends where it starts (a contour drawn as one
/// run of off-curve points) in two, so no piece is closed on itself.
fn split_loops(pieces: Vec<(u32, Vec<DVec2>)>) -> Vec<(u32, Vec<DVec2>)> {
    let mut out = Vec::with_capacity(pieces.len() + 1);
    for (degree, pts) in pieces {
        if pts[0] != pts[pts.len() - 1] {
            out.push((degree, pts));
            continue;
        }
        match degree {
            // [a, c1 … cn, a] splits at the implied point between c_k and
            // c_k+1 into two runs of the same B-spline.
            2 if pts.len() >= 4 => {
                let k = (pts.len() - 1) / 2;
                let m = (pts[k] + pts[k + 1]) * 0.5;
                let mut first = pts[..=k].to_vec();
                first.push(m);
                let mut second = vec![m];
                second.extend_from_slice(&pts[k + 1..]);
                out.push((2, first));
                out.push((2, second));
            }
            // A closed cubic splits at its middle (de Casteljau); both
            // halves share the computed point exactly.
            3 => {
                let (a, b, c, d) = (pts[0], pts[1], pts[2], pts[3]);
                let (ab, bc, cd) = ((a + b) * 0.5, (b + c) * 0.5, (c + d) * 0.5);
                let (abc, bcd) = ((ab + bc) * 0.5, (bc + cd) * 0.5);
                let m = (abc + bcd) * 0.5;
                out.push((3, vec![a, ab, abc, m]));
                out.push((3, vec![m, bcd, cd, d]));
            }
            // A closed line or single quadratic has no area.
            _ => {}
        }
    }
    out
}

/// Signed area of a closed contour from its control polygons (positive
/// counter-clockwise). Only the sign and whether it vanishes are used: a
/// contour with no area is dropped.
fn contour_area(pieces: &[Piece]) -> f64 {
    let mut twice = 0.0;
    for piece in pieces {
        for w in piece.points.windows(2) {
            twice += w[0].perp_dot(w[1]);
        }
    }
    twice * 0.5
}
