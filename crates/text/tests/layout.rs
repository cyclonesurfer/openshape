//! Text layout: closed, exact contours that reproduce the font's outlines,
//! shaped widths, and refusals (ADR 0015).

use glam::DVec2;
use openshape_text::{lay_out, Font, Layout, LayoutError, Piece};
use skrifa::instance::{LocationRef, Size};
use skrifa::outline::{DrawSettings, OutlinePen};
use skrifa::raw::tables::os2::SelectionFlags;
use skrifa::raw::TableProvider;
use skrifa::{FontRef, MetadataProvider};

/// Point of a clamped B-spline with uniform interior knots on [0, 1]
/// (de Boor's algorithm, Piegl & Tiller A3.1), the convention of sketch
/// control splines.
fn de_boor(degree: usize, ctrl: &[DVec2], u: f64) -> DVec2 {
    let n = ctrl.len() - 1;
    let spans = n + 1 - degree;
    let mut knots = vec![0.0; degree + 1];
    for i in 1..spans {
        knots.push(i as f64 / spans as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    let mut k = degree;
    while k < n && u >= knots[k + 1] {
        k += 1;
    }
    let mut d: Vec<DVec2> = (0..=degree).map(|j| ctrl[j + k - degree]).collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = j + k - degree;
            let den = knots[i + degree + 1 - r] - knots[i];
            let a = if den == 0.0 {
                0.0
            } else {
                (u - knots[i]) / den
            };
            d[j] = d[j - 1] * (1.0 - a) + d[j] * a;
        }
    }
    d[degree]
}

fn piece_at(p: &Piece, u: f64) -> DVec2 {
    de_boor(p.degree as usize, &p.points, u)
}

/// Every contour closes exactly and its pieces join exactly.
fn assert_closed(layout: &Layout) {
    for (i, c) in layout.contours.iter().enumerate() {
        assert!(c.len() >= 2, "contour {i} has {} pieces", c.len());
        for (j, p) in c.iter().enumerate() {
            let next = &c[(j + 1) % c.len()];
            assert_eq!(p.end(), next.start(), "contour {i} piece {j} does not join");
            assert!((1..=3).contains(&p.degree));
            assert!(p.points.len() > p.degree as usize);
            assert_ne!(
                p.start(),
                p.end(),
                "contour {i} piece {j} is closed on itself"
            );
            for w in p.points.windows(2) {
                assert_ne!(w[0], w[1], "contour {i} piece {j} repeats a point");
            }
            for q in &p.points {
                assert!(q.is_finite());
            }
        }
    }
}

/// Signed area by Green's theorem over the pieces, sampled finely.
fn area(c: &[Piece]) -> f64 {
    let mut twice = 0.0;
    for p in c {
        let n = 64;
        let mut prev = p.start();
        for k in 1..=n {
            let q = piece_at(p, k as f64 / n as f64);
            twice += prev.perp_dot(q);
            prev = q;
        }
    }
    twice * 0.5
}

/// The font's line, read straight from its tables: the box height
/// (ascender − descender, font units) and the descender. The Liberation
/// fonts do not ask for their typographic metrics, so the line is hhea's.
fn line(font: Font) -> (f64, f64) {
    let face = FontRef::new(font.data()).unwrap();
    let os2 = face.os2().unwrap();
    assert!(!os2
        .fs_selection()
        .contains(SelectionFlags::USE_TYPO_METRICS));
    let hhea = face.hhea().unwrap();
    let ascender = f64::from(hhea.ascender().to_i16());
    let descender = f64::from(hhea.descender().to_i16());
    (ascender - descender, descender)
}

/// One glyph's outline straight from the font (font units, unhinted).
fn raw_outline(font: Font, ch: char) -> Raw {
    let face = FontRef::new(font.data()).unwrap();
    let gid = face.charmap().map(ch).unwrap();
    let mut raw = Raw::default();
    face.outline_glyphs()
        .get(gid)
        .unwrap()
        .draw(
            DrawSettings::unhinted(Size::unscaled(), LocationRef::default()),
            &mut raw,
        )
        .unwrap();
    raw
}

/// The glyph outlines of a string straight from the font (font units),
/// each segment as its Bézier control points.
/// Contours are closed (OpenType `glyf`): a pen that is not brought back
/// to the start before `close` gets the closing line.
#[derive(Default)]
struct Raw {
    segs: Vec<Vec<DVec2>>,
    start: DVec2,
    last: DVec2,
}

impl OutlinePen for Raw {
    fn move_to(&mut self, x: f32, y: f32) {
        self.last = DVec2::new(x.into(), y.into());
        self.start = self.last;
    }
    fn line_to(&mut self, x: f32, y: f32) {
        let p = DVec2::new(x.into(), y.into());
        self.segs.push(vec![self.last, p]);
        self.last = p;
    }
    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        let p = DVec2::new(x.into(), y.into());
        self.segs
            .push(vec![self.last, DVec2::new(x1.into(), y1.into()), p]);
        self.last = p;
    }
    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        let p = DVec2::new(x.into(), y.into());
        self.segs.push(vec![
            self.last,
            DVec2::new(x1.into(), y1.into()),
            DVec2::new(x2.into(), y2.into()),
            p,
        ]);
        self.last = p;
    }
    fn close(&mut self) {
        if self.last != self.start {
            self.segs.push(vec![self.last, self.start]);
            self.last = self.start;
        }
    }
}

fn bezier(ctrl: &[DVec2], t: f64) -> DVec2 {
    let mut d = ctrl.to_vec();
    for r in 1..d.len() {
        for j in 0..d.len() - r {
            d[j] = d[j] * (1.0 - t) + d[j + 1] * t;
        }
    }
    d[0]
}

#[test]
fn every_face_lays_out_closed_exact_contours() {
    for font in Font::ALL {
        let layout = lay_out("OpenShape 2026 & Co. — ÅØ gjy%@", font).unwrap();
        assert!(layout.aspect > 5.0 && layout.aspect < 30.0, "{font:?}");
        assert_closed(&layout);
        // Letters stay near their box (accents may rise a little above it).
        for c in &layout.contours {
            for p in c {
                for q in &p.points {
                    assert!(q.x > -0.1 && q.x < layout.aspect + 0.1, "{font:?} {q}");
                    assert!(q.y > -0.1 && q.y < 1.25, "{font:?} {q}");
                }
            }
        }
    }
}

#[test]
fn pieces_reproduce_the_font_outline() {
    // Each piece, sampled, lies on the font's own segments: the joined
    // quadratic runs are the same curve as the segments they replace.
    let font = Font::Sans;
    let (height, descender) = line(font);
    for ch in ['O', 'S', 'g', '8', '@', 'e', 'Q', 'B', '&'] {
        let layout = lay_out(&ch.to_string(), font).unwrap();
        let raw = raw_outline(font, ch);
        let raw_box: Vec<Vec<DVec2>> = raw
            .segs
            .iter()
            .map(|s| {
                s.iter()
                    .map(|p| DVec2::new(p.x / height, (p.y - descender) / height))
                    .collect()
            })
            .collect();
        let raw_pts: Vec<DVec2> = raw_box
            .iter()
            .flat_map(|s| (0..=200).map(move |k| bezier(s, k as f64 / 200.0)))
            .collect();
        let mut worst = 0.0f64;
        let mut count = 0;
        for c in &layout.contours {
            for p in c {
                for k in 0..=50 {
                    let q = piece_at(p, k as f64 / 50.0);
                    let d = raw_pts
                        .iter()
                        .map(|r| r.distance(q))
                        .fold(f64::INFINITY, f64::min);
                    worst = worst.max(d);
                    count += 1;
                }
            }
        }
        assert!(count > 0);
        // Raw samples are 1/200 apart along each segment: a point on the
        // outline is within half a sample spacing of one.
        assert!(worst < 2e-3, "{ch}: a piece leaves the outline by {worst}");
        // And the pieces are fewer than the font's segments for curvy
        // letters: runs were joined.
        let pieces: usize = layout.contours.iter().map(Vec::len).sum();
        assert!(pieces <= raw.segs.len(), "{ch}");
    }
}

#[test]
fn joined_runs_meet_the_font_segments_at_their_knots() {
    // At the knot parameters a joined run passes exactly through the font's
    // segment ends (implied on-curve points), to rounding.
    let font = Font::Serif;
    let (height, descender) = line(font);
    let layout = lay_out("o", font).unwrap();
    let raw = raw_outline(font, 'o');
    let ends: Vec<DVec2> = raw
        .segs
        .iter()
        .map(|s| {
            let p = s[s.len() - 1];
            DVec2::new(p.x / height, (p.y - descender) / height)
        })
        .collect();
    let mut checked = 0;
    for c in &layout.contours {
        for p in c.iter().filter(|p| p.degree == 2 && p.points.len() > 3) {
            let spans = p.points.len() - 2;
            for i in 1..spans {
                let q = piece_at(p, i as f64 / spans as f64);
                let d = ends
                    .iter()
                    .map(|e| e.distance(q))
                    .fold(f64::INFINITY, f64::min);
                assert!(d < 1e-12, "knot {i} of a run is {d} from every segment end");
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "'o' has joined runs");
}

#[test]
fn letters_have_their_holes() {
    let font = Font::Sans;
    for (text, outer, holes) in [
        ("O", 1, 1),
        ("B", 1, 2),
        ("8", 1, 2),
        ("l", 1, 0),
        ("i", 2, 0),
    ] {
        let layout = lay_out(text, font).unwrap();
        let areas: Vec<f64> = layout.contours.iter().map(|c| area(c)).collect();
        // TrueType: outer contours clockwise (negative), holes anticlockwise.
        let n_outer = areas.iter().filter(|a| **a < 0.0).count();
        let n_holes = areas.iter().filter(|a| **a > 0.0).count();
        assert_eq!((n_outer, n_holes), (outer, holes), "{text}: {areas:?}");
        let total: f64 = areas.iter().sum();
        assert!(total < 0.0, "{text}: ink has area");
    }
}

#[test]
fn a_monospaced_line_is_as_wide_as_its_characters() {
    // Every Liberation Mono glyph advances by the same width, so the box's
    // aspect is that width times the character count over the line height,
    // read independently from the font's hmtx and hhea tables.
    for font in [Font::Mono, Font::MonoBold] {
        let (height, _) = line(font);
        let face = FontRef::new(font.data()).unwrap();
        let m = face.charmap().map('m').unwrap();
        let advance = f64::from(face.hmtx().unwrap().advance(m).unwrap());
        for text in ["a", "abc", "Hello, world", "    x"] {
            let layout = lay_out(text, font).unwrap();
            let expected = advance * text.chars().count() as f64 / height;
            assert_eq!(layout.aspect, expected, "{font:?} {text:?}");
        }
    }
}

#[test]
fn kerning_closes_up_awkward_pairs() {
    let font = Font::Sans;
    let a = lay_out("A", font).unwrap().aspect;
    let v = lay_out("V", font).unwrap().aspect;
    let av = lay_out("AV", font).unwrap().aspect;
    assert!(av < a + v - 1e-3, "AV {av} vs A+V {}", a + v);
}

#[test]
fn spaces_widen_the_box_without_contours() {
    let font = Font::Sans;
    let x = lay_out("x", font).unwrap();
    let spaced = lay_out("x   ", font).unwrap();
    assert!(spaced.aspect > x.aspect);
    assert_eq!(spaced.contours, x.contours);
    let blank = lay_out("   ", font).unwrap();
    assert!(blank.aspect > 0.0);
    assert!(blank.contours.is_empty());
}

#[test]
fn layout_is_deterministic() {
    for font in Font::ALL {
        let a = lay_out("Determinism 123", font).unwrap();
        let b = lay_out("Determinism 123", font).unwrap();
        assert_eq!(a, b);
    }
}

#[test]
fn what_cannot_be_laid_out_is_refused_by_name() {
    let font = Font::Sans;
    assert_eq!(lay_out("", font), Err(LayoutError::Empty));
    assert_eq!(lay_out("two\nlines", font), Err(LayoutError::LineBreak));
    assert_eq!(lay_out("a\u{2028}b", font), Err(LayoutError::LineBreak));
    assert_eq!(lay_out("tab\there", font), Err(LayoutError::Control));
    assert_eq!(
        lay_out(&"x".repeat(openshape_text::MAX_CHARS + 1), font),
        Err(LayoutError::TooLong)
    );
    let err = lay_out("Part 漢", font).unwrap_err();
    assert_eq!(
        err,
        LayoutError::Missing {
            font: "Liberation Sans",
            ch: '漢',
            code: 0x6F22
        }
    );
    assert_eq!(
        err.to_string(),
        "Liberation Sans has no character '漢' (U+6F22)"
    );
}

#[test]
fn font_names_round_trip() {
    for font in Font::ALL {
        assert_eq!(Font::from_name(font.name()), Some(font));
    }
    assert_eq!(Font::from_name("Comic Sans"), None);
    assert_eq!(Font::DEFAULT, Font::Sans);
}
