//! Text in the solver (ADR 0015): the box's four corners kept a rectangle of
//! the text's aspect ratio by linear equations, read the right way round or
//! mirrored; degrees of freedom; redundancy; refusals; outlines mapped into
//! the box exactly; mirror and pattern copies.

mod common;

use common::{pos, B};
use openshape_sketch_solver::{
    resolve_curves, solve, Constraint, DVec2, DistanceKind, Geometry, ResolvedCurve, SolveError,
    SolveOptions, SolveStatus, TextPiece,
};

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// A square contour in box units: the letter "□" of a font with aspect 1.
fn square() -> Vec<Vec<TextPiece>> {
    let (a, b, c, d) = (
        DVec2::new(0.2, 0.2),
        DVec2::new(0.8, 0.2),
        DVec2::new(0.8, 0.8),
        DVec2::new(0.2, 0.8),
    );
    let line = |p: DVec2, q: DVec2| TextPiece {
        degree: 1,
        points: vec![p, q],
    };
    vec![vec![
        line(a, b),
        TextPiece {
            degree: 2,
            points: vec![b, DVec2::new(0.9, 0.5), c],
        },
        line(c, d),
        line(d, a),
    ]]
}

/// A text on its own four corners, placed at `bl` with height `h` and turned
/// by `angle`, with the box drawn as four construction lines:
/// (text, corners, lines bottom, right, top, left).
fn text(
    b: &mut B,
    bl: DVec2,
    h: f64,
    angle: f64,
    aspect: f64,
    mirrored: bool,
    outline: Vec<Vec<TextPiece>>,
) -> (u32, [u32; 4], [u32; 4]) {
    let up = DVec2::from_angle(angle).rotate(DVec2::Y) * h;
    let across = DVec2::new(up.y, -up.x) * aspect * if mirrored { -1.0 } else { 1.0 };
    let at = [bl, bl + across, bl + across + up, bl + up];
    let corners = at.map(|p| b.point(p.x, p.y));
    let t = b.entity(Geometry::Text {
        corners,
        content: "□".into(),
        font: "Test".into(),
        aspect,
        mirrored,
        outline,
    });
    let lines = [0, 1, 2, 3].map(|i| {
        let l = b.line_between(corners[i], corners[(i + 1) % 4]);
        b.sketch.entities.last_mut().unwrap().construction = true;
        l
    });
    (t, corners, lines)
}

/// How far the box is from its equations (metres).
fn box_error(
    s: &openshape_sketch_solver::Sketch<u32>,
    c: [u32; 4],
    aspect: f64,
    mirrored: bool,
) -> f64 {
    let [bl, br, tr, tl] = c.map(|k| pos(s, k));
    let u = tl - bl;
    let r = DVec2::new(u.y, -u.x) * if mirrored { -aspect } else { aspect };
    ((br - bl) - r)
        .length()
        .max(((tr - tl) - (br - bl)).length())
}

#[test]
fn a_text_box_has_four_freedoms_and_stays_a_rectangle() {
    for mirrored in [false, true] {
        let mut b = B::new();
        let (_, c, _) = text(
            &mut b,
            DVec2::new(0.01, 0.02),
            0.005,
            0.3,
            3.5,
            mirrored,
            square(),
        );
        // Disturb the box: the solver restores it, nearest to where it was.
        let p = b.sketch.entities.iter_mut().find(|e| e.id == c[2]).unwrap();
        if let Geometry::Point { at } = &mut p.geometry {
            *at += DVec2::new(0.001, -0.0007);
        }
        let s = solve(&b.build(), &opts());
        assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
        assert_eq!(s.dof, 4);
        assert!(box_error(&s.sketch, c, 3.5, mirrored) < 1e-12);
    }
}

#[test]
fn a_horizontal_bottom_a_fixed_corner_and_a_height_fix_the_text() {
    let mut b = B::new();
    let (_, c, lines) = text(
        &mut b,
        DVec2::new(0.0, 0.0),
        0.004,
        0.05,
        2.25,
        false,
        square(),
    );
    b.c(Constraint::Horizontal { line: lines[0] });
    b.fix(c[0]);
    b.c(Constraint::Length {
        line: lines[3],
        value: 0.010,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 0);
    let [bl, br, tr, tl] = c.map(|k| pos(&s.sketch, k));
    assert!((bl - DVec2::ZERO).length() < 1e-12);
    assert!((br - DVec2::new(0.0225, 0.0)).length() < 1e-12, "{br}");
    assert!((tr - DVec2::new(0.0225, 0.010)).length() < 1e-12, "{tr}");
    assert!((tl - DVec2::new(0.0, 0.010)).length() < 1e-12, "{tl}");
}

#[test]
fn right_angles_and_parallels_of_the_box_are_redundant() {
    let mut b = B::new();
    let (_, _, lines) = text(&mut b, DVec2::ZERO, 0.004, 0.0, 2.0, false, square());
    let perp = b.c(Constraint::Perpendicular {
        a: lines[0],
        b: lines[3],
    });
    let par = b.c(Constraint::Parallel {
        a: lines[0],
        b: lines[2],
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.dof, 4);
    let mut red = s.redundant.clone();
    red.sort();
    assert_eq!(red, vec![perp, par]);
}

#[test]
fn a_width_dimension_drives_the_height() {
    // The width follows from the height through the aspect ratio, so a
    // dimension on the bottom line sizes the text too.
    let mut b = B::new();
    let (_, c, lines) = text(&mut b, DVec2::ZERO, 0.004, 0.0, 4.0, false, square());
    b.c(Constraint::Horizontal { line: lines[0] });
    b.fix(c[0]);
    b.c(Constraint::Distance {
        a: c[0],
        b: c[1],
        kind: DistanceKind::Horizontal,
        value: 0.030,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 0);
    assert!((pos(&s.sketch, c[3]) - DVec2::new(0.0, 0.0075)).length() < 1e-12);
}

#[test]
fn width_and_height_dimensions_together_conflict() {
    let mut b = B::new();
    let (_, c, lines) = text(&mut b, DVec2::ZERO, 0.004, 0.0, 4.0, false, square());
    b.c(Constraint::Horizontal { line: lines[0] });
    b.fix(c[0]);
    b.c(Constraint::Length {
        line: lines[3],
        value: 0.005,
    });
    let width = b.c(Constraint::Length {
        line: lines[0],
        value: 0.030,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Failed(SolveError::Conflicting));
    assert!(s.conflicting.contains(&width), "{:?}", s.conflicting);
}

#[test]
fn the_outline_is_mapped_into_the_box_exactly_closed() {
    let mut b = B::new();
    let bl = DVec2::new(0.003, -0.002);
    let (t, c, _) = text(&mut b, bl, 0.004, 0.7, 1.0, false, square());
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    let curves: Vec<ResolvedCurve> = resolve_curves(&s.sketch)
        .unwrap()
        .into_iter()
        .filter(|(id, _, _)| *id == t)
        .map(|(_, construction, c)| {
            assert!(!construction);
            c
        })
        .collect();
    assert_eq!(curves.len(), 4);
    let ends = |c: &ResolvedCurve| match c {
        ResolvedCurve::Line { start, end } => (*start, *end),
        ResolvedCurve::Spline {
            points,
            degree,
            knots,
        } => {
            assert_eq!(*degree, 2);
            assert_eq!(knots, &vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
            (points[0], points[points.len() - 1])
        }
        other => panic!("{other:?}"),
    };
    for i in 0..4 {
        // Bit-for-bit: the contour is closed for the kernel.
        assert_eq!(ends(&curves[i]).1, ends(&curves[(i + 1) % 4]).0);
    }
    // Box point (0.2, 0.2) of a box 1 × 1 at bl.
    let [p0, p1, _, p3] = c.map(|k| pos(&s.sketch, k));
    let expected = p0 + (p1 - p0) * 0.2 + (p3 - p0) * 0.2;
    assert!((ends(&curves[0]).0 - expected).length() < 1e-15);
}

#[test]
fn real_letters_map_into_a_turned_box() {
    // A word laid out by the text crate keeps its contours closed and its
    // proportions in any box.
    let layout = openshape_text::lay_out("Gear 42", openshape_text::Font::Sans).unwrap();
    let outline: Vec<Vec<TextPiece>> = layout
        .contours
        .iter()
        .map(|c| {
            c.iter()
                .map(|p| TextPiece {
                    degree: p.degree,
                    points: p.points.clone(),
                })
                .collect()
        })
        .collect();
    let mut b = B::new();
    let (t, c, _) = text(
        &mut b,
        DVec2::new(0.05, 0.01),
        0.008,
        -0.4,
        layout.aspect,
        false,
        outline.clone(),
    );
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    let curves: Vec<ResolvedCurve> = resolve_curves(&s.sketch)
        .unwrap()
        .into_iter()
        .filter(|(id, _, _)| *id == t)
        .map(|(_, _, c)| c)
        .collect();
    let pieces: usize = outline.iter().map(Vec::len).sum();
    assert_eq!(curves.len(), pieces);
    let ends = |c: &ResolvedCurve| match c {
        ResolvedCurve::Line { start, end } => (*start, *end),
        ResolvedCurve::Spline { points, .. } => (points[0], points[points.len() - 1]),
        other => panic!("{other:?}"),
    };
    let mut k = 0;
    for contour in &outline {
        let n = contour.len();
        for i in 0..n {
            assert_eq!(ends(&curves[k + i]).1, ends(&curves[k + (i + 1) % n]).0);
        }
        k += n;
    }
    // Every curve point lies within the box's reach.
    let [bl, br, _, tl] = c.map(|k| pos(&s.sketch, k));
    let (x, y) = ((br - bl).normalize(), (tl - bl).normalize());
    let (w, h) = ((br - bl).length(), (tl - bl).length());
    for curve in &curves {
        let (a, z) = ends(curve);
        for q in [a, z] {
            let (u, v) = ((q - bl).dot(x), (q - bl).dot(y));
            assert!(u > -0.05 * h && u < w + 0.05 * h, "{u} of {w}");
            assert!(v > -0.05 * h && v < 1.25 * h, "{v} of {h}");
        }
    }
}

#[test]
fn constraints_on_the_text_itself_are_refused() {
    let mut b = B::new();
    let (t, _, _) = text(&mut b, DVec2::ZERO, 0.004, 0.0, 2.0, false, square());
    let p = b.point(0.001, 0.001);
    b.c(Constraint::PointOnCurve { point: p, curve: t });
    let s = solve(&b.build(), &opts());
    match s.status {
        SolveStatus::Failed(SolveError::Invalid(m)) => {
            assert!(m.contains("is text; constrain its box"), "{m}")
        }
        other => panic!("{other:?}"),
    }
}

/// A change that breaks a text, named by what it breaks.
type Breakage = (&'static str, Box<dyn Fn(&mut Geometry<u32>)>);

#[test]
fn malformed_text_is_invalid() {
    let cases: Vec<Breakage> = vec![
        (
            "aspect ratio",
            Box::new(|g| {
                if let Geometry::Text { aspect, .. } = g {
                    *aspect = 0.0;
                }
            }),
        ),
        (
            "four distinct corner points",
            Box::new(|g| {
                if let Geometry::Text { corners, .. } = g {
                    corners[2] = corners[0];
                }
            }),
        ),
        (
            "does not meet the next",
            Box::new(|g| {
                if let Geometry::Text { outline, .. } = g {
                    outline[0][1].points[2].x += 1e-12;
                }
            }),
        ),
        (
            "degree 4",
            Box::new(|g| {
                if let Geometry::Text { outline, .. } = g {
                    outline[0][1].degree = 4;
                }
            }),
        ),
        (
            "fewer than two pieces",
            Box::new(|g| {
                if let Geometry::Text { outline, .. } = g {
                    outline[0].truncate(1);
                }
            }),
        ),
        (
            "out of range",
            Box::new(|g| {
                if let Geometry::Text { outline, .. } = g {
                    outline[0][0].points[1].y = f64::NAN;
                }
            }),
        ),
    ];
    for (why, break_it) in cases {
        let mut b = B::new();
        let (t, _, _) = text(&mut b, DVec2::ZERO, 0.004, 0.0, 2.0, false, square());
        let e = b.sketch.entities.iter_mut().find(|e| e.id == t).unwrap();
        break_it(&mut e.geometry);
        let s = solve(&b.build(), &opts());
        match s.status {
            SolveStatus::Failed(SolveError::Invalid(m)) => assert!(m.contains(why), "{why}: {m}"),
            other => panic!("{why}: {other:?}"),
        }
        assert!(resolve_curves(&b.build()).is_err(), "{why}");
    }
}

#[test]
fn a_collapsed_box_is_never_a_solution() {
    let mut b = B::new();
    let (_, c, _) = text(&mut b, DVec2::ZERO, 0.004, 0.0, 2.0, false, square());
    b.coincident(c[0], c[3]);
    let s = solve(&b.build(), &opts());
    assert!(matches!(s.status, SolveStatus::Failed(_)), "{:?}", s.status);
    if let SolveStatus::Failed(SolveError::Degenerate(m)) = &s.status {
        assert!(m.contains("text box"), "{m}");
    }
}

#[test]
fn mirror_images_of_text_are_mirrored_text() {
    let mut b = B::new();
    let (axis, _, _) = b.line(0.0, -1.0, 0.0, 1.0);
    b.c(Constraint::Vertical { line: axis });
    let (t, _, lines) = text(
        &mut b,
        DVec2::new(0.01, 0.0),
        0.004,
        0.0,
        2.0,
        false,
        square(),
    );
    let (ti, ci, lines_i) = text(
        &mut b,
        DVec2::new(-0.01, 0.0),
        0.004,
        0.0,
        2.0,
        true,
        square(),
    );
    let mut pairs = vec![(t, ti)];
    pairs.extend(lines.iter().copied().zip(lines_i));
    b.c(Constraint::Mirror {
        axis,
        pairs: pairs.clone(),
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    // The axis line's 4 (less vertical: 3) and the source text's 4.
    assert_eq!(s.dof, 3 + 4);
    // The image's bottom-left is the reflection of the source's.
    let p = pos(&s.sketch, ci[0]);
    assert!((p - DVec2::new(-0.01, 0.0)).length() < 1e-12, "{p}");
    // An unmirrored image is not a mirror image of text.
    let mut wrong = b.build();
    if let Some(e) = wrong.entities.iter_mut().find(|e| e.id == ti) {
        if let Geometry::Text { mirrored, .. } = &mut e.geometry {
            *mirrored = false;
        }
    }
    match solve(&wrong, &opts()).status {
        SolveStatus::Failed(SolveError::Invalid(m)) => {
            assert!(m.contains("not a mirrored copy"), "{m}")
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn pattern_copies_of_text_follow_their_source() {
    let mut b = B::new();
    let (dir, d0, d1) = b.line(0.0, 0.0, 0.02, 0.0);
    b.fix(d0);
    b.fix(d1);
    let (t, c, _) = text(
        &mut b,
        DVec2::new(0.0, 0.01),
        0.004,
        0.0,
        2.0,
        false,
        square(),
    );
    let (t1, c1, _) = text(
        &mut b,
        DVec2::new(0.02, 0.01),
        0.004,
        0.0,
        2.0,
        false,
        square(),
    );
    let (t2, c2, _) = text(
        &mut b,
        DVec2::new(0.04, 0.01),
        0.004,
        0.0,
        2.0,
        false,
        square(),
    );
    b.c(Constraint::LinearPattern {
        direction: dir,
        copies: vec![(t, t1, 1), (t, t2, 2)],
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    // Only the source text is free (the copies' boxes are implied).
    assert_eq!(s.dof, 4);
    for k in 0..4 {
        let d1 = pos(&s.sketch, c1[k]) - pos(&s.sketch, c[k]);
        let d2 = pos(&s.sketch, c2[k]) - pos(&s.sketch, c[k]);
        assert!((d1 - DVec2::new(0.02, 0.0)).length() < 1e-12);
        assert!((d2 - DVec2::new(0.04, 0.0)).length() < 1e-12);
    }
}
