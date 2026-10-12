//! The derived offset constraint (design note 0007 §10): images follow
//! their sources at the distance, add no degrees of freedom and no
//! redundancy, the distance drives them as a dimension, malformed offsets
//! are refused with a message naming the problem, and an offset that the
//! geometry has made invalid (a piece turned back, pieces crossing, a
//! smooth joint no longer smooth) fails the solve instead of being
//! returned.

mod common;

use common::gen::Rng;
use common::{assert_satisfied, offset_violation, pos, B};
use opendatum_sketch_solver::{
    solve, Constraint, ConstraintEntry, DVec2, Geometry, Sketch, SolveError, SolveOptions,
    SolveStatus, TangentKind,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

fn v(x: f64, y: f64) -> DVec2 {
    DVec2::new(x, y)
}

/// A closed polygon through `corners` (counter-clockwise): its lines, each
/// with its own ends, joined by coincident constraints. Returns the lines.
fn polygon(b: &mut B, corners: &[DVec2]) -> Vec<u32> {
    let n = corners.len();
    let mut lines = Vec::new();
    let mut ends = Vec::new();
    for k in 0..n {
        let (a, z) = (corners[k], corners[(k + 1) % n]);
        let (l, s, e) = b.line(a.x, a.y, z.x, z.y);
        lines.push(l);
        ends.push((s, e));
    }
    for k in 0..n {
        b.coincident(ends[k].1, ends[(k + 1) % n].0);
    }
    lines
}

/// The offset of a closed convex polygon by `d` (outwards when `out`):
/// corners where the offset lines meet. Adds the image points and lines
/// and returns the pieces.
fn offset_polygon(
    b: &mut B,
    corners: &[DVec2],
    lines: &[u32],
    d: f64,
    out: bool,
) -> Vec<(u32, u32, bool)> {
    let n = corners.len();
    // Counter-clockwise: the inside is on the left.
    let side = if out { -1.0 } else { 1.0 };
    let shifted = |k: usize| {
        let (a, z) = (corners[k], corners[(k + 1) % n]);
        let u = (z - a).normalize();
        let nrm = DVec2::new(-u.y, u.x) * side * d;
        (a + nrm, u)
    };
    let joins: Vec<DVec2> = (0..n)
        .map(|k| {
            let (p, u) = shifted((k + n - 1) % n);
            let (q, w) = shifted(k);
            p + u * ((q - p).perp_dot(w) / u.perp_dot(w))
        })
        .collect();
    let points: Vec<u32> = joins.iter().map(|j| b.point(j.x, j.y)).collect();
    (0..n)
        .map(|k| {
            let img = b.line_between(points[k], points[(k + 1) % n]);
            (lines[k], img, !out)
        })
        .collect()
}

fn offset_of(s: &Sketch<u32>, id: u32) -> (f64, Vec<(u32, u32, bool)>) {
    match &s.constraint(id).unwrap().constraint {
        Constraint::Offset {
            distance, pieces, ..
        } => (*distance, pieces.clone()),
        c => panic!("{c:?}"),
    }
}

/// A 4 × 2 rectangle held by horizontal and vertical constraints, offset
/// by `d` outwards or inwards. Returns the sketch, the rectangle's lines
/// and the offset's key.
fn rectangle(d: f64, out: bool) -> (Sketch<u32>, Vec<u32>, u32) {
    let mut b = B::new();
    let corners = [v(0.0, 0.0), v(4.0, 0.0), v(4.0, 2.0), v(0.0, 2.0)];
    let lines = polygon(&mut b, &corners);
    for (k, l) in lines.iter().enumerate() {
        b.c(if k % 2 == 0 {
            Constraint::Horizontal { line: *l }
        } else {
            Constraint::Vertical { line: *l }
        });
    }
    let pieces = offset_polygon(&mut b, &corners, &lines, d, out);
    let o = b.c(Constraint::Offset {
        distance: d,
        pieces,
        smooth: Vec::new(),
    });
    (b.build(), lines, o)
}

fn without_offset(s: &Sketch<u32>, o: u32) -> Sketch<u32> {
    let mut s = s.clone();
    s.constraints.retain(|c| c.id != o);
    s
}

#[test]
fn a_rectangle_offset_follows_it_with_no_freedom_added() {
    for out in [true, false] {
        let (mut s, _, o) = rectangle(0.5, out);
        // Images placed off their spots are pulled onto them.
        for e in s.entities.iter_mut().skip(16) {
            if let Geometry::Point { at } = &mut e.geometry {
                *at += v(0.03, -0.02);
            }
        }
        let alone = solve(&without_offset(&s, o), &opts());
        let out_s = solve(&s, &opts());
        assert_eq!(out_s.status, SolveStatus::Solved, "{:?}", out_s.redundant);
        // The rectangle: x, y, width, height; the images add nothing.
        assert_eq!(out_s.dof, 4);
        assert!(alone.dof > out_s.dof);
        assert_satisfied(&out_s.sketch);
        let (d, pieces) = offset_of(&out_s.sketch, o);
        assert!(offset_violation(&out_s.sketch, d, &pieces) < 1e-9);
        // The image corners are the rectangle's grown or shrunk by d.
        let corner = if out { v(-0.5, -0.5) } else { v(0.5, 0.5) };
        let img = pieces[0].1;
        let Geometry::Line { start, .. } = out_s.sketch.entity(img).unwrap().geometry else {
            panic!()
        };
        assert!((pos(&out_s.sketch, start) - corner).length() < 1e-9);
    }
}

#[test]
fn the_distance_is_a_dimension_and_drags_carry_the_images() {
    let (mut s, lines, o) = rectangle(0.5, true);
    // Changing the distance moves the images.
    for c in s.constraints.iter_mut().filter(|c| c.id == o) {
        assert_eq!(c.constraint.value(), Some(0.5));
        assert!(c.constraint.set_value(0.25));
    }
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    let (d, pieces) = offset_of(&out.sketch, o);
    assert_eq!(d, 0.25);
    assert!(offset_violation(&out.sketch, d, &pieces) < 1e-9);
    // Only the images moved: the rectangle is where it was.
    for l in &lines {
        assert_eq!(common::line(&out.sketch, *l), common::line(&s, *l));
    }
    // Dragging a corner of the source: the images follow.
    let Geometry::Line { start, .. } = out.sketch.entity(lines[0]).unwrap().geometry else {
        panic!()
    };
    let dragged = solve(
        &out.sketch,
        &SolveOptions {
            drag: vec![(start, v(-1.0, -0.5))],
            ..opts()
        },
    );
    assert_eq!(dragged.status, SolveStatus::Solved);
    assert!((pos(&dragged.sketch, start) - v(-1.0, -0.5)).length() < 1e-6);
    assert!(offset_violation(&dragged.sketch, d, &pieces) < 1e-9);
    assert_satisfied(&dragged.sketch);
}

/// Appends a constraint with a fresh key.
fn push(s: &mut Sketch<u32>, constraint: Constraint<u32>) -> u32 {
    let id = s
        .entities
        .iter()
        .map(|e| e.id)
        .chain(s.constraints.iter().map(|c| c.id))
        .max()
        .unwrap_or(0)
        + 1;
    s.constraints.push(ConstraintEntry { id, constraint });
    id
}

#[test]
fn a_dimension_on_an_image_of_a_dimensioned_source_is_redundant() {
    let (mut s, lines, o) = rectangle(0.5, true);
    let len = push(
        &mut s,
        Constraint::Length {
            line: lines[0],
            value: 4.0,
        },
    );
    let (_, pieces) = offset_of(&s, o);
    let on_image = push(
        &mut s,
        Constraint::Length {
            line: pieces[0].1,
            value: 5.0,
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Redundant);
    assert_eq!(out.redundant, vec![on_image]);
    // A dimension on the image alone drives the source.
    s.constraints.retain(|c| c.id != len);
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    let (s0, s1) = common::line(&out.sketch, lines[0]);
    assert!(((s1 - s0).length() - 4.0).abs() < 1e-9);
}

#[test]
fn a_chain_with_a_smooth_joint_and_free_ends() {
    // A line along y = 0 from (−3, 0) to (0, 0), continued tangentially by
    // the quarter arc about (0, 1) to (1, 1), then a line leaning on to
    // (2, 3): a smooth joint, a corner and two free ends. Offset to the
    // right (below and outside) by 0.25.
    let mut b = B::new();
    let (l0, p0, p1) = b.line(-3.0, 0.0, 0.0, 0.0);
    let (arc, _c, a0, a1) = b.arc(v(0.0, 1.0), 1.0, -std::f64::consts::FRAC_PI_2, 0.0);
    let (l1, q0, _q1) = b.line(1.0, 1.0, 2.0, 3.0);
    b.coincident(p1, a0);
    b.coincident(a1, q0);
    b.c(Constraint::Tangent {
        a: l0,
        b: arc,
        kind: TangentKind::AtPoint { point: a0 },
    });
    b.c(Constraint::Horizontal { line: l0 });
    b.fix(p0);
    let d = 0.25;
    // Images: (−3, −0.25) → (0, −0.25) → arc of radius 1.25 about (0, 1) to
    // where it meets the last line's offset (its right, by 0.25), nearest
    // the corner → on along that offset to the normal at (2, 3). Placed
    // roughly; the solve puts them exactly.
    let u = v(1.0, 2.0).normalize();
    let right = v(u.y, -u.x) * d;
    let along = v(1.0, 1.0) + right;
    let w = along - v(0.0, 1.0);
    let t = -w.dot(u) + (w.dot(u).powi(2) - (w.length_squared() - 1.25 * 1.25)).sqrt();
    let corner = along + u * t;
    let free = v(2.0, 3.0) + right;
    let j0 = b.point(-3.0, -0.25);
    let j1 = b.point(0.0, -0.25);
    let j2 = b.point(corner.x + 0.01, corner.y - 0.02);
    let j3 = b.point(free.x - 0.01, free.y);
    let i0 = b.line_between(j0, j1);
    let Geometry::Arc { center, .. } = b.sketch.entity(arc).unwrap().geometry else {
        panic!()
    };
    let i1 = b.entity(Geometry::Arc {
        center,
        start: j1,
        end: j2,
    });
    let i2 = b.line_between(j2, j3);
    let pieces = vec![(l0, i0, false), (arc, i1, false), (l1, i2, false)];
    let o = b.c(Constraint::Offset {
        distance: d,
        pieces: pieces.clone(),
        smooth: vec![j1],
    });
    let s = b.build();
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    assert_satisfied(&out.sketch);
    // The sources' freedom only: the first line's length, the arc's radius
    // and sweep, the last line's far end. The images add none.
    assert_eq!(out.dof, 5);
    // Without the offset the images are free (four points, the image arc's
    // own equation tying its ends).
    assert_eq!(solve(&without_offset(&s, o), &opts()).dof, 5 + 2 * 4 - 1);
    assert!(((corner - v(0.0, 1.0)).length() - 1.25).abs() < 1e-12);
    for (p, want) in [
        (j0, v(-3.0, -0.25)),
        (j1, v(0.0, -0.25)),
        (j2, corner),
        (j3, free),
    ] {
        assert!(
            (pos(&out.sketch, p) - want).length() < 1e-9,
            "{p} {:?} {want:?}",
            pos(&out.sketch, p)
        );
    }
    // Making the smooth joint a corner breaks the offset: the solve fails
    // and names it, rather than returning a broken offset.
    let mut broken = out.sketch.clone();
    broken.constraints.retain(|c| {
        !matches!(
            c.constraint,
            Constraint::Tangent { .. } | Constraint::Horizontal { .. }
        )
    });
    let dragged = solve(
        &broken,
        &SolveOptions {
            drag: vec![(p1, v(0.0, -0.5))],
            ..opts()
        },
    );
    match dragged.status {
        SolveStatus::Failed(SolveError::Degenerate(m)) => {
            assert!(m.contains("tangentially"), "{m}")
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_circle_offset_inwards() {
    let mut b = B::new();
    let (k, c) = b.circle(1.0, 2.0, 3.0);
    b.c(Constraint::Radius {
        curve: k,
        value: 3.0,
    });
    let img = b.entity(Geometry::Circle {
        center: c,
        radius: 2.5,
    });
    let o = b.c(Constraint::Offset {
        distance: 1.0,
        pieces: vec![(k, img, true)],
        smooth: Vec::new(),
    });
    let s = b.build();
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert_eq!(out.dof, 2);
    let (d, pieces) = offset_of(&out.sketch, o);
    assert!(offset_violation(&out.sketch, d, &pieces) < 1e-9);
    // Too far inwards: no circle left.
    let mut s2 = out.sketch.clone();
    for c in s2.constraints.iter_mut().filter(|c| c.id == o) {
        c.constraint.set_value(3.5);
    }
    assert!(matches!(solve(&s2, &opts()).status, SolveStatus::Failed(_)));
}

#[test]
fn too_large_a_distance_fails_naming_the_offset() {
    // Inwards by more than half the rectangle's height: the short sides'
    // offsets would point the other way.
    let (mut s, _, o) = rectangle(0.5, false);
    for c in s.constraints.iter_mut().filter(|c| c.id == o) {
        c.constraint.set_value(1.5);
    }
    match solve(&s, &opts()).status {
        SolveStatus::Failed(SolveError::Degenerate(m)) => {
            assert!(m.contains("too large"), "{m}");
            assert!(m.contains(&format!("offset {o}")), "{m}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn offset_pieces_that_cross_fail() {
    // Two separate crossing lines offset to their left with free ends:
    // their images cross.
    let mut b = B::new();
    let (l0, ..) = b.line(0.0, 0.0, 4.0, 0.0);
    let (l1, ..) = b.line(2.0, -2.0, 2.0, 2.0);
    let (i0, ..) = b.line(0.0, 0.5, 4.0, 0.5);
    let (i1, ..) = b.line(1.5, -2.0, 1.5, 2.0);
    b.c(Constraint::Offset {
        distance: 0.5,
        pieces: vec![(l0, i0, true), (l1, i1, true)],
        smooth: Vec::new(),
    });
    match solve(&b.build(), &opts()).status {
        SolveStatus::Failed(SolveError::Degenerate(m)) => assert!(m.contains("cross"), "{m}"),
        other => panic!("{other:?}"),
    }
}

/// `true` if solving `s` fails as invalid input with a message containing
/// `needle`.
fn refused(s: &Sketch<u32>, needle: &str) {
    match solve(s, &opts()).status {
        SolveStatus::Failed(SolveError::Invalid(m)) => assert!(m.contains(needle), "{m}"),
        other => panic!("expected a refusal mentioning {needle:?}, got {other:?}"),
    }
}

#[test]
fn malformed_offsets_are_refused() {
    let base = || {
        let mut b = B::new();
        let (l0, _, e0) = b.line(0.0, 0.0, 4.0, 0.0);
        let (l1, s1, _) = b.line(4.0, 0.0, 4.0, 2.0);
        b.coincident(e0, s1);
        let (k, c) = b.circle(10.0, 0.0, 1.0);
        let (arc, ..) = b.arc(v(20.0, 0.0), 1.0, 0.0, 1.0);
        (b, l0, l1, k, c, arc)
    };
    let with = |b: &mut B, pieces: Vec<(u32, u32, bool)>, smooth: Vec<u32>, d: f64| {
        b.c(Constraint::Offset {
            distance: d,
            pieces,
            smooth,
        });
        b.build()
    };
    // No pieces; a non-positive distance; a curve its own offset.
    let (mut b, l0, ..) = base();
    refused(&with(&mut b, vec![], vec![], 1.0), "no curves");
    let (mut b, l0b, ..) = base();
    let (i, ..) = b.line(0.0, 1.0, 4.0, 1.0);
    refused(&with(&mut b, vec![(l0b, i, true)], vec![], 0.0), "zero");
    let (mut b, ..) = base();
    refused(
        &with(&mut b, vec![(l0, l0, true)], vec![], 1.0),
        "its own offset",
    );
    // Different kinds; an image circle with a centre of its own.
    let (mut b, l0, _, k, ..) = base();
    refused(
        &with(&mut b, vec![(l0, k, true)], vec![], 1.0),
        "not two lines",
    );
    let (mut b, _, _, k, ..) = base();
    let (k2, _) = b.circle(10.0, 0.0, 0.5);
    refused(
        &with(&mut b, vec![(k, k2, true)], vec![], 0.5),
        "share its centre",
    );
    // An image end that is a source point.
    let (mut b, _, l1, ..) = base();
    let Geometry::Line { end, .. } = b.sketch.entity(l1).unwrap().geometry else {
        panic!()
    };
    let p = b.point(3.0, 0.0);
    let i = b.line_between(p, end);
    refused(
        &with(&mut b, vec![(l1, i, true)], vec![], 1.0),
        "both an image and a source",
    );
    // Images meeting where their sources do not.
    let (mut b, l0, _, _, _, arc) = base();
    let (p, q, r) = (b.point(0.0, 1.0), b.point(4.0, 1.0), b.point(5.0, 1.0));
    let i0 = b.line_between(p, q);
    let i1 = b.line_between(q, r);
    let Geometry::Arc { center, .. } = b.sketch.entity(arc).unwrap().geometry else {
        panic!()
    };
    let ia = b.entity(Geometry::Arc {
        center,
        start: q,
        end: r,
    });
    refused(
        &with(&mut b, vec![(l0, i0, true), (arc, ia, true)], vec![], 1.0),
        "meet where their sources do not",
    );
    let _ = i1;
    // A smooth point that ends no image; a point ending three images.
    let (mut b, l0, l1, ..) = base();
    let (i0, ..) = b.line(0.0, 1.0, 3.0, 1.0);
    let stray = b.point(9.0, 9.0);
    refused(
        &with(&mut b, vec![(l0, i0, true)], vec![stray], 1.0),
        "not an end",
    );
    let (mut b, l0, l1b, ..) = base();
    let (l2, s2, _) = b.line(4.0, 0.0, 6.0, 0.0);
    let Geometry::Line { end: e0, .. } = b.sketch.entity(l0).unwrap().geometry else {
        panic!()
    };
    b.coincident(e0, s2);
    let (p, q, r, t) = (
        b.point(0.0, 1.0),
        b.point(3.0, 1.0),
        b.point(3.0, 2.0),
        b.point(6.0, 1.0),
    );
    let (i0, i1, i2) = (
        b.line_between(p, q),
        b.line_between(q, r),
        b.line_between(q, t),
    );
    refused(
        &with(
            &mut b,
            vec![(l0, i0, true), (l1b, i1, true), (l2, i2, true)],
            vec![],
            1.0,
        ),
        "more than two",
    );
    let _ = l1;
}

#[test]
fn deleting_a_piece_keeps_the_rest_offset() {
    let (s, _, o) = rectangle(0.5, true);
    let c = s.constraint(o).unwrap().constraint.clone();
    let Constraint::Offset { pieces, .. } = &c else {
        panic!()
    };
    let gone = pieces[1].1;
    let kept = c.without(|k| k == gone).unwrap();
    let mut s2 = s.clone();
    s2.entities.retain(|e| e.id != gone);
    for e in s2.constraints.iter_mut().filter(|e| e.id == o) {
        e.constraint = kept.clone();
    }
    let out = solve(&s2, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    let (d, pieces) = offset_of(&out.sketch, o);
    assert_eq!(pieces.len(), 3);
    assert!(offset_violation(&out.sketch, d, &pieces) < 1e-9);
    // Deleting every piece's source deletes the offset.
    let all: Vec<u32> = pieces.iter().map(|(s, ..)| *s).collect();
    assert!(kept.without(|k| all.contains(&k)).is_none());
}

/// A convex polygon with `n` corners at random radii about `c`.
fn random_polygon(rng: &mut Rng, n: usize) -> Vec<DVec2> {
    let c = v(rng.range(-5.0, 5.0), rng.range(-5.0, 5.0));
    let r = rng.range(1.0, 4.0);
    (0..n)
        .map(|k| {
            let a = std::f64::consts::TAU * (k as f64 + rng.range(-0.2, 0.2)) / n as f64;
            c + DVec2::from_angle(a) * r
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Random convex polygons offset outwards, or inwards by less than
    /// their inradius allows, with images placed off their spots: the
    /// offset adds no freedom and no redundancy and holds after the solve;
    /// a polygon dragged afterwards keeps it holding, or the solve fails
    /// naming the offset.
    #[test]
    fn random_polygon_offsets(seed in any::<u64>(), n in 3usize..9, out in any::<bool>(), frac in 0.01..0.3f64) {
        let mut rng = Rng(seed);
        let corners = random_polygon(&mut rng, n);
        let mut b = B::new();
        let lines = polygon(&mut b, &corners);
        let size = corners.iter().map(|c| (*c - corners[0]).length()).fold(0.0, f64::max);
        let d = frac * size * 0.2;
        let alone = solve(&b.build(), &opts());
        prop_assume!(alone.status == SolveStatus::Solved);
        let pieces = offset_polygon(&mut b, &corners, &lines, d, out);
        let o = b.c(Constraint::Offset { distance: d, pieces, smooth: Vec::new() });
        let mut s = b.build();
        let k = s.entities.len() - 2 * n;
        for e in s.entities.iter_mut().skip(k) {
            if let Geometry::Point { at } = &mut e.geometry {
                *at += v(rng.range(-0.01, 0.01), rng.range(-0.01, 0.01)) * d;
            }
        }
        let res = solve(&s, &opts());
        match &res.status {
            SolveStatus::Solved => {
                prop_assert_eq!(res.dof, alone.dof);
                let (d, pieces) = offset_of(&res.sketch, o);
                prop_assert!(offset_violation(&res.sketch, d, &pieces) < 1e-9);
                assert_satisfied(&res.sketch);
            }
            // Inwards, a short side can be swallowed even below the bound.
            SolveStatus::Failed(SolveError::Degenerate(m)) if !out => {
                prop_assert!(m.contains("offset"), "{}", m);
            }
            other => prop_assert!(false, "{:?}", other),
        }
    }
}
