//! Ellipses and elliptical arcs in the solver (design note 0008): degrees
//! of freedom, points on them, tangency to lines and at shared ends,
//! derived images, refusals, dragging, collapse, and constructive property
//! tests, all checked with the independent formulas of `common`.

mod common;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use common::gen::{generate_with, perturb};
use common::{assert_satisfied, bits, ell, pos, resolved, size, worst, B};
use opendatum_sketch_solver::{
    solve, Constraint, DVec2, EntityState, Geometry, ResolvedCurve, SolveError, SolveOptions,
    SolveStatus, TangentKind,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// An ellipse with its own points: (ellipse, centre, major end, minor end).
fn ellipse(b: &mut B, c: DVec2, a: DVec2, ratio: f64) -> (u32, u32, u32, u32) {
    let kc = b.point(c.x, c.y);
    let km = b.point(c.x + a.x, c.y + a.y);
    let n = c + a.perp() * ratio;
    let kn = b.point(n.x, n.y);
    let k = b.entity(Geometry::Ellipse {
        center: kc,
        major: km,
        minor: kn,
    });
    (k, kc, km, kn)
}

/// An elliptical arc from eccentric anomaly `t0` to `t1` (counter-clockwise
/// in the plane, the minor axis left of the major): (arc, centre, major
/// end, minor end, start, end).
fn arc(b: &mut B, c: DVec2, a: DVec2, ratio: f64, t0: f64, t1: f64) -> [u32; 6] {
    let bb = a.perp() * ratio;
    let at = |t: f64| c + a * t.cos() + bb * t.sin();
    let kc = b.point(c.x, c.y);
    let km = b.point(c.x + a.x, c.y + a.y);
    let kn = b.point(c.x + bb.x, c.y + bb.y);
    let (p0, p1) = (at(t0), at(t1));
    let ks = b.point(p0.x, p0.y);
    let ke = b.point(p1.x, p1.y);
    let k = b.entity(Geometry::EllipticalArc {
        center: kc,
        major: km,
        minor: kn,
        start: ks,
        end: ke,
    });
    [k, kc, km, kn, ks, ke]
}

#[test]
fn an_ellipse_has_five_freedoms_and_an_arc_seven() {
    let mut b = B::new();
    let (k, ..) = ellipse(&mut b, DVec2::new(0.1, 0.2), DVec2::new(0.5, 0.2), 0.6);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 5);
    assert_eq!(s.entities[&k], EntityState::Free);

    let mut b = B::new();
    let [k, ..] = arc(&mut b, DVec2::ZERO, DVec2::new(0.4, 0.0), 0.5, 0.3, 2.5);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 7);
    assert_eq!(s.entities[&k], EntityState::Free);
    assert!(s.redundant.is_empty());
}

#[test]
fn axes_that_are_not_perpendicular_are_made_so() {
    let mut b = B::new();
    let kc = b.point(0.0, 0.0);
    let km = b.point(1.0, 0.0);
    let kn = b.point(0.2, 0.5);
    let k = b.entity(Geometry::Ellipse {
        center: kc,
        major: km,
        minor: kn,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    let (_, a, bb) = ell(&s.sketch, k);
    assert!(a.normalize().dot(bb.normalize()).abs() < 1e-12);
    assert_satisfied(&s.sketch);
}

#[test]
fn a_fully_constrained_ellipse_through_its_axis_lines() {
    let mut b = B::new();
    let (_, kc, km, kn) = ellipse(&mut b, DVec2::new(0.1, 0.1), DVec2::new(0.9, 0.3), 0.5);
    let major = b.line_between(kc, km);
    let minor = b.line_between(kc, kn);
    b.fix(kc);
    b.c(Constraint::Horizontal { line: major });
    b.c(Constraint::Length {
        line: major,
        value: 0.8,
    });
    b.c(Constraint::Length {
        line: minor,
        value: 0.3,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert_eq!(s.dof, 0);
    assert!(s.redundant.is_empty());
    assert!(s.entities.values().all(|&e| e == EntityState::Constrained));
    assert_satisfied(&s.sketch);
    // A perpendicular constraint on the axes is implied by the ellipse.
    b.c(Constraint::Perpendicular { a: major, b: minor });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant.len(), 1);
}

#[test]
fn points_on_ellipses_and_arcs() {
    let mut b = B::new();
    let (k, ..) = ellipse(&mut b, DVec2::ZERO, DVec2::new(0.6, 0.0), 0.5);
    let p = b.point(0.9, 0.7);
    b.c(Constraint::PointOnCurve { point: p, curve: k });
    let [k2, ..] = arc(
        &mut b,
        DVec2::new(2.0, 0.0),
        DVec2::new(0.0, 0.5),
        0.4,
        0.0,
        2.0,
    );
    let q = b.point(2.5, 1.0);
    b.c(Constraint::PointOnCurve {
        point: q,
        curve: k2,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 5 + 2 - 1 + 7 + 2 - 1);
    assert_satisfied(&s.sketch);
    // The centre cannot lie on its own ellipse.
    let mut b2 = b.clone();
    let centre = match b2.sketch.entity(k).unwrap().geometry {
        Geometry::Ellipse { center, .. } => center,
        _ => unreachable!(),
    };
    b2.c(Constraint::PointOnCurve {
        point: centre,
        curve: k,
    });
    let s = solve(&b2.build(), &opts());
    assert!(
        matches!(&s.status, SolveStatus::Failed(SolveError::Invalid(m)) if m.contains("centre")),
        "{:?}",
        s.status
    );
}

/// The line's signed distance (left positive) from the ellipse's centre
/// and the ellipse's extent towards the line, measured by dense sampling
/// (independently of both closed forms).
fn line_gap(s: &opendatum_sketch_solver::Sketch<u32>, l: u32, k: u32) -> (f64, f64) {
    let (s0, s1) = common::line(s, l);
    let n = (s1 - s0).perp().normalize();
    let (c, a, bb) = ell(s, k);
    let extent = (0..200_000)
        .map(|i| {
            let t = TAU * i as f64 / 200_000.0;
            n.dot(a * t.cos() + bb * t.sin())
        })
        .fold(f64::NEG_INFINITY, f64::max);
    (n.dot(c - s0), extent)
}

#[test]
fn lines_become_tangent_on_either_side() {
    for center_left in [true, false] {
        let mut b = B::new();
        let (k, ..) = ellipse(&mut b, DVec2::ZERO, DVec2::new(0.7, 0.4), 0.45);
        let y = if center_left { -0.9 } else { 0.9 };
        let (l, ..) = b.line(-1.5, y, 1.5, y + 0.2);
        b.c(Constraint::Tangent {
            a: l,
            b: k,
            kind: TangentKind::LineCircle { center_left },
        });
        let s = solve(&b.build(), &opts());
        assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
        assert_satisfied(&s.sketch);
        let (signed, extent) = line_gap(&s.sketch, l, k);
        let want = if center_left { extent } else { -extent };
        // Dense sampling resolves the support to ~(Δθ)²/2 of the semi-axes.
        assert!((signed - want).abs() < 1e-9, "{signed} vs {want}");
    }
}

#[test]
fn arcs_continue_tangentially_into_lines_and_other_arcs() {
    let mut b = B::new();
    let [k, .., ke] = arc(&mut b, DVec2::ZERO, DVec2::new(0.5, 0.1), 0.6, -1.0, 1.2);
    let (l, ls, _) = b.line(0.2, 0.6, -0.4, 1.2);
    b.coincident(ke, ls);
    b.c(Constraint::Tangent {
        a: k,
        b: l,
        kind: TangentKind::AtPoint { point: ls },
    });
    // A second elliptical arc leaving the first one's start.
    let [k2, .., ks2, _] = arc(
        &mut b,
        DVec2::new(1.0, -0.6),
        DVec2::new(0.3, 0.0),
        0.7,
        2.0,
        4.0,
    );
    let start = match b.sketch.entity(k).unwrap().geometry {
        Geometry::EllipticalArc { start, .. } => start,
        _ => unreachable!(),
    };
    b.coincident(start, ks2);
    b.c(Constraint::Tangent {
        a: k,
        b: k2,
        kind: TangentKind::AtPoint { point: start },
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!(s.redundant.is_empty());
    assert_satisfied(&s.sketch);
    // An arc end that is not at the point, or a full ellipse, is refused.
    let mut b2 = b.clone();
    let (e, ..) = ellipse(&mut b2, DVec2::new(3.0, 0.0), DVec2::new(0.3, 0.0), 0.5);
    b2.c(Constraint::Tangent {
        a: e,
        b: l,
        kind: TangentKind::AtPoint { point: ls },
    });
    let s = solve(&b2.build(), &opts());
    assert!(matches!(
        s.status,
        SolveStatus::Failed(SolveError::Invalid(_))
    ));
}

#[test]
fn arcs_of_one_ellipse_meet_at_a_coincidence_without_redundancy() {
    // The pieces of a split arc (design note 0008 §8): two arcs on the same
    // centre and axis points, the first's end coincident with the second's
    // start. Each end is on the ellipse; the joint only once.
    let mut b = B::new();
    let [k1, kc, km, kn, _, e1] = arc(&mut b, DVec2::ZERO, DVec2::new(0.6, 0.2), 0.5, 0.3, 1.4);
    let a = DVec2::new(0.6, 0.2);
    let bb = a.perp() * 0.5;
    let at = |t: f64| a * t.cos() + bb * t.sin();
    let (p, q) = (at(1.4), at(2.9));
    let s2 = b.point(p.x, p.y);
    let e2 = b.point(q.x, q.y);
    let k2 = b.entity(Geometry::EllipticalArc {
        center: kc,
        major: km,
        minor: kn,
        start: s2,
        end: e2,
    });
    b.coincident(e1, s2);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!(s.redundant.is_empty(), "{:?}", s.redundant);
    // The ellipse's 5, the first arc's ends, the second's far end.
    assert_eq!(s.dof, 5 + 3);
    assert_satisfied(&s.sketch);
    for k in [k1, k2] {
        assert!(matches!(
            resolved(&s.sketch, k),
            ResolvedCurve::EllipticalArc { .. }
        ));
    }
    // Dragged, the joint stays on the ellipse.
    let mut o = opts();
    o.drag = vec![(s2, DVec2::new(-0.3, 0.6))];
    let s = solve(&s.sketch, &o);
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!(common::dist_to_ellipse(&s.sketch, k1, pos(&s.sketch, e1)) < 1e-12);
    assert!((pos(&s.sketch, e1) - pos(&s.sketch, s2)).length() < 1e-12);

    // An end made coincident with an axis end lies on the ellipse through
    // it: not redundant either.
    let mut b = B::new();
    let [_, _, km, _, ks, _] = arc(&mut b, DVec2::ZERO, DVec2::new(0.6, 0.0), 0.5, 0.2, 2.0);
    b.coincident(ks, km);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!(s.redundant.is_empty());
    assert_eq!(s.dof, 5 + 1);
    assert!((pos(&s.sketch, ks) - pos(&s.sketch, km)).length() < 1e-12);

    // A point-on-curve constraint at the joint says it again: redundant.
    let mut b = B::new();
    let [k1, _, _, _, _, e1] = arc(&mut b, DVec2::ZERO, DVec2::new(0.6, 0.2), 0.5, 0.3, 1.4);
    b.c(Constraint::PointOnCurve {
        point: e1,
        curve: k1,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.redundant.len(), 1, "{:?}", s.status);
}

#[test]
fn what_ellipses_do_not_take_is_refused() {
    let mut b = B::new();
    let (k, kc, ..) = ellipse(&mut b, DVec2::ZERO, DVec2::new(0.5, 0.0), 0.5);
    let (circle, _) = b.circle(1.5, 0.0, 0.3);
    let (k2, ..) = ellipse(&mut b, DVec2::new(-2.0, 0.0), DVec2::new(0.5, 0.0), 0.5);
    let (l, ..) = b.line(-1.0, 1.0, 1.0, 1.0);
    let refused = [
        Constraint::Radius {
            curve: k,
            value: 0.5,
        },
        Constraint::Diameter {
            curve: k,
            value: 0.5,
        },
        Constraint::Equal { a: k, b: k2 },
        Constraint::Concentric { a: k, b: circle },
        Constraint::Normal { line: l, curve: k },
        Constraint::Tangent {
            a: k,
            b: circle,
            kind: TangentKind::CircleCircle { internal: false },
        },
        Constraint::Tangent {
            a: k,
            b: k2,
            kind: TangentKind::CircleCircle { internal: true },
        },
        Constraint::Distance {
            a: kc,
            b: k,
            kind: opendatum_sketch_solver::DistanceKind::Aligned,
            value: 0.2,
        },
    ];
    for c in refused {
        let mut b2 = b.clone();
        b2.c(c.clone());
        let s = solve(&b2.build(), &opts());
        assert!(
            matches!(s.status, SolveStatus::Failed(SolveError::Invalid(_))),
            "{c:?}: {:?}",
            s.status
        );
        assert_eq!(bits(&s.sketch), bits(&b2.build()), "{c:?}: input kept");
    }
    // Bad definitions.
    let mut b3 = B::new();
    let p = b3.point(0.0, 0.0);
    let q = b3.point(1.0, 0.0);
    b3.entity(Geometry::Ellipse {
        center: p,
        major: q,
        minor: q,
    });
    let s = solve(&b3.build(), &opts());
    assert!(matches!(
        s.status,
        SolveStatus::Failed(SolveError::Invalid(_))
    ));
}

#[test]
fn mirrored_and_patterned_ellipses_add_no_freedom() {
    let mut b = B::new();
    let (axis, ..) = b.line(0.0, -2.0, 0.0, 2.0);
    let [src, ..] = arc(
        &mut b,
        DVec2::new(0.8, 0.3),
        DVec2::new(0.4, 0.2),
        0.5,
        0.2,
        2.8,
    );
    let (full, ..) = ellipse(&mut b, DVec2::new(1.5, -1.0), DVec2::new(0.3, 0.1), 0.6);
    let (dir, ..) = b.line(2.5, -1.5, 3.0, -1.5);
    let free = solve(&b.build(), &opts()).dof;
    // The image: the reflection, its ends crossed over.
    let reflect = |p: DVec2| DVec2::new(-p.x, p.y);
    let s0 = b.build();
    let pts: Vec<DVec2> = s0
        .entity(src)
        .unwrap()
        .geometry
        .referenced_points()
        .iter()
        .map(|&q| reflect(pos(&s0, q)))
        .collect();
    let keys: Vec<u32> = pts
        .iter()
        .map(|p| b.point(p.x + 0.05, p.y - 0.03))
        .collect();
    let img = b.entity(Geometry::EllipticalArc {
        center: keys[0],
        major: keys[1],
        minor: keys[2],
        start: keys[4],
        end: keys[3],
    });
    let copy_pts: Vec<DVec2> = {
        let s = b.build();
        s.entity(full)
            .unwrap()
            .geometry
            .referenced_points()
            .iter()
            .map(|&q| pos(&s, q) + DVec2::new(0.5, 0.02))
            .collect()
    };
    let ck: Vec<u32> = copy_pts.iter().map(|p| b.point(p.x, p.y)).collect();
    let copy = b.entity(Geometry::Ellipse {
        center: ck[0],
        major: ck[1],
        minor: ck[2],
    });
    b.c(Constraint::Mirror {
        axis,
        pairs: vec![(src, img)],
    });
    b.c(Constraint::LinearPattern {
        direction: dir,
        copies: vec![(full, copy, 1)],
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!(s.redundant.is_empty());
    assert_eq!(s.dof, free, "images add no degrees of freedom");
    assert_satisfied(&s.sketch);
    // The image is the same arc, reflected: same sweep, opposite sense of θ.
    let (
        ResolvedCurve::EllipticalArc { sweep_param: a, .. },
        ResolvedCurve::EllipticalArc { sweep_param: m, .. },
    ) = (resolved(&s.sketch, src), resolved(&s.sketch, img))
    else {
        panic!("not elliptical arcs");
    };
    assert!((a + m).abs() < 1e-9, "{a} vs {m}");
}

#[test]
fn dragging_an_axis_end_keeps_the_axes_perpendicular() {
    let mut b = B::new();
    let [k, kc, km, kn, ..] = arc(&mut b, DVec2::ZERO, DVec2::new(0.5, 0.0), 0.5, 0.5, 2.5);
    b.fix(kc);
    let target = DVec2::new(0.4, 0.4);
    let s = solve(
        &b.build(),
        &SolveOptions {
            drag: vec![(km, target)],
            ..SolveOptions::default()
        },
    );
    assert_eq!(s.status, SolveStatus::Solved);
    // The dragged point moves last (weight 1e6 against 1 for the others,
    // design note 0005 §7), so it lands within about 1e-6 of the other
    // points' correction: here the minor end turns through 45°.
    let miss = (pos(&s.sketch, km) - target).length();
    assert!(miss < 1e-6, "{miss}");
    assert_satisfied(&s.sketch);
    // Perpendicular to the model's angular tolerance (the solver's
    // contract, 1e-9 rad), not to rounding: Newton stops once it holds.
    let (_, a, bb) = ell(&s.sketch, k);
    let cos = a.dot(bb) / (a.length() * bb.length());
    assert!(cos.abs() < 1e-9, "{cos:e}");
    assert!(pos(&s.sketch, kn).is_finite());
}

#[test]
fn an_axis_held_at_zero_length_is_a_conflict_not_a_result() {
    let mut b = B::new();
    let (_, kc, _, kn) = ellipse(&mut b, DVec2::ZERO, DVec2::new(0.5, 0.0), 0.5);
    let minor = b.line_between(kc, kn);
    b.c(Constraint::Length {
        line: minor,
        value: 0.3,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    // Coincident centre and minor end: the ellipse collapses.
    b.coincident(kc, kn);
    let input = b.build();
    let s = solve(&input, &opts());
    assert!(matches!(s.status, SolveStatus::Failed(_)), "{:?}", s.status);
    assert_eq!(bits(&s.sketch), bits(&input));
}

#[test]
fn resolved_arcs_run_counter_clockwise_either_handedness() {
    for flip in [false, true] {
        let mut b = B::new();
        let c = DVec2::new(0.2, -0.1);
        let a = DVec2::new(0.5, 0.0);
        let bb = if flip {
            DVec2::new(0.0, -0.3)
        } else {
            DVec2::new(0.0, 0.3)
        };
        let kc = b.point(c.x, c.y);
        let km = b.point(c.x + a.x, c.y + a.y);
        let kn = b.point(c.x + bb.x, c.y + bb.y);
        // From the right vertex to the top vertex, counter-clockwise.
        let ks = b.point(c.x + 0.5, c.y);
        let ke = b.point(c.x, c.y + 0.3);
        let k = b.entity(Geometry::EllipticalArc {
            center: kc,
            major: km,
            minor: kn,
            start: ks,
            end: ke,
        });
        let ResolvedCurve::EllipticalArc {
            start_param,
            sweep_param,
            major,
            minor,
            ..
        } = resolved(&b.build(), k)
        else {
            panic!()
        };
        let at = |t: f64| c + major * t.cos() + minor * t.sin();
        assert!((at(start_param) - DVec2::new(0.7, -0.1)).length() < 1e-15);
        assert!((at(start_param + sweep_param) - DVec2::new(0.2, 0.2)).length() < 1e-15);
        assert!((sweep_param.abs() - FRAC_PI_2).abs() < 1e-15);
        assert_eq!(sweep_param > 0.0, !flip);
        // Its middle is in the first quadrant: the counter-clockwise way.
        let mid = at(start_param + 0.5 * sweep_param) - c;
        assert!(mid.x > 0.0 && mid.y > 0.0);
    }
    // An arc whose ends meet is degenerate.
    let mut b = B::new();
    let [k, .., ks, ke] = arc(&mut b, DVec2::ZERO, DVec2::new(0.5, 0.0), 0.5, 0.0, PI);
    let at = b.sketch.point(ks).unwrap();
    for e in b.sketch.entities.iter_mut() {
        if e.id == ke {
            e.geometry = Geometry::Point { at };
        }
    }
    assert!(opendatum_sketch_solver::resolve_curves(&b.build()).is_err());
    let _ = k;
}

#[test]
fn degenerate_ellipses_do_not_resolve() {
    let message = |c: DVec2, a: DVec2, ratio: f64, bend: DVec2| {
        let mut b = B::new();
        let (_, _, _, kn) = ellipse(&mut b, c, a, ratio);
        for e in b.sketch.entities.iter_mut() {
            if e.id == kn {
                e.geometry = Geometry::Point { at: c + bend };
            }
        }
        match opendatum_sketch_solver::resolve_curves(&b.build()) {
            Err(SolveError::Degenerate(m)) => m,
            other => panic!("{other:?}"),
        }
    };
    let (c, a) = (DVec2::new(0.1, 0.2), DVec2::new(0.5, 0.0));
    assert!(message(c, a, 0.5, DVec2::new(0.25, 0.0)).contains("parallel axes"));
    assert!(message(c, a, 0.5, DVec2::ZERO).contains("zero length"));
}

proptest! {
    #![proptest_config(ProptestConfig {
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// The generator's sketches with ellipses satisfy their constraints.
    #[test]
    fn generated_sketches_with_ellipses_hold(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
    ) {
        let scale = 10f64.powi(e);
        let g = generate_with(seed, motifs, scale, DVec2::new(ox, oy), true);
        let (w, what) = worst(&g.sketch, size(&g.sketch));
        prop_assert!(w <= 1e-11, "generator: violation {:e} at {}", w, what);
    }

    /// Satisfied sketches with ellipses stay where they are.
    #[test]
    fn satisfied_sketches_with_ellipses_stay_put(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
    ) {
        let scale = 10f64.powi(e);
        let g = generate_with(seed, motifs, scale, DVec2::ZERO, true);
        let s = solve(&g.sketch, &opts());
        prop_assert!(
            matches!(s.status, SolveStatus::Solved | SolveStatus::Redundant),
            "status {:?}, conflicting {:?}",
            s.status, s.conflicting
        );
        assert_satisfied(&s.sketch);
    }

    /// Perturbed by up to 0.01 of the generator's unit in each coordinate
    /// (0.25–0.5 % of the sketch size, 1.5–10 % of an ellipse's semi-axes),
    /// they solve and every constraint holds.
    ///
    /// Larger perturbations of sketches with ellipses fail more often than
    /// those of the other motifs: random fixes and dimensions often
    /// over-determine an ellipse (four of an arc's five points fixed with
    /// its centre free: three equations on two unknowns), and the
    /// least-squares problem of an over-determined system has local minima
    /// near its solution. Measured (design note 0008 §9): 1 in 300,000 up
    /// to this size and up to 0.02, about 8 in 100,000 up to 0.05, each an
    /// over-determined or near-singular sketch; every failure says so (the
    /// next property).
    #[test]
    fn perturbed_sketches_with_ellipses_solve(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
        amount in 0.0f64..0.01,
    ) {
        let scale = 10f64.powi(e);
        let g = generate_with(seed, motifs, scale, DVec2::new(ox, oy), true);
        let input = perturb(&g.sketch, seed, amount, scale);
        let s = solve(&input, &opts());
        prop_assert!(
            matches!(s.status, SolveStatus::Solved | SolveStatus::Redundant),
            "status {:?}, conflicting {:?}, residual {:e}",
            s.status, s.conflicting, s.max_residual
        );
        prop_assert!(s.conflicting.is_empty());
        assert_satisfied(&s.sketch);
    }

    /// Larger perturbations may fail, but never with a wrong answer.
    #[test]
    fn large_perturbations_with_ellipses_never_lie(
        seed in any::<u64>(),
        motifs in 1usize..9,
        amount in 0.01f64..0.5,
    ) {
        let g = generate_with(seed, motifs, 1.0, DVec2::ZERO, true);
        let input = perturb(&g.sketch, seed, amount, 1.0);
        let s = solve(&input, &opts());
        match &s.status {
            SolveStatus::Solved | SolveStatus::Redundant => assert_satisfied(&s.sketch),
            SolveStatus::Failed(_) => prop_assert_eq!(bits(&s.sketch), bits(&input)),
        }
    }
}
