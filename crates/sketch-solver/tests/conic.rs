//! Conic arcs in the solver (design note 0009): degrees of freedom, the ρ
//! dimension across elliptic, parabolic and hyperbolic values, points on
//! conics checked against Lee's implicit equation, tangency at the ends,
//! derived images, refusals, dragging, collapse, and constructive property
//! tests, all checked with the independent formulas of `common`.

mod common;

use common::gen::{generate_with_conics, perturb};
use common::{
    assert_satisfied, bits, cnc, dist_to_conic, pos, resolved, size, violation, worst, B,
};
use openshape_sketch_solver::{
    measure, resolve_curves, solve, Constraint, DVec2, EntityState, Geometry, ResolvedCurve,
    SolveError, SolveOptions, SolveStatus, TangentKind,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// A conic with its own points: [conic, start, control, end].
fn conic(b: &mut B, p0: DVec2, p1: DVec2, p2: DVec2, rho: f64) -> [u32; 4] {
    let k0 = b.point(p0.x, p0.y);
    let k1 = b.point(p1.x, p1.y);
    let k2 = b.point(p2.x, p2.y);
    let k = b.entity(Geometry::Conic {
        start: k0,
        control: k1,
        end: k2,
        rho,
    });
    [k, k0, k1, k2]
}

/// The conic's point at `t`, written independently of the solver.
fn at(p: [DVec2; 3], rho: f64, t: f64) -> DVec2 {
    let w = rho / (1.0 - rho);
    let (a, b, c) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t) * w, t * t);
    (p[0] * a + p[1] * b + p[2] * c) / (a + b + c)
}

const P: [DVec2; 3] = [
    DVec2::new(0.0, 0.0),
    DVec2::new(0.3, 0.8),
    DVec2::new(1.0, 0.1),
];

#[test]
fn a_conic_has_seven_freedoms_and_six_with_its_rho_dimensioned() {
    let mut b = B::new();
    let [k, ..] = conic(&mut b, P[0], P[1], P[2], 0.4);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 7);
    assert_eq!(s.entities[&k], EntityState::Free);
    b.c(Constraint::Rho {
        curve: k,
        value: 0.4,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 6);
    // The stored ρ already is the dimension's: it comes back bit for bit.
    assert_eq!(cnc(&s.sketch, k).1.to_bits(), 0.4f64.to_bits());
}

#[test]
fn the_rho_dimension_sets_the_shoulder_point() {
    for value in [0.05, 0.2, 0.5, 0.8, 0.97] {
        let mut b = B::new();
        let [k, k0, k1, k2] = conic(&mut b, P[0], P[1], P[2], 0.5);
        for q in [k0, k1, k2] {
            b.fix(q);
        }
        b.c(Constraint::Rho { curve: k, value });
        let s = solve(&b.build(), &opts());
        assert_eq!(s.status, SolveStatus::Solved, "{value}: {:?}", s.status);
        assert_eq!(s.dof, 0);
        assert!(s.entities.values().all(|&e| e == EntityState::Constrained));
        let (p, rho) = cnc(&s.sketch, k);
        assert_eq!(rho, value, "the dimension's value, exactly");
        // The shoulder C(½) is ρ of the way from the chord's midpoint to
        // the control point (Liming; P&T §7.3).
        let m = (p[0] + p[2]) * 0.5;
        let shoulder = at(p, rho, 0.5);
        assert!((shoulder - (m + (p[1] - m) * value)).length() < 1e-14);
        let probe = Constraint::Rho {
            curve: k,
            value: 0.0,
        };
        assert!((measure(&s.sketch, &probe).unwrap() - value).abs() <= 2.0 * f64::EPSILON);
    }
}

#[test]
fn points_on_conics_satisfy_lees_implicit_equation() {
    for rho in [0.2, 0.5, 0.75] {
        let mut b = B::new();
        let [k, ..] = conic(&mut b, P[0], P[1], P[2], rho);
        let q = b.point(0.9, 0.9);
        b.c(Constraint::PointOnCurve { point: q, curve: k });
        let s = solve(&b.build(), &opts());
        assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
        assert_eq!(s.dof, 8, "a point on the curve removes one freedom");
        let d = dist_to_conic(&s.sketch, k, pos(&s.sketch, q));
        assert!(d < 1e-12, "{d:e}");
        assert_satisfied(&s.sketch);
    }
}

#[test]
fn a_point_through_which_the_conic_must_pass_sets_its_rho() {
    // Fixed triangle; the conic of ρ = 0.37 passes through x. Starting
    // from ρ = 0.5 with x fixed on the curve, ρ is found.
    let x = at(P, 0.37, 0.3);
    let mut b = B::new();
    let [k, k0, k1, k2] = conic(&mut b, P[0], P[1], P[2], 0.5);
    for q in [k0, k1, k2] {
        b.fix(q);
    }
    let kx = b.point(x.x, x.y);
    b.fix(kx);
    b.c(Constraint::PointOnCurve {
        point: kx,
        curve: k,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert_eq!(s.dof, 0);
    assert!((cnc(&s.sketch, k).1 - 0.37).abs() < 1e-10);
    assert_eq!(s.entities[&k], EntityState::Constrained);
}

#[test]
fn ends_continue_tangentially_into_lines_and_arcs() {
    let mut b = B::new();
    let [k, k0, _, k2] = conic(&mut b, P[0], P[1], P[2], 0.6);
    // A line off the end, bent away from the tangent; an arc into the
    // start, its centre off the normal.
    let (l, ls, _) = b.line(1.0, 0.1, 1.6, 0.0);
    b.coincident(k2, ls);
    b.c(Constraint::Tangent {
        a: k,
        b: l,
        kind: TangentKind::AtPoint { point: ls },
    });
    let (arc, _, _, ae) = b.arc(DVec2::new(-0.3, 0.2), 0.36, 2.0, -0.6);
    b.coincident(ae, k0);
    b.c(Constraint::Tangent {
        a: arc,
        b: k,
        kind: TangentKind::AtPoint { point: ae },
    });
    let s = solve(&b.build(), &opts());
    assert!(
        matches!(s.status, SolveStatus::Solved),
        "{:?} {:?}",
        s.status,
        s.conflicting
    );
    assert_satisfied(&s.sketch);
    let (p, _) = cnc(&s.sketch, k);
    let (a, e) = common::line(&s.sketch, l);
    assert!(
        (p[2] - p[1])
            .normalize()
            .perp_dot((e - a).normalize())
            .abs()
            < 1e-9
    );
}

#[test]
fn what_conics_do_not_take_is_refused() {
    let invalid = |b: &B, what: &str| {
        let s = solve(&b.build(), &opts());
        match &s.status {
            SolveStatus::Failed(SolveError::Invalid(m)) => {
                assert!(m.contains(what), "expected {what:?} in {m:?}")
            }
            other => panic!("expected Invalid ({what}), got {other:?}"),
        }
    };
    let base = |rho: f64| {
        let mut b = B::new();
        let ks = conic(&mut b, P[0], P[1], P[2], rho);
        (b, ks)
    };
    for rho in [0.0, 1.0, -0.2, 1.5, f64::NAN, 1e-7] {
        invalid(&base(rho).0, "ρ must lie in");
    }
    let (mut b, [k, ..]) = base(0.5);
    b.c(Constraint::Radius {
        curve: k,
        value: 1.0,
    });
    invalid(&b, "must be a circle or an arc");
    let (mut b, _) = base(0.5);
    let (circle, _) = b.circle(2.0, 0.0, 0.5);
    b.c(Constraint::Rho {
        curve: circle,
        value: 0.5,
    });
    invalid(&b, "must be a conic");
    let (mut b, [k2, ..]) = base(0.5);
    b.c(Constraint::Rho {
        curve: k2,
        value: 1.0,
    });
    invalid(&b, "ρ must lie in");
    let (mut b, [k3, _, k1, _]) = base(0.5);
    b.c(Constraint::PointOnCurve {
        point: k1,
        curve: k3,
    });
    invalid(&b, "control point");
    let (mut b, [k4, ..]) = base(0.5);
    let (l, ..) = b.line(0.0, 1.0, 1.0, 1.0);
    b.c(Constraint::Tangent {
        a: l,
        b: k4,
        kind: TangentKind::LineCircle { center_left: true },
    });
    invalid(&b, "away from its ends");
    let (mut b, [k5, ..]) = base(0.5);
    let [k6, ..] = conic(
        &mut b,
        P[0] + DVec2::X,
        P[1] + DVec2::X,
        P[2] + DVec2::X,
        0.5,
    );
    b.c(Constraint::Equal { a: k5, b: k6 });
    invalid(&b, "equal needs");
    // A tangency at the control point is not at an end.
    let (mut b, [k7, _, k1, _]) = base(0.5);
    let l2 = b.point(1.0, 1.0);
    let line = b.line_between(k1, l2);
    b.c(Constraint::Tangent {
        a: k7,
        b: line,
        kind: TangentKind::AtPoint { point: k1 },
    });
    invalid(&b, "not an end point");
    // Repeated points.
    let mut b = B::new();
    let k0 = b.point(0.0, 0.0);
    let k1 = b.point(1.0, 1.0);
    b.entity(Geometry::Conic {
        start: k0,
        control: k1,
        end: k0,
        rho: 0.5,
    });
    invalid(&b, "distinct");
}

#[test]
fn mirrored_and_patterned_conics_add_no_freedom_and_keep_rho() {
    let mut b = B::new();
    let (axis, ..) = b.line(-0.2, -2.0, -0.2, 2.0);
    let [src, ..] = conic(&mut b, P[0], P[1], P[2], 0.3);
    let (dir, ..) = b.line(0.0, -1.5, 0.5, -1.5);
    let free = solve(&b.build(), &opts()).dof;
    let reflect = |p: DVec2| DVec2::new(-0.4 - p.x, p.y);
    // Images a little off, with another ρ: the derived constraints put
    // them right.
    let [img, ..] = conic(
        &mut b,
        reflect(P[0]) + DVec2::new(0.02, 0.0),
        reflect(P[1]),
        reflect(P[2]) - DVec2::new(0.0, 0.03),
        0.6,
    );
    let shift = DVec2::new(0.5, 0.0);
    let [copy, ..] = conic(&mut b, P[0] + shift, P[1] + shift, P[2] + shift, 0.45);
    b.c(Constraint::Mirror {
        axis,
        pairs: vec![(src, img)],
    });
    b.c(Constraint::LinearPattern {
        direction: dir,
        copies: vec![(src, copy, 1)],
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!(s.redundant.is_empty());
    assert_eq!(s.dof, free, "images add no degrees of freedom");
    assert_satisfied(&s.sketch);
    let rho = cnc(&s.sketch, src).1;
    assert_eq!(cnc(&s.sketch, img).1, rho);
    assert_eq!(cnc(&s.sketch, copy).1, rho);
    // The source keeps its ρ: the images moved (they weigh nothing).
    assert!((rho - 0.3).abs() < 1e-12, "{rho}");
    // Dimensioning the source's ρ moves the images' with it.
    b.c(Constraint::Rho {
        curve: src,
        value: 0.7,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert_eq!(s.dof, free - 1);
    for k in [src, img, copy] {
        assert_eq!(cnc(&s.sketch, k).1, 0.7);
    }
    assert_satisfied(&s.sketch);
}

#[test]
fn dragging_the_control_point_keeps_a_free_rho() {
    let mut b = B::new();
    let [k, k0, k1, k2] = conic(&mut b, P[0], P[1], P[2], 0.35);
    b.fix(k0);
    b.fix(k2);
    let q = b.point(0.6, 0.6);
    b.c(Constraint::PointOnCurve { point: q, curve: k });
    let target = DVec2::new(0.6, 1.2);
    let s = solve(
        &b.build(),
        &SolveOptions {
            drag: vec![(k1, target)],
            ..SolveOptions::default()
        },
    );
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert!((pos(&s.sketch, k1) - target).length() < 1e-6);
    assert_satisfied(&s.sketch);
    // Without a point on it, the drag leaves ρ where it was.
    let mut b = B::new();
    let [k, _, k1, _] = conic(&mut b, P[0], P[1], P[2], 0.35);
    let s = solve(
        &b.build(),
        &SolveOptions {
            drag: vec![(k1, target)],
            ..SolveOptions::default()
        },
    );
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(cnc(&s.sketch, k).1, 0.35);
}

#[test]
fn a_conic_forced_flat_is_a_conflict_not_a_result() {
    // The control point pulled onto the chord by a point-on-line.
    let mut b = B::new();
    let [_, k0, k1, k2] = conic(&mut b, P[0], P[1], P[2], 0.5);
    let chord = b.line_between(k0, k2);
    b.fix(k0);
    b.fix(k2);
    let on = b.c(Constraint::PointOnCurve {
        point: k1,
        curve: chord,
    });
    let input = b.build();
    let s = solve(&input, &opts());
    assert_eq!(
        s.status,
        SolveStatus::Failed(SolveError::Conflicting),
        "{:?}",
        s.status
    );
    assert_eq!(s.conflicting, vec![on]);
    assert_eq!(bits(&s.sketch), bits(&input));
    // Its ends brought together.
    let mut b = B::new();
    let [_, k0, _, k2] = conic(&mut b, P[0], P[1], P[2], 0.5);
    let meet = b.coincident(k0, k2);
    let input = b.build();
    let s = solve(&input, &opts());
    assert_eq!(s.status, SolveStatus::Failed(SolveError::Conflicting));
    assert_eq!(s.conflicting, vec![meet]);
    // Made to pass through its control point: ρ would have to reach 1.
    let mut b = B::new();
    let [k, k0, k1, k2] = conic(&mut b, P[0], P[1], P[2], 0.5);
    for q in [k0, k1, k2] {
        b.fix(q);
    }
    let x = b.point(P[1].x, P[1].y);
    b.fix(x);
    b.c(Constraint::PointOnCurve { point: x, curve: k });
    let input = b.build();
    let s = solve(&input, &opts());
    assert!(matches!(s.status, SolveStatus::Failed(_)), "{:?}", s.status);
    assert_eq!(bits(&s.sketch), bits(&input));
}

#[test]
fn a_flat_conic_in_the_input_is_degenerate() {
    let mut b = B::new();
    conic(&mut b, P[0], (P[0] + P[2]) * 0.5, P[2], 0.5);
    let s = solve(&b.build(), &opts());
    match &s.status {
        SolveStatus::Failed(SolveError::Degenerate(m)) => assert!(m.contains("is flat"), "{m}"),
        other => panic!("{other:?}"),
    }
    match resolve_curves(&b.build()) {
        Err(SolveError::Degenerate(m)) => assert!(m.contains("is flat"), "{m}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn conics_resolve_to_their_weight() {
    let mut b = B::new();
    let [k, ..] = conic(&mut b, P[0], P[1], P[2], 0.25);
    match resolved(&b.build(), k) {
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => {
            assert_eq!([start, control, end], P);
            assert!((weight - 1.0 / 3.0).abs() < 1e-15);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn the_checker_measures_off_curve_points() {
    // The independent checker itself: zero on the curve, the offset (to
    // first order) beside it, and large on the other side of the chord.
    let mut b = B::new();
    let [k, ..] = conic(&mut b, P[0], P[1], P[2], 0.4);
    let s = b.build();
    let x = at(P, 0.4, 0.6);
    assert!(dist_to_conic(&s, k, x) < 1e-14);
    let n = (at(P, 0.4, 0.6 + 1e-6) - at(P, 0.4, 0.6 - 1e-6))
        .perp()
        .normalize();
    let d = dist_to_conic(&s, k, x + n * 1e-4);
    assert!((d - 1e-4).abs() < 1e-6, "{d}");
    // The point of the full conic's other arc, beyond the chord.
    assert!(dist_to_conic(&s, k, (P[0] + P[2]) * 0.5 - (P[1] - P[0]) * 0.2) > 0.01);
    let c = Constraint::PointOnCurve { point: 1, curve: k };
    assert!(violation(&s, &c).0 < 1e-14);
}

proptest! {
    #![proptest_config(ProptestConfig {
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// The generator's sketches with conics satisfy their constraints.
    #[test]
    fn generated_sketches_with_conics_hold(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
    ) {
        let scale = 10f64.powi(e);
        let g = generate_with_conics(seed, motifs, scale, DVec2::new(ox, oy));
        let (w, what) = worst(&g.sketch, size(&g.sketch));
        prop_assert!(w <= 1e-11, "generator: violation {:e} at {}", w, what);
    }

    /// Satisfied sketches with conics stay where they are.
    #[test]
    fn satisfied_sketches_with_conics_stay_put(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
    ) {
        let scale = 10f64.powi(e);
        let g = generate_with_conics(seed, motifs, scale, DVec2::ZERO);
        let s = solve(&g.sketch, &opts());
        prop_assert!(
            matches!(s.status, SolveStatus::Solved | SolveStatus::Redundant),
            "status {:?}, conflicting {:?}",
            s.status, s.conflicting
        );
        assert_satisfied(&s.sketch);
    }

    /// Perturbed by up to 0.01 of the generator's unit, they solve and
    /// every constraint holds.
    #[test]
    fn perturbed_sketches_with_conics_solve(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
        amount in 0.0f64..0.01,
    ) {
        let scale = 10f64.powi(e);
        let g = generate_with_conics(seed, motifs, scale, DVec2::new(ox, oy));
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
    fn large_perturbations_with_conics_never_lie(
        seed in any::<u64>(),
        motifs in 1usize..9,
        amount in 0.01f64..0.5,
    ) {
        let g = generate_with_conics(seed, motifs, 1.0, DVec2::ZERO);
        let input = perturb(&g.sketch, seed, amount, 1.0);
        let s = solve(&input, &opts());
        match &s.status {
            SolveStatus::Solved | SolveStatus::Redundant => assert_satisfied(&s.sketch),
            SolveStatus::Failed(_) => prop_assert_eq!(bits(&s.sketch), bits(&input)),
        }
    }
}
