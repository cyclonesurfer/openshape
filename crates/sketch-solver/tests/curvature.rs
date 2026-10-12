//! Curvature continuity (G2) in the solver (design note 0010): a spline
//! joined to every kind of curve it accepts, at either end, checked against
//! the independent end curvatures of `common`; degrees of freedom;
//! redundancy with a tangency; conflicts; refusals; scale invariance; and a
//! property test over random joints.

mod common;

use common::{assert_satisfied, away_curvature, pos, size, violation, B};
use opendatum_sketch_solver::{
    solve, Constraint, DVec2, Geometry, Sketch, SolveError, SolveOptions, SolveStatus, TangentKind,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// A control-point spline on its own points: (spline, points).
fn control(b: &mut B, pts: &[DVec2], degree: u32) -> (u32, Vec<u32>) {
    let control: Vec<u32> = pts.iter().map(|p| b.point(p.x, p.y)).collect();
    let k = b.entity(Geometry::ControlSpline {
        control: control.clone(),
        degree,
        periodic: false,
    });
    (k, control)
}

/// An interpolated spline through its own points, parameters by chord
/// length: (spline, points).
fn through(b: &mut B, pts: &[DVec2]) -> (u32, Vec<u32>) {
    let keys: Vec<u32> = pts.iter().map(|p| b.point(p.x, p.y)).collect();
    let mut acc = vec![0.0];
    for w in pts.windows(2) {
        acc.push(acc.last().unwrap() + (w[1] - w[0]).length());
    }
    let total = *acc.last().unwrap();
    let mut params: Vec<f64> = acc.iter().map(|a| a / total).collect();
    *params.last_mut().unwrap() = 1.0;
    let k = b.entity(Geometry::InterpolatedSpline {
        through: keys.clone(),
        params,
        periodic: false,
    });
    (k, keys)
}

/// The curves a spline can be joined to, and which of their ends.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Other {
    Line,
    ArcStart,
    ArcEnd,
    EllipseStart,
    EllipseEnd,
    /// An elliptical arc whose minor axis is right of its major axis.
    MirroredEllipseEnd,
    ConicStart,
    ConicEnd,
    Quadratic,
    Cubic,
    Through,
}

const OTHERS: [Other; 11] = [
    Other::Line,
    Other::ArcStart,
    Other::ArcEnd,
    Other::EllipseStart,
    Other::EllipseEnd,
    Other::MirroredEllipseEnd,
    Other::ConicStart,
    Other::ConicEnd,
    Other::Quadratic,
    Other::Cubic,
    Other::Through,
];

/// Builds the other curve so that the end the joint takes is near the
/// origin and the curve lies left of x = 0 (scaled by `k`, shifted by
/// `o`): (curve, the point at the joint).
fn other(b: &mut B, kind: Other, k: f64, o: DVec2) -> (u32, u32) {
    let at = |x: f64, y: f64| o + DVec2::new(x, y) * k;
    match kind {
        Other::Line => {
            let (l, _, e) = {
                let (p, q) = (at(-1.0, -0.1), at(0.0, 0.0));
                b.line(p.x, p.y, q.x, q.y)
            };
            (l, e)
        }
        Other::ArcStart | Other::ArcEnd => {
            // Arcs run counter-clockwise: one starting at the top of its
            // circle heads left, one ending at the bottom arrives from the
            // left.
            let (c, a0, a1) = if kind == Other::ArcStart {
                (at(0.0, -0.8), 1.5, 2.3)
            } else {
                (at(0.0, 0.8), -2.3, -1.5)
            };
            let (arc, _, s, e) = b.arc(c, 0.8 * k, a0, a1);
            (arc, if kind == Other::ArcStart { s } else { e })
        }
        Other::EllipseStart | Other::EllipseEnd | Other::MirroredEllipseEnd => {
            // (centre, first semi-axis, second semi-axis as a multiple of
            // the first turned left, start and end anomalies, joint at the
            // start). With the second axis on the right (a negative
            // multiple) increasing anomaly runs clockwise, so the arc runs
            // from the larger anomaly to the smaller.
            let (c, a, ratio, t0, t1, at_start) = match kind {
                Other::EllipseStart => {
                    (at(0.0, -0.6), DVec2::new(0.9, 0.0) * k, 0.6, 1.5, 2.3, true)
                }
                Other::EllipseEnd => (
                    at(0.0, 0.6),
                    DVec2::new(0.9, 0.1) * k,
                    0.6,
                    -2.3,
                    -1.5,
                    false,
                ),
                _ => (
                    at(0.0, 0.6),
                    DVec2::new(0.9, 0.1) * k,
                    -0.6,
                    2.3,
                    1.5,
                    false,
                ),
            };
            let bb = a.perp() * ratio;
            let e_at = |t: f64| c + a * t.cos() + bb * t.sin();
            let kc = b.point(c.x, c.y);
            let km = b.point(c.x + a.x, c.y + a.y);
            let kn = b.point(c.x + bb.x, c.y + bb.y);
            let (p0, p1) = (e_at(t0), e_at(t1));
            let ks = b.point(p0.x, p0.y);
            let ke = b.point(p1.x, p1.y);
            let e = b.entity(Geometry::EllipticalArc {
                center: kc,
                major: km,
                minor: kn,
                start: ks,
                end: ke,
            });
            (e, if at_start { ks } else { ke })
        }
        Other::ConicStart | Other::ConicEnd => {
            let (p0, p1, p2) = (at(-1.0, -0.6), at(-0.6, 0.1), at(0.0, 0.0));
            let (p0, p2) = if kind == Other::ConicStart {
                (p2, p0)
            } else {
                (p0, p2)
            };
            let k0 = b.point(p0.x, p0.y);
            let k1 = b.point(p1.x, p1.y);
            let k2 = b.point(p2.x, p2.y);
            let c = b.entity(Geometry::Conic {
                start: k0,
                control: k1,
                end: k2,
                rho: 0.45,
            });
            (c, if kind == Other::ConicStart { k0 } else { k2 })
        }
        Other::Quadratic | Other::Cubic => {
            let pts = [at(-1.2, -0.3), at(-0.8, 0.2), at(-0.4, -0.1), at(0.0, 0.0)];
            let (s, keys) = control(b, &pts, if kind == Other::Quadratic { 2 } else { 3 });
            (s, keys[3])
        }
        Other::Through => {
            let pts = [
                at(-1.2, -0.3),
                at(-0.7, 0.15),
                at(-0.3, -0.05),
                at(0.0, 0.0),
            ];
            let (s, keys) = through(b, &pts);
            (s, keys[3])
        }
    }
}

/// A cubic spline leaving the origin to the right, roughly (but not) G2
/// with whatever it is joined to: (spline, its first point).
fn leaving(b: &mut B, k: f64, o: DVec2, interpolated: bool) -> (u32, u32) {
    let at = |x: f64, y: f64| o + DVec2::new(x, y) * k;
    let pts = [
        at(0.02, 0.01),
        at(0.4, 0.15),
        at(0.8, -0.1),
        at(1.2, 0.2),
        at(1.6, 0.0),
    ];
    let (s, keys) = if interpolated {
        through(b, &pts)
    } else {
        control(b, &pts, 3)
    };
    (s, keys[0])
}

/// The two-curve sketch joined at a coincident pair, with the curvature
/// constraint or not: (builder, spline, other, the spline's joint point).
fn joint(kind: Other, k: f64, o: DVec2, interpolated: bool, g2: bool) -> (B, u32, u32, u32) {
    let mut b = B::new();
    let (c, q) = other(&mut b, kind, k, o);
    let (s, p) = leaving(&mut b, k, o, interpolated);
    b.coincident(p, q);
    if g2 {
        b.c(Constraint::Curvature {
            a: s,
            b: c,
            point: p,
        });
    }
    (b, s, c, p)
}

#[test]
fn a_spline_joins_every_kind_of_curve_curvature_continuously() {
    for interpolated in [false, true] {
        for kind in OTHERS {
            let name = format!("{kind:?} (interpolated: {interpolated})");
            let (b, s, c, p) = joint(kind, 1.0, DVec2::ZERO, interpolated, true);
            let sol = solve(&b.build(), &opts());
            assert_eq!(sol.status, SolveStatus::Solved, "{name}: {:?}", sol.status);
            assert_satisfied(&sol.sketch);
            // The independent check, at both curves' own formulas: one
            // tangent line, opposite away tangents (a smooth joint, the
            // branch nearest the start), equal curvature through it.
            let (ts, ks) = away_curvature(&sol.sketch, s, p);
            let (tc, kc) = away_curvature(&sol.sketch, c, p);
            assert!(ts.dot(tc) < -1.0 + 1e-12, "{name}: not smooth: {ts} {tc}");
            assert!(
                (ks + kc).abs() <= 1e-9 * (1.0 + ks.abs()),
                "{name}: {ks} vs {kc}"
            );
            // Before solving it was not G2 (the test means something).
            let (ts0, ks0) = away_curvature(&b.sketch, s, p);
            let (tc0, kc0) = away_curvature(&b.sketch, c, p);
            // The spline adapts; a curve that is not a spline keeps its
            // curvature at the joint (design note 0010 §4.1; measured:
            // arcs within 1 %, conics 4 %, elliptical arcs 4.9 %), even
            // where the spline must turn its own the other way.
            if !matches!(
                kind,
                Other::Line | Other::Quadratic | Other::Cubic | Other::Through
            ) {
                assert!((kc / kc0 - 1.0).abs() <= 0.06, "{name}: {kc0} became {kc}");
            }
            assert!(
                ts0.perp_dot(tc0).abs() > 1e-3 || (ks0 + kc0).abs() > 1e-3,
                "{name}: already G2"
            );
        }
    }
}

#[test]
fn curvature_continuity_takes_two_freedoms() {
    for kind in OTHERS {
        let (b, ..) = joint(kind, 1.0, DVec2::ZERO, false, false);
        let free = solve(&b.build(), &opts());
        let (b, ..) = joint(kind, 1.0, DVec2::ZERO, false, true);
        let g2 = solve(&b.build(), &opts());
        assert_eq!(g2.status, SolveStatus::Solved, "{kind:?}");
        assert_eq!(free.dof, g2.dof + 2, "{kind:?}");
    }
}

#[test]
fn a_tangency_on_the_same_joint_is_redundant() {
    let (mut b, s, c, p) = joint(Other::ArcEnd, 1.0, DVec2::ZERO, false, true);
    let t = b.c(Constraint::Tangent {
        a: c,
        b: s,
        kind: TangentKind::AtPoint { point: p },
    });
    let sol = solve(&b.build(), &opts());
    assert_eq!(sol.status, SolveStatus::Redundant);
    // One of the two is implied by the other.
    let g2 = b.sketch.constraints[1].id;
    assert_eq!(sol.redundant.len(), 1, "{:?}", sol.redundant);
    assert!(sol.redundant.contains(&t) || sol.redundant.contains(&g2));
    assert_satisfied(&sol.sketch);
}

#[test]
fn a_curvature_the_rest_forbids_is_a_conflict_not_a_result() {
    // A fixed spline and an arc of fixed radius: the end curvature cannot
    // be both.
    let (mut b, s, c, _) = joint(Other::ArcEnd, 1.0, DVec2::ZERO, false, true);
    let pts = b.sketch.entity(s).unwrap().geometry.referenced_points();
    for q in pts {
        b.fix(q);
    }
    let r = 0.8;
    b.c(Constraint::Radius { curve: c, value: r });
    let sol = solve(&b.build(), &opts());
    match sol.status {
        SolveStatus::Failed(SolveError::Conflicting) => {}
        other => panic!("expected a conflict, got {other:?}"),
    }
    // The input comes back untouched.
    assert_eq!(sol.sketch, b.build());
}

#[test]
fn what_curvature_continuity_does_not_take_is_refused() {
    let invalid = |s: Sketch<u32>, needle: &str| {
        let sol = solve(&s, &opts());
        match &sol.status {
            SolveStatus::Failed(SolveError::Invalid(m)) => {
                assert!(m.contains(needle), "{m}")
            }
            other => panic!("expected '{needle}', got {other:?}"),
        }
    };
    // No spline: a line and an arc.
    let mut b = B::new();
    let (l, _, e) = b.line(-1.0, 0.0, 0.0, 0.0);
    let (a, _, s, _) = b.arc(
        DVec2::new(0.0, 1.0),
        1.0,
        -std::f64::consts::FRAC_PI_2,
        -1.0,
    );
    b.coincident(e, s);
    b.c(Constraint::Curvature {
        a: l,
        b: a,
        point: e,
    });
    invalid(b.build(), "needs a spline");
    // A degree-1 spline.
    let mut b = B::new();
    let (l, _, e) = b.line(-1.0, 0.0, 0.0, 0.0);
    let (sp, keys) = control(
        &mut b,
        &[DVec2::ZERO, DVec2::new(0.5, 0.2), DVec2::new(1.0, 0.0)],
        1,
    );
    b.coincident(e, keys[0]);
    b.c(Constraint::Curvature {
        a: sp,
        b: l,
        point: keys[0],
    });
    invalid(b.build(), "degree 1");
    // A joint away from the spline's ends.
    let mut b = B::new();
    let (l, _, e) = b.line(-1.0, 0.0, 0.0, 0.0);
    let (sp, keys) = control(
        &mut b,
        &[
            DVec2::new(-0.5, 0.5),
            DVec2::ZERO,
            DVec2::new(0.5, 0.2),
            DVec2::new(1.0, 0.0),
        ],
        3,
    );
    b.c(Constraint::Curvature {
        a: sp,
        b: l,
        point: keys[1],
    });
    let _ = e;
    invalid(b.build(), "not an end point");
    // A circle has no ends.
    let mut b = B::new();
    let (ci, cc) = b.circle(0.0, 1.0, 1.0);
    let (sp, keys) = control(
        &mut b,
        &[DVec2::ZERO, DVec2::new(0.5, 0.2), DVec2::new(1.0, 0.0)],
        2,
    );
    b.c(Constraint::Curvature {
        a: sp,
        b: ci,
        point: keys[0],
    });
    let _ = cc;
    invalid(b.build(), "circle");
    // The same curve twice.
    let mut b = B::new();
    let (sp, keys) = control(
        &mut b,
        &[DVec2::ZERO, DVec2::new(0.5, 0.2), DVec2::new(1.0, 0.0)],
        2,
    );
    b.c(Constraint::Curvature {
        a: sp,
        b: sp,
        point: keys[0],
    });
    let sol = solve(&b.build(), &opts());
    assert!(matches!(
        sol.status,
        SolveStatus::Failed(SolveError::Invalid(_))
    ));
}

#[test]
fn joints_solve_alike_at_any_scale() {
    for kind in OTHERS {
        let (b, s, c, p) = joint(kind, 1.0, DVec2::ZERO, false, true);
        let unit = solve(&b.build(), &opts());
        for (k, o) in [
            (1e-4, DVec2::new(3e-4, -2e-4)),
            (1e3, DVec2::new(-500.0, 2000.0)),
        ] {
            let (b2, ..) = joint(kind, k, o, false, true);
            let sol = solve(&b2.build(), &opts());
            assert_eq!(sol.status, SolveStatus::Solved, "{kind:?} at {k}");
            assert_satisfied(&sol.sketch);
            // Same shape: every point maps to the unit solution's.
            for e in &unit.sketch.entities {
                if let Geometry::Point { at } = e.geometry {
                    let q = (pos(&sol.sketch, e.id) - o) / k;
                    assert!((q - at).length() < 1e-9, "{kind:?} at {k}: {q} vs {at}");
                }
            }
            let _ = (s, c, p);
        }
    }
}

#[test]
fn a_cusp_is_never_the_answer() {
    // The spline leaves the joint back along the line it meets, a cusp:
    // the solve either turns it round into a smooth joint or fails; a
    // cusp with matching curvature is never accepted.
    for interpolated in [false, true] {
        let mut b = B::new();
        let (l, _, e) = b.line(-1.0, -0.1, 0.0, 0.0);
        let pts = [
            DVec2::new(0.0, 0.0),
            DVec2::new(-0.4, 0.05),
            DVec2::new(-0.8, 0.3),
            DVec2::new(-1.2, 0.2),
        ];
        let (sp, keys) = if interpolated {
            through(&mut b, &pts)
        } else {
            control(&mut b, &pts, 3)
        };
        b.coincident(e, keys[0]);
        b.c(Constraint::Curvature {
            a: sp,
            b: l,
            point: keys[0],
        });
        let sol = solve(&b.build(), &opts());
        if let SolveStatus::Solved = sol.status {
            assert_satisfied(&sol.sketch);
            let (ts, _) = away_curvature(&sol.sketch, sp, keys[0]);
            let (tl, _) = away_curvature(&sol.sketch, l, keys[0]);
            assert!(ts.dot(tl) < 0.0, "a cusp: {ts} {tl}");
        } else {
            assert!(
                matches!(sol.status, SolveStatus::Failed(_)),
                "{:?}",
                sol.status
            );
        }
    }
}

/// The curvature constraint's own violation, relative to the sketch size.
fn g2_violation(s: &Sketch<u32>, c: &Constraint<u32>) -> f64 {
    violation(s, c).0 / size(s)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// Random joints, perturbed: solved with every constraint holding (the
    /// curvature one by the independent formulas), or reported as failed;
    /// never a wrong result.
    #[test]
    fn random_joints_are_continuous_or_reported(
        kind in 0usize..OTHERS.len(),
        interpolated in any::<bool>(),
        k in 0.05f64..20.0,
        ox in -5.0f64..5.0,
        oy in -5.0f64..5.0,
        jiggle in prop::collection::vec(-0.15f64..0.15, 32),
    ) {
        let kind = OTHERS[kind];
        let o = DVec2::new(ox, oy);
        let (mut b, ..) = joint(kind, k, o, interpolated, true);
        // Perturb every point by up to 15% of the scale.
        let mut j = jiggle.iter().cycle();
        for e in &mut b.sketch.entities {
            if let Geometry::Point { at } = &mut e.geometry {
                *at += DVec2::new(*j.next().unwrap(), *j.next().unwrap()) * k;
            }
        }
        let sol = solve(&b.build(), &opts());
        match sol.status {
            SolveStatus::Solved | SolveStatus::Redundant => {
                assert_satisfied(&sol.sketch);
                for c in &sol.sketch.constraints {
                    if let Constraint::Curvature { .. } = c.constraint {
                        let v = g2_violation(&sol.sketch, &c.constraint);
                        prop_assert!(v <= 1e-9, "{kind:?}: {v:e}");
                    }
                }
            }
            SolveStatus::Failed(_) => {
                prop_assert_eq!(&sol.sketch, &b.build());
            }
        }
    }
}

#[test]
fn corners_up_to_nearly_reversed_solve_smooth() {
    // A spline leaving a line at a corner of 10° to 179° (measured from
    // straight on): solved to a smooth joint, the spline turned round.
    for deg in [10, 30, 60, 90, 120, 150, 170, 179] {
        for interpolated in [false, true] {
            let a = (deg as f64).to_radians();
            let mut b = B::new();
            let (l, _, e) = b.line(-1.0, 0.0, 0.0, 0.0);
            let d = DVec2::new(a.cos(), a.sin());
            let n = d.perp();
            let pts = [
                DVec2::ZERO,
                d * 0.4 + n * 0.05,
                d * 0.8 - n * 0.1,
                d * 1.2 + n * 0.1,
            ];
            let (sp, keys) = if interpolated {
                through(&mut b, &pts)
            } else {
                control(&mut b, &pts, 3)
            };
            b.coincident(e, keys[0]);
            b.c(Constraint::Curvature {
                a: sp,
                b: l,
                point: keys[0],
            });
            let sol = solve(&b.build(), &opts());
            let name = format!("{deg}° (interpolated: {interpolated})");
            assert_eq!(sol.status, SolveStatus::Solved, "{name}");
            assert_satisfied(&sol.sketch);
            let (ts, ks) = away_curvature(&sol.sketch, sp, keys[0]);
            let (tl, _) = away_curvature(&sol.sketch, l, keys[0]);
            assert!(ts.dot(tl) < -1.0 + 1e-12, "{name}: {ts} {tl}");
            assert!(ks.abs() <= 1e-9, "{name}: {ks}");
        }
    }
}

/// A small deterministic generator for the measured statistics below
/// (uniform in [−0.15, 0.15)).
fn jiggle(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64 / (1u64 << 53) as f64) * 0.3 - 0.15
}

#[test]
fn random_joints_mostly_solve_and_keep_the_other_curve() {
    // Every kind of joint, every point moved by up to 15 % of the sketch
    // size, solved with and without the constraint. Measured over 1,000
    // joints of each kind (design note 0010 §6): 30 of 11,000 fail, all
    // reported; where the other curve is not a spline its curvature at
    // the joint ends within 25 % of where it ends without the constraint
    // in 98.3 % of joints (the rest are joints the spline could not reach
    // alone, where every curve gives). Here 100 of each, with bounds that
    // catch a regression.
    let mut seed = 777u64;
    let (mut failed, mut kept, mut compared) = (0, 0, 0);
    for kind in OTHERS {
        for i in 0..100 {
            let interpolated = i % 2 == 0;
            let deltas: Vec<DVec2> = (0..12)
                .map(|_| DVec2::new(jiggle(&mut seed), jiggle(&mut seed)))
                .collect();
            let mut out = Vec::new();
            for g2 in [false, true] {
                let (mut b, _, c, p) = joint(kind, 1.0, DVec2::ZERO, interpolated, g2);
                let mut d = deltas.iter();
                for e in &mut b.sketch.entities {
                    if let Geometry::Point { at } = &mut e.geometry {
                        *at += *d.next().unwrap();
                    }
                }
                let start = b.build();
                let sol = solve(&start, &opts());
                match sol.status {
                    SolveStatus::Failed(_) => assert_eq!(sol.sketch, start),
                    _ => assert_satisfied(&sol.sketch),
                }
                out.push((sol, c, p));
            }
            let (free, c, p) = &out[0];
            let g2 = &out[1].0;
            if matches!(g2.status, SolveStatus::Failed(_)) {
                failed += 1;
                continue;
            }
            if matches!(
                kind,
                Other::Line | Other::Quadratic | Other::Cubic | Other::Through
            ) || matches!(free.status, SolveStatus::Failed(_))
            {
                continue;
            }
            let (_, k0) = away_curvature(&free.sketch, *c, *p);
            let (_, k1) = away_curvature(&g2.sketch, *c, *p);
            compared += 1;
            if (k1 / k0 - 1.0).abs() <= 0.25 {
                kept += 1;
            }
        }
    }
    assert!(failed <= 11, "{failed} of 1100 failed");
    assert!(kept * 100 >= compared * 95, "kept {kept} of {compared}");
}
