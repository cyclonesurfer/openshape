//! Branch stability (design note 0005 §6, §8): changing dimensions keeps
//! the sketch on the branch it is on, and undoing a change returns to the
//! start.

mod common;

use common::{assert_satisfied, pos, B};
use opendatum_sketch_solver::{
    solve, Constraint, DVec2, Sketch, SolveOptions, SolveStatus, TangentKind,
};

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// A triangle from three lines with coincident corners: base fixed along
/// the x axis, three side lengths dimensioned. Returns (sketch, the three
/// length constraints' keys, the corner points a, b, c).
fn triangle(apex: DVec2) -> (Sketch<u32>, [u32; 3], [u32; 3]) {
    let mut b = B::new();
    let (base, a0, a1) = b.line(0.0, 0.0, 4.0, 0.0);
    let (right, b0, b1) = b.line(4.0, 0.0, apex.x, apex.y);
    let (left, c0, c1) = b.line(apex.x, apex.y, 0.0, 0.0);
    b.coincident(a1, b0);
    b.coincident(b1, c0);
    b.coincident(c1, a0);
    b.fix(a0);
    b.c(Constraint::Horizontal { line: base });
    let mut dims = [0; 3];
    for (k, l) in [base, right, left].into_iter().enumerate() {
        let (s, e) = common::line(&b.sketch, l);
        dims[k] = b.c(Constraint::Length {
            line: l,
            value: (e - s).length(),
        });
    }
    (b.build(), dims, [a0, b0, c0])
}

fn orientation(s: &Sketch<u32>, [a, b, c]: [u32; 3]) -> f64 {
    let (pa, pb, pc) = (pos(s, a), pos(s, b), pos(s, c));
    (pb - pa).perp_dot(pc - pa)
}

fn set_value(s: &mut Sketch<u32>, key: u32, v: f64) {
    for c in s.constraints.iter_mut() {
        if c.id == key {
            assert!(c.constraint.set_value(v));
        }
    }
}

fn value(s: &Sketch<u32>, key: u32) -> f64 {
    s.constraint(key).unwrap().constraint.value().unwrap()
}

#[test]
fn triangle_never_flips_under_repeated_length_changes() {
    for apex in [
        DVec2::new(1.0, 3.0),
        DVec2::new(2.5, -2.0),
        DVec2::new(3.5, 0.6),
    ] {
        let (mut s, dims, corners) = triangle(apex);
        let sign = orientation(&s, corners).signum();
        // Changes of ±20 % to each side in turn, many times; each change
        // is solved from the previous result (as in an editing session).
        let mut seed = 0x1234_5678_u64;
        for step in 0..60 {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let k = (seed >> 33) as usize % 3;
            let factor = if (seed >> 20) & 1 == 0 {
                1.2
            } else {
                1.0 / 1.2
            };
            let mut next = s.clone();
            let v = value(&next, dims[k]);
            set_value(&mut next, dims[k], v * factor);
            let r = solve(&next, &opts());
            match r.status {
                SolveStatus::Solved => {
                    assert_satisfied(&r.sketch);
                    let o = orientation(&r.sketch, corners);
                    assert_eq!(o.signum(), sign, "flipped at step {step}");
                    s = r.sketch;
                }
                // A change that breaks the triangle inequality cannot be
                // solved; the sketch stays as it was.
                SolveStatus::Failed(_) => {
                    let l: Vec<f64> = dims.iter().map(|&d| value(&next, d)).collect();
                    let longest = l.iter().cloned().fold(0.0, f64::max);
                    assert!(
                        2.0 * longest >= l.iter().sum::<f64>() * (1.0 - 1e-9),
                        "a valid triangle {l:?} failed at step {step}: {:?}",
                        r.status
                    );
                    assert_eq!(r.sketch, next);
                }
                SolveStatus::Redundant => panic!("nothing is redundant here"),
            }
        }
    }
}

#[test]
fn undoing_a_dimension_change_returns_to_the_start() {
    let (s0, dims, corners) = triangle(DVec2::new(1.0, 3.0));
    let start = solve(&s0, &opts());
    assert_eq!(start.status, SolveStatus::Solved);
    for (k, &dim) in dims.iter().enumerate() {
        for factor in [1.2, 0.8, 1.05] {
            let mut changed = start.sketch.clone();
            let v = value(&changed, dim);
            set_value(&mut changed, dim, v * factor);
            let mid = solve(&changed, &opts());
            assert_eq!(mid.status, SolveStatus::Solved);
            let mut back = mid.sketch.clone();
            set_value(&mut back, dim, v);
            let end = solve(&back, &opts());
            assert_eq!(end.status, SolveStatus::Solved);
            // The triangle is fully determined, so the way back ends where
            // the way out started, to rounding.
            for p in corners {
                let d = (pos(&end.sketch, p) - pos(&start.sketch, p)).length();
                assert!(d <= 1e-12, "side {k} ×{factor}: moved {d:e}");
            }
            assert!(orientation(&end.sketch, corners) > 0.0);
        }
    }
}

#[test]
fn large_angle_changes_keep_the_sense_of_the_angle() {
    // Two lines from a fixed common point; the angle from the first to the
    // second is driven across 0 and through ±π/2 towards ±π, in steps and
    // in one jump; the signed angle always matches the dimension.
    let mut b = B::new();
    let (la, a0, _) = b.line(0.0, 0.0, 2.0, 0.0);
    let (lb, b0, _) = b.line(0.0, 0.0, 1.0, 1.0);
    b.coincident(a0, b0);
    b.fix(a0);
    b.c(Constraint::Horizontal { line: la });
    b.c(Constraint::Length {
        line: la,
        value: 2.0,
    });
    b.c(Constraint::Length {
        line: lb,
        value: 2f64.sqrt(),
    });
    let angle = b.c(Constraint::Angle {
        a: la,
        b: lb,
        value: std::f64::consts::FRAC_PI_4,
    });
    let base = b.build();
    let signed = |s: &Sketch<u32>| {
        let (a0, a1) = common::line(s, la);
        let (b0, b1) = common::line(s, lb);
        let (da, db) = (a1 - a0, b1 - b0);
        da.perp_dot(db).atan2(da.dot(db))
    };
    for target in [3.0, -3.0, 2.0, -0.2, 0.0, 1.5, -1.5, std::f64::consts::PI] {
        let mut s = base.clone();
        set_value(&mut s, angle, target);
        let r = solve(&s, &opts());
        assert_eq!(r.status, SolveStatus::Solved, "target {target}");
        let got = signed(&r.sketch);
        let diff = (got - target + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI;
        assert!(diff.abs() < 1e-9, "target {target}, got {got}");
    }
    // Stepwise, each from the previous result.
    let mut s = base;
    for i in 0..=24 {
        let target = std::f64::consts::FRAC_PI_4 - i as f64 * 0.15;
        set_value(&mut s, angle, target);
        let r = solve(&s, &opts());
        assert_eq!(r.status, SolveStatus::Solved, "target {target}");
        assert!((signed(&r.sketch) - target).abs() < 1e-9);
        s = r.sketch;
    }
}

#[test]
fn an_arc_keeps_its_side_when_its_radius_changes() {
    // A fillet-like arc between two perpendicular lines; growing and
    // shrinking its radius keeps the centre inside the corner.
    let mut b = B::new();
    let (h, h0, h1) = b.line(-2.0, 0.0, -0.5, 0.0);
    let (v, v0, v1) = b.line(0.0, 0.5, 0.0, 2.0);
    let (arc, c, s, e) = b.arc(
        DVec2::new(-0.5, 0.5),
        0.5,
        -std::f64::consts::FRAC_PI_2,
        0.0,
    );
    b.coincident(h1, s);
    b.coincident(e, v0);
    b.fix(h0);
    b.fix(v1);
    b.c(Constraint::Horizontal { line: h });
    b.c(Constraint::Vertical { line: v });
    b.c(Constraint::Tangent {
        a: h,
        b: arc,
        kind: TangentKind::AtPoint { point: h1 },
    });
    b.c(Constraint::Tangent {
        a: arc,
        b: v,
        kind: TangentKind::AtPoint { point: v0 },
    });
    let rad = b.c(Constraint::Radius {
        curve: arc,
        value: 0.5,
    });
    let mut sk = b.build();
    for r in [0.6, 0.8, 1.0, 0.7, 0.3, 0.1, 0.5] {
        set_value(&mut sk, rad, r);
        let out = solve(&sk, &opts());
        assert!(
            matches!(out.status, SolveStatus::Solved | SolveStatus::Redundant),
            "radius {r}: {:?}",
            out.status
        );
        assert_satisfied(&out.sketch);
        let cc = pos(&out.sketch, c);
        assert!(
            (cc - DVec2::new(-r, r)).length() < 1e-9,
            "radius {r}: centre {cc}"
        );
        sk = out.sketch;
    }
}
