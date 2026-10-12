//! The derived pattern constraints (design note 0007 §11): copies follow
//! their sources exactly, add no degrees of freedom and no redundancy, a
//! circular pattern's angle is a dimension that moves only the copies,
//! copies of a hub point share it, and malformed patterns are refused with
//! a message naming the problem.

mod common;

use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

use common::gen::{generate, perturb, Rng};
use common::{assert_satisfied, pattern_violation, pos, size, B};
use openshape_sketch_solver::{
    solve, Constraint, ConstraintEntry, DVec2, Entity, EntityState, Geometry, Sketch, SolveError,
    SolveOptions, SolveStatus,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// Adds copies `1..count` of each of `sources` to `s`, each point moved by
/// `map` (a point and the copy number) and then by up to `jitter` in each
/// coordinate: same kind, new points (one per source point and copy,
/// shared as the sources share them, except points in `shared`, which
/// every copy uses itself). Keys start after `next`. Returns the copies.
#[allow(clippy::too_many_arguments)]
fn add_copies(
    s: &mut Sketch<u32>,
    sources: &[u32],
    count: u32,
    map: &dyn Fn(DVec2, u32) -> DVec2,
    shared: &[u32],
    jitter: f64,
    rng: &mut Rng,
    next: &mut u32,
) -> Vec<(u32, u32, u32)> {
    let mut copies = Vec::new();
    for k in 1..count {
        let mut points: BTreeMap<u32, u32> = shared.iter().map(|&p| (p, p)).collect();
        for &src in sources {
            let g = s.entity(src).expect("source").geometry.clone();
            for p in g.referenced_points() {
                if let Entry::Vacant(slot) = points.entry(p) {
                    let at = map(s.point(p).unwrap(), k)
                        + DVec2::new(rng.range(-jitter, jitter), rng.range(-jitter, jitter));
                    *next += 1;
                    s.entities.push(Entity {
                        id: *next,
                        geometry: Geometry::Point { at },
                        construction: false,
                    });
                    slot.insert(*next);
                }
            }
            let image = match &g {
                Geometry::Point { at } => Geometry::Point {
                    at: map(*at, k) + DVec2::splat(rng.range(-jitter, jitter)),
                },
                Geometry::Circle { center, radius } => Geometry::Circle {
                    center: points[center],
                    radius: radius * (1.0 + rng.range(0.0, jitter)),
                },
                other => other.map_keys(|p| points[&p]),
            };
            *next += 1;
            s.entities.push(Entity {
                id: *next,
                geometry: image,
                construction: false,
            });
            copies.push((src, *next, k));
        }
    }
    copies
}

fn constrain(s: &mut Sketch<u32>, next: &mut u32, constraint: Constraint<u32>) -> u32 {
    *next += 1;
    s.constraints.push(ConstraintEntry {
        id: *next,
        constraint,
    });
    *next
}

fn max_key(s: &Sketch<u32>) -> u32 {
    s.entities
        .iter()
        .map(|e| e.id)
        .chain(s.constraints.iter().map(|c| c.id))
        .max()
        .unwrap_or(0)
}

fn translate(v: DVec2) -> impl Fn(DVec2, u32) -> DVec2 {
    move |p, k| p + v * f64::from(k)
}

fn rotate(c: DVec2, angle: f64) -> impl Fn(DVec2, u32) -> DVec2 {
    move |p, k| c + DVec2::from_angle(f64::from(k) * angle).rotate(p - c)
}

/// One source of every kind, a fixed direction line from (3, 0) to (3.5,
/// 0.25) and a fixed centre at (−1, −1); the copies are added off their
/// places. Returns the sketch, the sources, the line and the centre.
fn every_kind() -> (Sketch<u32>, Vec<u32>, u32, u32) {
    let mut b = B::new();
    let (direction, d0, d1) = b.line(3.0, 0.0, 3.5, 0.25);
    b.fix(d0);
    b.fix(d1);
    let centre = b.point(-1.0, -1.0);
    b.fix(centre);
    let (line, ..) = b.line(0.2, 1.0, 1.0, 1.6);
    let (arc, ..) = b.arc(DVec2::new(0.5, 2.0), 0.4, 0.3, 2.0);
    let (circle, _) = b.circle(1.5, 2.2, 0.3);
    let point = b.point(-0.4, 0.7);
    let control: Vec<u32> = [(0.0, 3.0), (0.4, 3.4), (0.8, 3.1), (1.2, 3.5)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    let control = b.entity(Geometry::ControlSpline { control, degree: 2 });
    let through: Vec<u32> = [(1.5, 3.0), (1.9, 3.3), (2.3, 3.1)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    let through = b.entity(Geometry::InterpolatedSpline {
        through,
        params: vec![0.0, 0.4, 1.0],
    });
    (
        b.build(),
        vec![line, arc, circle, point, control, through],
        direction,
        centre,
    )
}

/// The sources' freedom in [`every_kind`]: line 4, arc 5, circle 3, point
/// 2, control spline 8, interpolated spline 6.
const SOURCES_DOF: usize = 28;

#[test]
fn linear_copies_of_every_kind_follow_their_sources_exactly() {
    let (mut s, sources, direction, _) = every_kind();
    let mut next = max_key(&s);
    let step = DVec2::new(0.5, 0.25);
    let copies = add_copies(
        &mut s,
        &sources,
        4,
        &translate(step),
        &[],
        0.05,
        &mut Rng(7),
        &mut next,
    );
    let pattern = constrain(
        &mut s,
        &mut next,
        Constraint::LinearPattern {
            direction,
            copies: copies.clone(),
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    assert_eq!(out.dof, SOURCES_DOF);
    assert_satisfied(&out.sketch);
    assert!(pattern_violation(&out.sketch, &copies, translate(step)) < 1e-12);
    for &(src, img, _) in &copies {
        assert_eq!(out.entities[&src], EntityState::Free);
        assert_eq!(out.entities[&img], EntityState::Free);
    }
    // Without the pattern the copies are free geometry of their own.
    let mut alone = s.clone();
    alone.constraints.retain(|c| c.id != pattern);
    assert!(solve(&alone, &opts()).dof > SOURCES_DOF + 3 * 20);
}

#[test]
fn circular_copies_of_every_kind_follow_their_sources_exactly() {
    let (mut s, sources, _, centre) = every_kind();
    let mut next = max_key(&s);
    let angle = 72f64.to_radians();
    let c = DVec2::new(-1.0, -1.0);
    let copies = add_copies(
        &mut s,
        &sources,
        5,
        &rotate(c, angle),
        &[],
        0.05,
        &mut Rng(9),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::CircularPattern {
            centre,
            angle,
            copies: copies.clone(),
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    assert_eq!(out.dof, SOURCES_DOF);
    assert_satisfied(&out.sketch);
    assert!(pattern_violation(&out.sketch, &copies, rotate(c, angle)) < 1e-12);
}

#[test]
fn the_angle_is_a_dimension_that_moves_only_the_copies() {
    let (mut s, sources, _, centre) = every_kind();
    let mut next = max_key(&s);
    let c = DVec2::new(-1.0, -1.0);
    let angle = 0.5;
    let copies = add_copies(
        &mut s,
        &sources,
        3,
        &rotate(c, angle),
        &[],
        0.0,
        &mut Rng(1),
        &mut next,
    );
    let pattern = constrain(
        &mut s,
        &mut next,
        Constraint::CircularPattern {
            centre,
            angle,
            copies: copies.clone(),
        },
    );
    let settled = solve(&s, &opts()).sketch;
    let mut changed = settled.clone();
    for c in changed.constraints.iter_mut().filter(|c| c.id == pattern) {
        assert_eq!(c.constraint.value(), Some(0.5));
        assert!(c.constraint.set_value(-0.75));
    }
    let out = solve(&changed, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(pattern_violation(&out.sketch, &copies, rotate(c, -0.75)) < 1e-12);
    // Every source point is exactly where it was.
    for &src in &sources {
        for p in settled.entity(src).unwrap().geometry.referenced_points() {
            assert_eq!(pos(&out.sketch, p), pos(&settled, p));
        }
    }
}

#[test]
fn the_direction_line_drives_the_copies() {
    let mut b = B::new();
    let (direction, d0, _) = b.line(0.0, -1.0, 1.0, -1.0);
    b.fix(d0);
    b.c(Constraint::Horizontal { line: direction });
    let length = b.c(Constraint::Length {
        line: direction,
        value: 1.0,
    });
    // The source runs along the step, so moving its free end would help
    // the ordinary minimum-norm step (by −6/4 of the line's change).
    let (line, l0, _) = b.line(0.0, 0.0, 0.5, 0.0);
    b.fix(l0);
    b.c(Constraint::Horizontal { line });
    let mut s = b.build();
    let mut next = max_key(&s);
    let copies = add_copies(
        &mut s,
        &[line],
        3,
        &translate(DVec2::X),
        &[],
        0.01,
        &mut Rng(4),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::LinearPattern {
            direction,
            copies: copies.clone(),
        },
    );
    for c in s.constraints.iter_mut().filter(|c| c.id == length) {
        c.constraint.set_value(2.5);
    }
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    // Only the source line's length is free.
    assert_eq!(out.dof, 1);
    assert!(pattern_violation(&out.sketch, &copies, translate(DVec2::X * 2.5)) < 1e-12);
    // The step moved the line's end and the copies, not the source: its
    // free end is where it was.
    let Geometry::Line { end, .. } = out.sketch.entity(line).unwrap().geometry else {
        unreachable!()
    };
    assert!((out.sketch.point(end).unwrap() - DVec2::new(0.5, 0.0)).length() < 1e-12);
}

#[test]
fn spokes_share_the_centre() {
    // A line from the centre out, turned by a quarter three times: every
    // copy starts at the centre point itself.
    let mut b = B::new();
    let centre = b.point(0.0, 0.0);
    let tip = b.point(1.0, 0.2);
    let spoke = b.line_between(centre, tip);
    let mut s = b.build();
    let mut next = max_key(&s);
    let angle = std::f64::consts::FRAC_PI_2;
    let copies = add_copies(
        &mut s,
        &[spoke],
        4,
        &rotate(DVec2::ZERO, angle),
        &[centre],
        0.02,
        &mut Rng(5),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::CircularPattern {
            centre,
            angle,
            copies: copies.clone(),
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    // The centre and the tip: the copies add nothing.
    assert_eq!(out.dof, 4);
    for &(_, img, _) in &copies {
        let Geometry::Line { start, .. } = out.sketch.entity(img).unwrap().geometry else {
            panic!()
        };
        assert_eq!(start, centre);
    }
    let c = pos(&out.sketch, centre);
    assert!(pattern_violation(&out.sketch, &copies, rotate(c, angle)) < 1e-12);
}

#[test]
fn a_dimension_on_a_copy_of_a_dimensioned_source_is_redundant() {
    let (mut s, sources, direction, _) = every_kind();
    let mut next = max_key(&s);
    let line = sources[0];
    let copies = add_copies(
        &mut s,
        &[line],
        2,
        &translate(DVec2::new(0.5, 0.25)),
        &[],
        0.0,
        &mut Rng(2),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::LinearPattern {
            direction,
            copies: copies.clone(),
        },
    );
    constrain(&mut s, &mut next, Constraint::Length { line, value: 1.0 });
    let on_copy = constrain(
        &mut s,
        &mut next,
        Constraint::Length {
            line: copies[0].1,
            value: 1.0,
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Redundant);
    assert_eq!(out.redundant, vec![on_copy]);
}

fn refused(s: &Sketch<u32>, needle: &str) {
    match solve(s, &opts()).status {
        SolveStatus::Failed(SolveError::Invalid(m)) => {
            assert!(m.contains(needle), "message {m:?} lacks {needle:?}")
        }
        other => panic!("expected a refusal mentioning {needle:?}, got {other:?}"),
    }
}

#[test]
fn malformed_patterns_are_refused() {
    let (mut s, sources, direction, centre) = every_kind();
    let mut next = max_key(&s);
    let step = DVec2::new(0.5, 0.25);
    let copies = add_copies(
        &mut s,
        &sources,
        3,
        &translate(step),
        &[],
        0.0,
        &mut Rng(3),
        &mut next,
    );
    let pattern = constrain(
        &mut s,
        &mut next,
        Constraint::LinearPattern {
            direction,
            copies: copies.clone(),
        },
    );
    let with = |c: Constraint<u32>| {
        let mut s = s.clone();
        s.constraints
            .iter_mut()
            .find(|e| e.id == pattern)
            .unwrap()
            .constraint = c;
        s
    };
    let linear =
        |copies: Vec<(u32, u32, u32)>| with(Constraint::LinearPattern { direction, copies });
    let circular = |angle: f64, copies: Vec<(u32, u32, u32)>| {
        with(Constraint::CircularPattern {
            centre,
            angle,
            copies,
        })
    };
    let (line, line_copy, _) = copies[0];
    let (arc, arc_copy, _) = copies[1];
    let line_copy2 = copies.iter().find(|c| c.0 == line && c.2 == 2).unwrap().1;
    refused(&linear(vec![]), "no copies");
    refused(&linear(vec![(line, line_copy, 0)]), "copy 1 or later");
    refused(&linear(vec![(line, line, 1)]), "its own image");
    refused(&linear(vec![(line, arc_copy, 1)]), "not the same kind");
    refused(
        &linear(vec![(line, line_copy, 1), (arc, line_copy, 2)]),
        "not the same kind",
    );
    refused(
        &linear(vec![(line, line_copy, 1), (line, line_copy, 2)]),
        "more than one copy",
    );
    refused(&linear(vec![(direction, line_copy, 1)]), "pattern's line");
    refused(&linear(vec![(line, direction, 1)]), "pattern's line");
    refused(
        &linear(vec![(line, line_copy, 1), (line_copy, line_copy2, 1)]),
        "both",
    );
    let Geometry::Line { end: d_end, .. } = s.entity(direction).unwrap().geometry else {
        unreachable!()
    };
    let point = sources[3];
    refused(
        &linear(vec![(point, d_end, 1)]),
        "end of the pattern's line",
    );
    refused(&circular(0.0, vec![(line, line_copy, 1)]), "non-zero");
    refused(
        &circular(7.0, vec![(line, line_copy, 1)]),
        "less than a turn",
    );
    refused(&circular(f64::NAN, vec![(line, line_copy, 1)]), "finite");
    refused(
        &circular(
            std::f64::consts::PI,
            vec![(line, line_copy, 1), (line, line_copy2, 2)],
        ),
        "whole number of turns",
    );
    refused(&circular(1.0, vec![(point, centre, 1)]), "centre");
    // A point its own copy away from the centre.
    let lone = sources[3];
    let mut s2 = s.clone();
    s2.constraints
        .iter_mut()
        .find(|e| e.id == pattern)
        .unwrap()
        .constraint = Constraint::CircularPattern {
        centre,
        angle: 1.0,
        copies: vec![(lone, lone, 1)],
    };
    refused(&s2, "its own image");
    // A copy sharing a point with its source in a linear pattern (only a
    // circular pattern's centre may be shared).
    let mut s3 = s.clone();
    let mut next3 = max_key(&s3);
    let Geometry::Line { start, .. } = s3.entity(line).unwrap().geometry else {
        unreachable!()
    };
    next3 += 1;
    let far = next3;
    s3.entities.push(Entity {
        id: far,
        geometry: Geometry::Point {
            at: DVec2::new(5.0, 5.0),
        },
        construction: false,
    });
    next3 += 1;
    s3.entities.push(Entity {
        id: next3,
        geometry: Geometry::Line { start, end: far },
        construction: false,
    });
    s3.constraints
        .iter_mut()
        .find(|e| e.id == pattern)
        .unwrap()
        .constraint = Constraint::LinearPattern {
        direction,
        copies: vec![(line, next3, 1)],
    };
    refused(&s3, "cannot be its own copy");
}

#[test]
fn deleting_copies_keeps_the_rest_and_the_line_takes_the_pattern() {
    let c = Constraint::LinearPattern {
        direction: 1u32,
        copies: vec![(2, 3, 1), (2, 4, 2), (5, 6, 1)],
    };
    assert_eq!(
        c.without(|k| k == 4),
        Some(Constraint::LinearPattern {
            direction: 1,
            copies: vec![(2, 3, 1), (5, 6, 1)],
        })
    );
    assert_eq!(
        c.without(|k| k == 2),
        Some(Constraint::LinearPattern {
            direction: 1,
            copies: vec![(5, 6, 1)],
        })
    );
    assert_eq!(c.without(|k| k == 1), None);
    let c = Constraint::CircularPattern {
        centre: 1u32,
        angle: 0.5,
        copies: vec![(2, 3, 1)],
    };
    assert_eq!(c.without(|k| k == 3), None);
    assert_eq!(c.without(|k| k == 1), None);
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Random constrained sketches with random entities patterned along a
    /// random line or about a random point (fixed or free), 2 to 5
    /// instances: the pattern adds no freedom and no redundancy beyond its
    /// line's or centre's own, and everything holds after the solve.
    #[test]
    fn random_patterns_add_no_freedom_and_no_redundancy(
        seed in any::<u64>(),
        motifs in 1usize..6,
        e in -3i32..=2,
        circular in any::<bool>(),
        fixed in any::<bool>(),
        count in 2u32..6,
        jitter in prop_oneof![Just(0.0), 0.0f64..0.2],
        amount in prop_oneof![Just(0.0), 0.0f64..0.02],
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::ZERO);
        let mut s = perturb(&g.sketch, seed, amount, scale);
        let mut rng = Rng(seed ^ 0x5A5A);
        let mut next = max_key(&s);
        let mut point = |s: &mut Sketch<u32>, at: DVec2| {
            next += 1;
            s.entities.push(Entity { id: next, geometry: Geometry::Point { at }, construction: false });
            next
        };
        // The line from a to b, or the centre a.
        let a = rng.pt() * (2.0 * scale);
        let b = a + rng.dir() * (scale * rng.range(0.5, 3.0));
        let a0 = point(&mut s, a);
        let mut held = vec![(a0, a)];
        let mut anchor = a0;
        if !circular {
            let a1 = point(&mut s, b);
            next += 1;
            anchor = next;
            s.entities.push(Entity {
                id: anchor,
                geometry: Geometry::Line { start: a0, end: a1 },
                construction: true,
            });
            held.push((a1, b));
        }
        if fixed {
            for (p, at) in held {
                constrain(&mut s, &mut next, Constraint::Fix { point: p, at });
            }
        }
        let before = solve(&s, &opts());
        prop_assume!(matches!(before.status, SolveStatus::Solved | SolveStatus::Redundant));
        let used: Vec<u32> = s.entities.iter().flat_map(|e| e.geometry.referenced_points()).collect();
        let sources: Vec<u32> = s
            .entities
            .iter()
            .filter(|e| e.id != anchor && e.id != a0)
            .filter(|e| e.geometry.is_curve() || !used.contains(&e.id))
            .map(|e| e.id)
            .filter(|_| rng.chance(0.5))
            .collect();
        prop_assume!(!sources.is_empty());
        let mut s = before.sketch.clone();
        let angle = std::f64::consts::TAU / f64::from(count) * if rng.chance(0.5) { 1.0 } else { -1.0 };
        let map: Box<dyn Fn(DVec2, u32) -> DVec2> = if circular {
            Box::new(rotate(s.point(a0).unwrap(), angle))
        } else {
            let Geometry::Line { start, end } = s.entity(anchor).unwrap().geometry else { unreachable!() };
            Box::new(translate(s.point(end).unwrap() - s.point(start).unwrap()))
        };
        let copies = add_copies(&mut s, &sources, count, &*map, &[], jitter * scale, &mut rng, &mut next);
        constrain(&mut s, &mut next, if circular {
            Constraint::CircularPattern { centre: a0, angle, copies }
        } else {
            Constraint::LinearPattern { direction: anchor, copies }
        });
        let after = solve(&s, &opts());
        prop_assert!(
            matches!(after.status, SolveStatus::Solved | SolveStatus::Redundant),
            "status {:?}, conflicting {:?}", after.status, after.conflicting
        );
        prop_assert_eq!(after.dof, before.dof);
        prop_assert_eq!(&after.redundant, &before.redundant);
        let scale_now = size(&after.sketch).max(1e-300);
        let (w, what) = common::worst(&after.sketch, scale_now);
        prop_assert!(w <= 1e-9, "violation {:e} at {}", w, what);
    }
}

/// A hexagon of six lines joined end to end (radius 0.3 about the origin)
/// and a direction line from (0, −1) to (1, −1) with its start fixed;
/// returns the sketch, the direction line, its end, the hexagon's lines
/// and their (start, end) points.
/// A sketch with its direction line, that line's end, the hexagon's lines
/// and their (start, end) points.
type HexagonAndLine = (Sketch<u32>, u32, u32, Vec<u32>, Vec<(u32, u32)>);

fn hexagon_and_line() -> HexagonAndLine {
    let mut b = B::new();
    let (direction, d0, d1) = b.line(0.0, -1.0, 1.0, -1.0);
    b.fix(d0);
    let pts: Vec<DVec2> = (0..6)
        .map(|i| DVec2::from_angle(f64::from(i) * std::f64::consts::TAU / 6.0) * 0.3)
        .collect();
    let mut lines = Vec::new();
    let mut ends = Vec::new();
    for i in 0..6 {
        let (a, c) = (pts[i], pts[(i + 1) % 6]);
        let (l, s0, s1) = b.line(a.x, a.y, c.x, c.y);
        lines.push(l);
        ends.push((s0, s1));
    }
    for i in 0..6 {
        b.coincident(ends[i].1, ends[(i + 1) % 6].0);
    }
    (b.build(), direction, d1, lines, ends)
}

#[test]
fn copies_only_their_pattern_holds_are_computed_and_follow_their_sources() {
    // A hundred hexagons: the copies are fixed explicitly by the pattern
    // (design note 0005 §5 step 2), so this solves like one hexagon.
    let (mut s, direction, d1, lines, ends) = hexagon_and_line();
    let mut next = max_key(&s);
    let copies = add_copies(
        &mut s,
        &lines,
        100,
        &translate(DVec2::X),
        &[],
        0.01,
        &mut Rng(9),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::LinearPattern {
            direction,
            copies: copies.clone(),
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(out.redundant.is_empty());
    // The hexagon's twelve coordinates and the line's free end.
    assert_eq!(out.dof, 14);
    assert!(pattern_violation(&out.sketch, &copies, translate(DVec2::X)) < 1e-12);
    // Free sources, free copies.
    assert!(copies
        .iter()
        .all(|c| out.entities[&c.1] == EntityState::Free));
    // Fix the hexagon and the line: everything is held, copies included.
    for &(p, _) in &ends {
        let at = s.point(p).unwrap();
        constrain(&mut s, &mut next, Constraint::Fix { point: p, at });
    }
    let at = s.point(d1).unwrap();
    constrain(&mut s, &mut next, Constraint::Fix { point: d1, at });
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert_eq!(out.dof, 0);
    assert!(out
        .entities
        .values()
        .all(|e| *e == EntityState::Constrained));
    assert!(pattern_violation(&out.sketch, &copies, translate(DVec2::X)) < 1e-12);
}

#[test]
fn a_copy_held_by_more_stays_in_the_system_and_a_copy_of_a_copy_follows() {
    let (mut s, direction, _, lines, ends) = hexagon_and_line();
    for &(p, _) in &ends {
        let at = s.point(p).unwrap();
        s.constraints.push(ConstraintEntry {
            id: 10_000 + p,
            constraint: Constraint::Fix { point: p, at },
        });
    }
    let mut next = max_key(&s).max(20_000);
    let copies = add_copies(
        &mut s,
        &lines,
        3,
        &translate(DVec2::X),
        &[],
        0.01,
        &mut Rng(3),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::LinearPattern {
            direction,
            copies: copies.clone(),
        },
    );
    // A circle of copies of copy 2's first line, turned about a free
    // point: copies of a copy.
    let centre = {
        next += 1;
        s.entities.push(Entity {
            id: next,
            geometry: Geometry::Point {
                at: DVec2::new(3.0, 2.0),
            },
            construction: false,
        });
        next
    };
    let of_copy = copies.iter().find(|c| c.2 == 2).unwrap().1;
    let turned = add_copies(
        &mut s,
        &[of_copy],
        4,
        &rotate(DVec2::new(3.0, 2.0), 1.0),
        &[],
        0.01,
        &mut Rng(5),
        &mut next,
    );
    constrain(
        &mut s,
        &mut next,
        Constraint::CircularPattern {
            centre,
            angle: 1.0,
            copies: turned.clone(),
        },
    );
    // A free point held on copy 1's first line's start: that copy point is
    // in two equations and is solved with the rest.
    let held = {
        next += 1;
        s.entities.push(Entity {
            id: next,
            geometry: Geometry::Point {
                at: DVec2::new(0.9, 0.1),
            },
            construction: false,
        });
        next
    };
    let first = copies.iter().find(|c| c.2 == 1).unwrap().1;
    let Geometry::Line { start, .. } = s.entity(first).unwrap().geometry else {
        unreachable!()
    };
    constrain(
        &mut s,
        &mut next,
        Constraint::Coincident { a: held, b: start },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(out.redundant.is_empty());
    // The direction line's free end and the free centre.
    assert_eq!(out.dof, 4);
    assert_satisfied(&out.sketch);
    let step = {
        let Geometry::Line { start, end } = out.sketch.entity(direction).unwrap().geometry else {
            unreachable!()
        };
        out.sketch.point(end).unwrap() - out.sketch.point(start).unwrap()
    };
    assert!(pattern_violation(&out.sketch, &copies, translate(step)) < 1e-12);
    let c = out.sketch.point(centre).unwrap();
    assert!(pattern_violation(&out.sketch, &turned, rotate(c, 1.0)) < 1e-12);
    assert!((out.sketch.point(held).unwrap() - out.sketch.point(start).unwrap()).length() < 1e-12);
    // The copies move with the line: free; so does the held point.
    assert_eq!(out.entities[&held], EntityState::Free);
    assert!(turned
        .iter()
        .all(|t| out.entities[&t.1] == EntityState::Free));
    // Dragging a copy of a copy turns it about the centre (the copies'
    // sources are fixed), moving the centre: the drag reaches through.
    let Geometry::Line { end: tip, .. } = out.sketch.entity(turned[0].1).unwrap().geometry else {
        unreachable!()
    };
    let from = out.sketch.point(tip).unwrap();
    let dragged = solve(
        &out.sketch,
        &SolveOptions {
            drag: vec![(tip, from + DVec2::new(0.05, 0.0))],
            ..opts()
        },
    );
    assert_eq!(dragged.status, SolveStatus::Solved);
    assert!((dragged.sketch.point(centre).unwrap() - c).length() > 1e-3);
    let c2 = dragged.sketch.point(centre).unwrap();
    let step2 = {
        let Geometry::Line { start, end } = dragged.sketch.entity(direction).unwrap().geometry
        else {
            unreachable!()
        };
        dragged.sketch.point(end).unwrap() - dragged.sketch.point(start).unwrap()
    };
    assert!(pattern_violation(&dragged.sketch, &copies, translate(step2)) < 1e-12);
    assert!(pattern_violation(&dragged.sketch, &turned, rotate(c2, 1.0)) < 1e-12);
}
