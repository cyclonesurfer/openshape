//! The derived mirror constraint (design note 0007 §9): images follow their
//! sources exactly, add no degrees of freedom and no redundancy, a dimension
//! on an image of a dimensioned curve is redundant, and malformed mirrors
//! are refused with a message naming the problem.

mod common;

use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

use common::gen::{generate, perturb, Rng};
use common::{assert_satisfied, size, B};
use openshape_sketch_solver::{
    solve, Constraint, ConstraintEntry, DVec2, Entity, EntityState, Geometry, Sketch, SolveError,
    SolveOptions, SolveStatus,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// The reflection of `p` about the line through `a` and `b`.
fn reflect(p: DVec2, a: DVec2, b: DVec2) -> DVec2 {
    let d = (b - a).normalize();
    let v = p - a;
    a + d * (2.0 * v.dot(d)) - v
}

/// Adds an image of each of `sources` to `s`: same kind, new points (one
/// per source point, shared as the sources share them, except points in
/// `shared`, which the image uses itself), placed at the reflection about
/// the line `a`–`b` moved by up to `jitter` in each coordinate. Keys start
/// at `next`. Returns the `(source, image)` pairs.
fn add_images(
    s: &mut Sketch<u32>,
    sources: &[u32],
    (a, b): (DVec2, DVec2),
    shared: &[u32],
    jitter: f64,
    rng: &mut Rng,
    next: &mut u32,
) -> Vec<(u32, u32)> {
    let mut map: BTreeMap<u32, u32> = shared.iter().map(|&p| (p, p)).collect();
    let mut key = || {
        *next += 1;
        *next
    };
    let mut pairs = Vec::new();
    for &src in sources {
        let g = s.entity(src).expect("source").geometry.clone();
        for p in g.referenced_points() {
            if let Entry::Vacant(slot) = map.entry(p) {
                let at = reflect(s.point(p).unwrap(), a, b)
                    + DVec2::new(rng.range(-jitter, jitter), rng.range(-jitter, jitter));
                let q = key();
                s.entities.push(Entity {
                    id: q,
                    geometry: Geometry::Point { at },
                    construction: false,
                });
                slot.insert(q);
            }
        }
        let image = match &g {
            Geometry::Point { at } => {
                let at = reflect(*at, a, b) + DVec2::splat(rng.range(-jitter, jitter));
                Geometry::Point { at }
            }
            Geometry::Arc { center, start, end } => Geometry::Arc {
                center: map[center],
                start: map[end],
                end: map[start],
            },
            Geometry::Circle { center, radius } => Geometry::Circle {
                center: map[center],
                radius: radius * (1.0 + rng.range(0.0, jitter)),
            },
            other => other.map_keys(|k| map[&k]),
        };
        let id = key();
        s.entities.push(Entity {
            id,
            geometry: image,
            construction: false,
        });
        pairs.push((src, id));
    }
    pairs
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

/// A fixed axis from (0, 0) to (2, 1) and one source of every kind; the
/// images are added off their places. Returns the sketch, the axis and the
/// pairs.
fn every_kind() -> (Sketch<u32>, u32, Vec<(u32, u32)>, u32) {
    let mut b = B::new();
    let (axis, a0, a1) = b.line(0.0, 0.0, 2.0, 1.0);
    b.fix(a0);
    b.fix(a1);
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
    let mut s = b.build();
    let mut next = max_key(&s);
    let mut rng = Rng(7);
    let axis_ends = (DVec2::ZERO, DVec2::new(2.0, 1.0));
    let pairs = add_images(
        &mut s,
        &[line, arc, circle, point, control, through],
        axis_ends,
        &[],
        0.05,
        &mut rng,
        &mut next,
    );
    let mirror = constrain(
        &mut s,
        &mut next,
        Constraint::Mirror {
            axis,
            pairs: pairs.clone(),
        },
    );
    (s, axis, pairs, mirror)
}

#[test]
fn images_of_every_kind_follow_their_sources_exactly() {
    let (s, _, pairs, _) = every_kind();
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    // The sources' freedom only: line 4, arc 5, circle 3, point 2, control
    // spline 8, interpolated spline 6.
    assert_eq!(out.dof, 28);
    assert_satisfied(&out.sketch);
    for &(src, img) in &pairs {
        assert_eq!(out.entities[&src], EntityState::Free);
        assert_eq!(out.entities[&img], EntityState::Free);
    }
    // The image circle's radius is its source's, exactly by substitution.
    let radius = |k: u32| match out.sketch.entity(k).unwrap().geometry {
        Geometry::Circle { radius, .. } => radius,
        _ => unreachable!(),
    };
    let (circle, circle_image) = pairs[2];
    assert_eq!(radius(circle), radius(circle_image));
}

#[test]
fn an_image_arc_runs_counter_clockwise_like_its_source() {
    let (s, _, pairs, _) = every_kind();
    let out = solve(&s, &opts());
    let curves = openshape_sketch_solver::resolve_curves(&out.sketch).unwrap();
    let sweep = |k: u32| {
        curves
            .iter()
            .find_map(|(id, _, c)| match c {
                openshape_sketch_solver::ResolvedCurve::Arc { sweep_angle, .. } if *id == k => {
                    Some(*sweep_angle)
                }
                _ => None,
            })
            .unwrap()
    };
    let (arc, image) = pairs[1];
    assert!((sweep(arc) - sweep(image)).abs() < 1e-9);
}

#[test]
fn fully_constrained_sources_fix_their_images() {
    let (mut s, _, pairs, _) = every_kind();
    let mut next = max_key(&s);
    let points: Vec<u32> = pairs
        .iter()
        .flat_map(|(src, _)| s.entity(*src).unwrap().geometry.referenced_points())
        .collect();
    for p in points {
        let at = s.point(p).unwrap();
        constrain(&mut s, &mut next, Constraint::Fix { point: p, at });
    }
    let (circle, _) = pairs[2];
    constrain(
        &mut s,
        &mut next,
        Constraint::Radius {
            curve: circle,
            value: 0.3,
        },
    );
    let point = pairs[3].0;
    let at = s.point(point).unwrap();
    constrain(&mut s, &mut next, Constraint::Fix { point, at });
    // Fixing all three points of an arc would be redundant with its own
    // equal-radius condition: its end is held by a vertical distance from
    // the centre instead.
    let (arc, arc_image) = pairs[1];
    let Geometry::Arc { center, end, .. } = s.entity(arc).unwrap().geometry else {
        unreachable!()
    };
    s.constraints
        .retain(|c| !matches!(c.constraint, Constraint::Fix { point, .. } if point == end));
    let dy = s.point(end).unwrap().y - s.point(center).unwrap().y;
    constrain(
        &mut s,
        &mut next,
        Constraint::Distance {
            a: center,
            b: end,
            kind: openshape_sketch_solver::DistanceKind::Vertical,
            value: dy.abs(),
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    assert_eq!(out.dof, 0);
    assert!(out
        .entities
        .values()
        .all(|&e| e == EntityState::Constrained));
    assert_eq!(out.entities[&arc_image], EntityState::Constrained);
    assert_satisfied(&out.sketch);
}

#[test]
fn a_dimension_on_an_image_of_a_dimensioned_curve_is_redundant() {
    let (mut s, _, pairs, mirror) = every_kind();
    let mut next = max_key(&s);
    let (arc, arc_image) = pairs[1];
    let (line, line_image) = pairs[0];
    constrain(
        &mut s,
        &mut next,
        Constraint::Radius {
            curve: arc,
            value: 0.4,
        },
    );
    let image_radius = constrain(
        &mut s,
        &mut next,
        Constraint::Radius {
            curve: arc_image,
            value: 0.4,
        },
    );
    let length = (DVec2::new(1.0, 1.6) - DVec2::new(0.2, 1.0)).length();
    constrain(
        &mut s,
        &mut next,
        Constraint::Length {
            line,
            value: length,
        },
    );
    let image_length = constrain(
        &mut s,
        &mut next,
        Constraint::Length {
            line: line_image,
            value: length,
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Redundant);
    assert_eq!(out.redundant, vec![image_radius, image_length]);
    assert!(!out.redundant.contains(&mirror));
    assert_satisfied(&out.sketch);
}

#[test]
fn a_dimension_on_an_image_drives_its_source() {
    let (mut s, _, pairs, _) = every_kind();
    let mut next = max_key(&s);
    let (line, line_image) = pairs[0];
    constrain(
        &mut s,
        &mut next,
        Constraint::Length {
            line: line_image,
            value: 0.25,
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    let len = |k: u32| match out.sketch.entity(k).unwrap().geometry {
        Geometry::Line { start, end } => {
            (out.sketch.point(end).unwrap() - out.sketch.point(start).unwrap()).length()
        }
        _ => unreachable!(),
    };
    assert!((len(line) - 0.25).abs() < 1e-9);
    assert!((len(line_image) - 0.25).abs() < 1e-9);
}

#[test]
fn dragging_an_image_moves_its_source_to_the_mirror_position() {
    let (s, _, pairs, _) = every_kind();
    let settled = solve(&s, &opts()).sketch;
    let (point, image) = pairs[3];
    let target = DVec2::new(1.0, -0.8);
    let out = solve(
        &settled,
        &SolveOptions {
            drag: vec![(image, target)],
            ..opts()
        },
    );
    assert_eq!(out.status, SolveStatus::Solved);
    let at = |k: u32| out.sketch.point(k).unwrap();
    // The drag pulls the image to the target, its source resisting with a
    // millionth of the weight (design note 0005 §7)...
    assert!((at(image) - target).length() < 1e-5);
    // ...and the source is exactly its mirror image.
    let expected = reflect(at(image), DVec2::ZERO, DVec2::new(2.0, 1.0));
    assert!((at(point) - expected).length() < 1e-9);
}

/// A horizontal axis fixed from (0, 0) to (4, 0) and a line from (1, 1) to
/// a point `p` that the image shares. `p` starts at `p_at`. Returns the
/// sketch, the axis, the line, its image, `p` and the axis's end.
fn shared(p_at: DVec2) -> (Sketch<u32>, u32, u32, u32, u32, u32) {
    let mut b = B::new();
    let (axis, a0, a1) = b.line(0.0, 0.0, 4.0, 0.0);
    b.fix(a0);
    b.fix(a1);
    let a = b.point(1.0, 1.0);
    let p = b.point(p_at.x, p_at.y);
    let line = b.line_between(a, p);
    let mut s = b.build();
    let mut next = max_key(&s);
    let pairs = add_images(
        &mut s,
        &[line],
        (DVec2::ZERO, DVec2::new(4.0, 0.0)),
        &[p],
        0.0,
        &mut Rng(1),
        &mut next,
    );
    (s, axis, line, pairs[0].1, p, a1)
}

#[test]
fn a_shared_point_is_held_on_the_axis() {
    let (mut s, axis, line, image, p, _) = shared(DVec2::new(2.0, 0.0));
    let mut next = max_key(&s);
    let mirror = constrain(
        &mut s,
        &mut next,
        Constraint::Mirror {
            axis,
            pairs: vec![(line, image)],
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    // a: 2, p: 2, less p's distance from the axis.
    assert_eq!(out.dof, 3);
    // Dragged off the axis, the shared point slides along it.
    let out = solve(
        &s,
        &SolveOptions {
            drag: vec![(p, DVec2::new(2.5, 0.7))],
            ..opts()
        },
    );
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(out.sketch.point(p).unwrap().y.abs() < 1e-9);
    assert_satisfied(&out.sketch);
    // The same point already held on the axis by a constraint: the mirror
    // implies it, so one of the two is redundant (the mirror tool removes
    // the point-on-line constraint it replaces).
    constrain(
        &mut s,
        &mut next,
        Constraint::PointOnCurve {
            point: p,
            curve: axis,
        },
    );
    s.constraints.rotate_right(1);
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Redundant);
    assert_eq!(out.redundant, vec![mirror]);
}

#[test]
fn a_shared_point_at_an_end_of_the_axis_adds_no_equation() {
    // The shared point is the axis's own end.
    let mut b = B::new();
    let (axis, a0, a1) = b.line(0.0, 0.0, 4.0, 0.0);
    b.fix(a0);
    b.fix(a1);
    let a = b.point(5.0, 1.0);
    let line = b.line_between(a, a1);
    let mut s = b.build();
    let mut next = max_key(&s);
    let pairs = add_images(
        &mut s,
        &[line],
        (DVec2::ZERO, DVec2::new(4.0, 0.0)),
        &[a1],
        0.0,
        &mut Rng(2),
        &mut next,
    );
    constrain(&mut s, &mut next, Constraint::Mirror { axis, pairs });
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert_eq!(out.dof, 2);

    // A point coincident with the axis's end.
    let (mut s, axis, line, image, p, a1) = shared(DVec2::new(4.0, 0.0));
    let mut next = max_key(&s);
    constrain(&mut s, &mut next, Constraint::Coincident { a: p, b: a1 });
    constrain(
        &mut s,
        &mut next,
        Constraint::Mirror {
            axis,
            pairs: vec![(line, image)],
        },
    );
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    assert_eq!(out.dof, 2);
    assert_satisfied(&out.sketch);
}

#[test]
fn coincident_shared_points_are_held_on_the_axis_once() {
    // Two lines ending at coincident points on the axis (a joint drawn as
    // two line ends, as a slot's centre line meets its arc's centre), each
    // shared with its image: holding both on the axis would say the same
    // twice. Found by the mirror tool's fuzz.
    let mut b = B::new();
    let (axis, a0, a1) = b.line(0.0, 0.0, 4.0, 0.0);
    b.fix(a0);
    b.fix(a1);
    let a = b.point(1.0, 1.0);
    let p = b.point(2.0, 0.0);
    let first = b.line_between(a, p);
    let c = b.point(3.0, 2.0);
    let q = b.point(2.0, 0.0);
    let second = b.line_between(c, q);
    b.coincident(p, q);
    let mut s = b.build();
    let mut next = max_key(&s);
    let pairs = add_images(
        &mut s,
        &[first, second],
        (DVec2::ZERO, DVec2::new(4.0, 0.0)),
        &[p, q],
        0.0,
        &mut Rng(3),
        &mut next,
    );
    constrain(&mut s, &mut next, Constraint::Mirror { axis, pairs });
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.redundant);
    // a, c: 4; the joint: 2, less its distance from the axis.
    assert_eq!(out.dof, 5);
    assert_satisfied(&out.sketch);
}

#[test]
fn deleting_a_pair_frees_only_that_image() {
    let (s, _, pairs, mirror) = every_kind();
    let c = s.constraint(mirror).unwrap().constraint.clone();
    let (line, _) = pairs[0];
    let shrunk = c.without(|k| k == line).unwrap();
    let Constraint::Mirror { pairs: kept, .. } = &shrunk else {
        unreachable!()
    };
    assert_eq!(kept.as_slice(), &pairs[1..]);
    let mut s2 = s.clone();
    s2.constraints.last_mut().unwrap().constraint = shrunk;
    let out = solve(&s2, &opts());
    assert_eq!(out.status, SolveStatus::Solved);
    // The line's image is a free line again.
    assert_eq!(out.dof, 28 + 4);
}

#[test]
fn a_curve_leaving_an_image_is_mirrored_back_onto_its_source() {
    // L from a to b, mirrored to L' from a' to b'; then M from b' to c is
    // mirrored to M' from b to c'. The pair (b', b) is (b, b') read
    // backwards: the same relation, counted once.
    let mut b = B::new();
    let (axis, a0, a1) = b.line(0.0, 0.0, 0.0, 4.0);
    b.fix(a0);
    b.fix(a1);
    let (l, _, pb) = b.line(1.0, 1.0, 2.0, 3.0);
    let mut s = b.build();
    let mut next = max_key(&s);
    let ends = (DVec2::ZERO, DVec2::new(0.0, 4.0));
    let first = add_images(&mut s, &[l], ends, &[], 0.0, &mut Rng(4), &mut next);
    let Geometry::Line { end: pb_image, .. } = s.entity(first[0].1).unwrap().geometry else {
        unreachable!()
    };
    let mut point = |s: &mut Sketch<u32>, at: DVec2| {
        next += 1;
        s.entities.push(Entity {
            id: next,
            geometry: Geometry::Point { at },
            construction: false,
        });
        next
    };
    let c = point(&mut s, DVec2::new(-3.0, 2.0));
    let c_image = point(&mut s, DVec2::new(3.0, 2.0));
    let mut next = max_key(&s);
    next += 1;
    let m = next;
    s.entities.push(Entity {
        id: m,
        geometry: Geometry::Line {
            start: pb_image,
            end: c,
        },
        construction: false,
    });
    next += 1;
    let m_image = next;
    s.entities.push(Entity {
        id: m_image,
        geometry: Geometry::Line {
            start: pb,
            end: c_image,
        },
        construction: false,
    });
    let mut pairs = first;
    pairs.push((m, m_image));
    constrain(&mut s, &mut next, Constraint::Mirror { axis, pairs });
    let out = solve(&s, &opts());
    assert_eq!(out.status, SolveStatus::Solved, "{:?}", out.status);
    // a, b and c are free; everything else follows.
    assert_eq!(out.dof, 6);
    assert_satisfied(&out.sketch);
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
fn malformed_mirrors_are_refused() {
    let (s, axis, pairs, mirror) = every_kind();
    let with = |pairs: Vec<(u32, u32)>, axis: u32| {
        let mut s = s.clone();
        s.constraints.last_mut().unwrap().constraint = Constraint::Mirror { axis, pairs };
        s
    };
    let (line, line_image) = pairs[0];
    let (arc, arc_image) = pairs[1];
    let (circle, circle_image) = pairs[2];
    let (control, control_image) = pairs[4];
    let (through, through_image) = pairs[5];
    refused(&with(vec![], axis), "no entities");
    refused(&with(vec![(line, line)], axis), "its own image");
    refused(&with(vec![(line, arc_image)], axis), "not the same kind");
    refused(&with(vec![(circle, arc_image)], axis), "not the same kind");
    refused(&with(vec![(arc, circle_image)], axis), "not the same kind");
    refused(&with(vec![(control, through_image)], axis), "same form");
    refused(&with(vec![(through, control_image)], axis), "same form");
    refused(
        &with(vec![(line, line_image), (arc, line_image)], axis),
        "more than one",
    );
    refused(&with(vec![(axis, line_image)], axis), "mirror line");
    refused(&with(vec![(line, axis)], axis), "mirror line");
    refused(&with(vec![(line, line_image)], arc), "must be a line");
    refused(
        &with(vec![(line, line_image), (line_image, arc)], axis),
        "not the same kind",
    );
    // A chain inside one mirror: the image's points are sources too.
    let mut s3 = s.clone();
    let mut next = max_key(&s3);
    let third = add_images(
        &mut s3,
        &[line_image],
        (DVec2::ZERO, DVec2::new(2.0, 1.0)),
        &[],
        0.0,
        &mut Rng(3),
        &mut next,
    )[0]
    .1;
    s3.constraints
        .iter_mut()
        .find(|c| c.id == mirror)
        .unwrap()
        .constraint = Constraint::Mirror {
        axis,
        pairs: vec![(line, line_image), (line_image, third)],
    };
    refused(&s3, "both an image and a source");
    // An image point with two sources: the second image line shares the
    // first's start but mirrors a different line.
    let mut s4 = s.clone();
    let mut next = max_key(&s4);
    let src = {
        next += 1;
        let p = next;
        s4.entities.push(Entity {
            id: p,
            geometry: Geometry::Point {
                at: DVec2::new(3.0, 3.0),
            },
            construction: false,
        });
        next += 1;
        let q = next;
        s4.entities.push(Entity {
            id: q,
            geometry: Geometry::Point {
                at: DVec2::new(3.5, 3.0),
            },
            construction: false,
        });
        next += 1;
        s4.entities.push(Entity {
            id: next,
            geometry: Geometry::Line { start: p, end: q },
            construction: false,
        });
        next
    };
    let Geometry::Line {
        start: shared_start,
        ..
    } = s4.entity(line_image).unwrap().geometry
    else {
        unreachable!()
    };
    next += 1;
    let lone = next;
    s4.entities.push(Entity {
        id: lone,
        geometry: Geometry::Point {
            at: DVec2::new(0.0, -3.0),
        },
        construction: false,
    });
    next += 1;
    s4.entities.push(Entity {
        id: next,
        geometry: Geometry::Line {
            start: shared_start,
            end: lone,
        },
        construction: false,
    });
    let second_image = next;
    s4.constraints
        .iter_mut()
        .find(|c| c.id == mirror)
        .unwrap()
        .constraint = Constraint::Mirror {
        axis,
        pairs: vec![(line, line_image), (src, second_image)],
    };
    refused(&s4, "image of both");
    // An end of the axis as an image.
    let mut s5 = s.clone();
    let Geometry::Line { end: axis_end, .. } = s5.entity(axis).unwrap().geometry else {
        unreachable!()
    };
    let lone_point = pairs[3].0;
    s5.constraints
        .iter_mut()
        .find(|c| c.id == mirror)
        .unwrap()
        .constraint = Constraint::Mirror {
        axis,
        pairs: vec![(lone_point, axis_end)],
    };
    refused(&s5, "end of the mirror line");
}

#[test]
fn images_across_mirrors_must_not_overlap_or_cycle() {
    let (s, axis, pairs, _) = every_kind();
    let (line, line_image) = pairs[0];
    let mut next = max_key(&s);
    // The same image in two mirrors.
    let mut s1 = s.clone();
    constrain(
        &mut s1,
        &mut next,
        Constraint::Mirror {
            axis,
            pairs: vec![(line, line_image)],
        },
    );
    refused(&s1, "already an image");
    // A cycle: the image mirrored back onto its source by a second mirror.
    let mut s2 = s.clone();
    constrain(
        &mut s2,
        &mut next,
        Constraint::Mirror {
            axis,
            pairs: vec![(line_image, line)],
        },
    );
    refused(&s2, "derived from its own image");
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Random constrained sketches with random curves mirrored about a
    /// random (fixed or free) axis: the mirror adds no freedom and no
    /// redundancy, and everything holds after the solve.
    #[test]
    fn random_mirrors_add_no_freedom_and_no_redundancy(
        seed in any::<u64>(),
        motifs in 1usize..7,
        e in -3i32..=2,
        fixed_axis in any::<bool>(),
        jitter in prop_oneof![Just(0.0), 0.0f64..0.2],
        amount in prop_oneof![Just(0.0), 0.0f64..0.02],
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::ZERO);
        let mut s = perturb(&g.sketch, seed, amount, scale);
        let mut rng = Rng(seed ^ 0xA5A5);
        let mut next = max_key(&s);
        // The axis: a new line, long enough not to be degenerate.
        let a = rng.pt() * (2.0 * scale);
        let b = a + rng.dir() * (scale * rng.range(0.5, 3.0));
        let mut point = |s: &mut Sketch<u32>, at: DVec2| {
            next += 1;
            s.entities.push(Entity { id: next, geometry: Geometry::Point { at }, construction: false });
            next
        };
        let (a0, a1) = (point(&mut s, a), point(&mut s, b));
        next += 1;
        let axis = next;
        s.entities.push(Entity {
            id: axis,
            geometry: Geometry::Line { start: a0, end: a1 },
            construction: true,
        });
        if fixed_axis {
            constrain(&mut s, &mut next, Constraint::Fix { point: a0, at: a });
            constrain(&mut s, &mut next, Constraint::Fix { point: a1, at: b });
        }
        let before = solve(&s, &opts());
        prop_assume!(matches!(before.status, SolveStatus::Solved | SolveStatus::Redundant));
        // Sources: about half the curves, and lone points.
        let used: Vec<u32> = s.entities.iter().flat_map(|e| e.geometry.referenced_points()).collect();
        let sources: Vec<u32> = s
            .entities
            .iter()
            .filter(|e| e.id != axis)
            .filter(|e| e.geometry.is_curve() || !used.contains(&e.id))
            .map(|e| e.id)
            .filter(|_| rng.chance(0.5))
            .collect();
        prop_assume!(!sources.is_empty());
        // Start from the solved sketch, as the mirror tool does.
        let mut s = before.sketch.clone();
        let ends = (s.point(a0).unwrap(), s.point(a1).unwrap());
        let pairs = add_images(&mut s, &sources, ends, &[], jitter * scale, &mut rng, &mut next);
        constrain(&mut s, &mut next, Constraint::Mirror { axis, pairs });
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
