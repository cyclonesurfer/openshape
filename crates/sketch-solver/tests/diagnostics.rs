//! Redundant and conflicting constraints (design note 0005 §5, §8): an
//! appended constraint implied by others is reported and the others are
//! not; an appended inconsistent one makes the solve fail as `Conflicting`
//! with it in the conflicting set, and the geometry is returned unchanged.

mod common;

use common::{assert_satisfied, bits, B};
use openshape_sketch_solver::{
    measure, solve, Constraint, DVec2, DistanceKind, EntityState, SolveError, SolveOptions,
    SolveStatus,
};
use std::f64::consts::{FRAC_PI_2, PI};

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// A fully constrained 3 × 2 rectangle: lines with their own end points,
/// coincident corners, H/V sides, width, height and a fixed corner.
/// Returns (builder, lines [bottom, right, top, left], corners in order
/// [bottom-left, bottom-right, top-right, top-left], line end points).
fn rectangle() -> (B, [u32; 4], [u32; 4], [u32; 8]) {
    let mut b = B::new();
    let (l0, p0, p1) = b.line(0.0, 0.0, 3.0, 0.0);
    let (l1, p2, p3) = b.line(3.0, 0.0, 3.0, 2.0);
    let (l2, p4, p5) = b.line(3.0, 2.0, 0.0, 2.0);
    let (l3, p6, p7) = b.line(0.0, 2.0, 0.0, 0.0);
    b.coincident(p1, p2);
    b.coincident(p3, p4);
    b.coincident(p5, p6);
    b.coincident(p7, p0);
    b.c(Constraint::Horizontal { line: l0 });
    b.c(Constraint::Vertical { line: l1 });
    b.c(Constraint::Horizontal { line: l2 });
    b.c(Constraint::Vertical { line: l3 });
    b.c(Constraint::Length {
        line: l0,
        value: 3.0,
    });
    b.c(Constraint::Length {
        line: l1,
        value: 2.0,
    });
    b.fix(p0);
    (
        b,
        [l0, l1, l2, l3],
        [p0, p1, p4, p6],
        [p0, p1, p2, p3, p4, p5, p6, p7],
    )
}

#[test]
fn the_base_rectangle_is_well_constrained() {
    let (b, ..) = rectangle();
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 0);
    assert!(s.redundant.is_empty());
    assert!(s.entities.values().all(|&e| e == EntityState::Constrained));
}

#[test]
fn appended_implied_constraints_are_reported_alone() {
    let (base, lines, corners, _) = rectangle();
    let implied: Vec<Constraint<u32>> = vec![
        // Opposite sides of a rectangle are parallel and equal.
        Constraint::Parallel {
            a: lines[0],
            b: lines[2],
        },
        Constraint::Equal {
            a: lines[1],
            b: lines[3],
        },
        Constraint::Perpendicular {
            a: lines[0],
            b: lines[1],
        },
        // The diagonal is already determined.
        Constraint::Distance {
            a: corners[0],
            b: corners[2],
            kind: DistanceKind::Aligned,
            value: 13f64.sqrt(),
        },
        // So is the top edge's length and the top-left corner's height.
        Constraint::Length {
            line: lines[2],
            value: 3.0,
        },
        Constraint::Distance {
            a: corners[0],
            b: corners[3],
            kind: DistanceKind::Vertical,
            value: 2.0,
        },
        Constraint::Angle {
            a: lines[0],
            b: lines[1],
            value: FRAC_PI_2,
        },
        Constraint::HorizontalPoints {
            a: corners[2],
            b: corners[3],
        },
        Constraint::Fix {
            point: corners[2],
            at: DVec2::new(3.0, 2.0),
        },
    ];
    for c in implied {
        let mut b = base.clone();
        let id = b.c(c.clone());
        let s = solve(&b.build(), &opts());
        assert_eq!(s.status, SolveStatus::Redundant, "{c:?}");
        assert_eq!(s.redundant, vec![id], "only the appended {c:?}");
        assert!(s.conflicting.is_empty());
        assert_eq!(s.dof, 0);
        assert_satisfied(&s.sketch);
    }
}

#[test]
fn redundancy_follows_constraint_order() {
    // A horizontal line, then the same condition on its end points: the
    // later one is reported. In the other order, the other one is.
    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 1.0, 0.1);
    b.c(Constraint::Horizontal { line: l });
    let second = b.c(Constraint::HorizontalPoints { a: p, b: q });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant, vec![second]);

    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 1.0, 0.1);
    b.c(Constraint::HorizontalPoints { a: p, b: q });
    let second = b.c(Constraint::Horizontal { line: l });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.redundant, vec![second]);

    // A triangle's three sides fix its shape; a fourth dimension (the
    // height over the base) is reported, not the earlier sides.
    let mut b = B::new();
    let (base, a0, a1) = b.line(0.0, 0.0, 4.0, 0.0);
    let (_, b0, b1) = b.line(4.0, 0.0, 1.0, 3.0);
    let (_, c0, c1) = b.line(1.0, 3.0, 0.0, 0.0);
    b.coincident(a1, b0);
    b.coincident(b1, c0);
    b.coincident(c1, a0);
    for (p, q) in [(a0, a1), (b0, b1), (c0, c1)] {
        let value = (common::pos(&b.sketch, q) - common::pos(&b.sketch, p)).length();
        b.c(Constraint::Distance {
            a: p,
            b: q,
            kind: DistanceKind::Aligned,
            value,
        });
    }
    let height = b.c(Constraint::Distance {
        a: c0,
        b: base,
        kind: DistanceKind::Aligned,
        value: 3.0,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant, vec![height]);
}

#[test]
fn appended_inconsistent_constraints_conflict_and_leave_geometry() {
    let (base, lines, corners, ends) = rectangle();
    let inconsistent: Vec<Constraint<u32>> = vec![
        // The diagonal cannot be 4 in a 3 × 2 rectangle.
        Constraint::Distance {
            a: corners[0],
            b: corners[2],
            kind: DistanceKind::Aligned,
            value: 4.0,
        },
        // Nor the top 2.5 long.
        Constraint::Length {
            line: lines[2],
            value: 2.5,
        },
        // Nor can adjacent sides be parallel.
        Constraint::Parallel {
            a: lines[0],
            b: lines[1],
        },
        // Nor the angle between them 60°.
        Constraint::Angle {
            a: lines[0],
            b: lines[1],
            value: PI / 3.0,
        },
        // Nor −90° (the wrong way round).
        Constraint::Angle {
            a: lines[0],
            b: lines[1],
            value: -FRAC_PI_2,
        },
        // The fixed corner fixed elsewhere.
        Constraint::Fix {
            point: corners[0],
            at: DVec2::new(0.5, 0.0),
        },
        // An end point of the bottom on the top.
        Constraint::PointOnCurve {
            point: ends[0],
            curve: lines[2],
        },
    ];
    for c in inconsistent {
        let mut b = base.clone();
        // Start away from the solution so that "unchanged" is meaningful.
        for e in b.sketch.entities.iter_mut() {
            if let openshape_sketch_solver::Geometry::Point { at } = &mut e.geometry {
                *at += DVec2::new(0.01, -0.02);
            }
        }
        let id = b.c(c.clone());
        let input = b.build();
        let s = solve(&input, &opts());
        assert_eq!(
            s.status,
            SolveStatus::Failed(SolveError::Conflicting),
            "{c:?}"
        );
        assert!(s.conflicting.contains(&id), "{c:?}: {:?}", s.conflicting);
        assert_eq!(s.sketch, input, "geometry unchanged");
        assert_eq!(bits(&s.sketch), bits(&input));
        assert!(s.max_residual.is_finite() && s.max_residual > 0.0);
    }
}

#[test]
fn conflicts_among_eliminated_constraints() {
    // Two fixes of one point at different places: both are in the set.
    let mut b = B::new();
    let p = b.point(0.0, 0.0);
    b.c(Constraint::Fix {
        point: p,
        at: DVec2::new(0.0, 0.0),
    });
    let f2 = b.c(Constraint::Fix {
        point: p,
        at: DVec2::new(1.0, 0.0),
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Failed(SolveError::Conflicting));
    // The later fix conflicts in x and is implied in y: it is reported as
    // conflicting only.
    assert_eq!(s.conflicting, vec![f2]);
    assert!(s.redundant.is_empty());
    assert_eq!(s.sketch, b.build());

    // A horizontal line between two points fixed at different heights.
    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 1.0, 1.0);
    b.fix(p);
    b.fix(q);
    let h = b.c(Constraint::Horizontal { line: l });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Failed(SolveError::Conflicting));
    assert_eq!(s.conflicting, vec![h]);

    // Concentric circles with different fixed radii but equal: conflict.
    let mut b = B::new();
    let (c1, _) = b.circle(0.0, 0.0, 1.0);
    let (c2, _) = b.circle(0.0, 0.0, 2.0);
    b.c(Constraint::Radius {
        curve: c1,
        value: 1.0,
    });
    b.c(Constraint::Radius {
        curve: c2,
        value: 2.0,
    });
    let eq = b.c(Constraint::Equal { a: c1, b: c2 });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Failed(SolveError::Conflicting));
    assert_eq!(s.conflicting, vec![eq]);
}

#[test]
fn a_consistent_duplicate_dimension_is_redundant() {
    // The same length twice, measured: redundant, not conflicting, even
    // though the values were computed independently.
    let mut b = B::new();
    let (l, p, _) = b.line(0.1, 0.2, 1.3, 0.9);
    b.fix(p);
    let c = Constraint::Length {
        line: l,
        value: 1.0,
    };
    let v = measure(&b.sketch, &c).unwrap();
    b.c(Constraint::Length { line: l, value: v });
    let dup = b.c(Constraint::Length { line: l, value: v });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant, vec![dup]);
}

#[test]
fn duplicates_within_the_model_tolerance_are_redundant() {
    // Two fixes of a point 5e-8 m apart (under the 1e-7 m tolerance): the
    // input already satisfies both, so it is the solution, and the second
    // is reported as implied.
    let mut b = B::new();
    let p = b.point(0.0, 0.0);
    b.c(Constraint::Fix {
        point: p,
        at: DVec2::new(0.0, 0.0),
    });
    let second = b.c(Constraint::Fix {
        point: p,
        at: DVec2::new(5e-8, 0.0),
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant, vec![second]);
    assert_eq!(s.sketch, b.build());
    assert!(s.max_residual <= 1e-7);

    // Two lengths of one line differing by 5e-8.
    let mut b = B::new();
    let (l, p, _) = b.line(0.0, 0.0, 1.0, 0.0);
    b.fix(p);
    b.c(Constraint::Length {
        line: l,
        value: 1.0,
    });
    let second = b.c(Constraint::Length {
        line: l,
        value: 1.0 + 5e-8,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant, vec![second]);
    assert_eq!(s.sketch, b.build());

    // A horizontal line between points fixed 5e-8 apart in height.
    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 1.0, 5e-8);
    b.fix(p);
    b.fix(q);
    let h = b.c(Constraint::Horizontal { line: l });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Redundant);
    assert_eq!(s.redundant, vec![h]);

    // Twice the tolerance is a conflict.
    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 1.0, 2e-7);
    b.fix(p);
    b.fix(q);
    let h = b.c(Constraint::Horizontal { line: l });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Failed(SolveError::Conflicting));
    assert_eq!(s.conflicting, vec![h]);
}
