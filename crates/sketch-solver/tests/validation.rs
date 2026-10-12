//! Input validation (design note 0005 §4, §8): every malformed sketch is
//! rejected as `Invalid` (or `Degenerate`) with a message naming the
//! problem, the input returned unchanged and nothing diagnosed; random
//! garbage never panics and never produces a NaN.

mod common;

use common::gen::Rng;
use common::{bits, violation, B};
use openshape_sketch_solver::{
    measure, resolve_curves, solve, Constraint, ConstraintEntry, DVec2, DistanceKind, Entity,
    EntityState, Geometry, ResolvedCurve, Sketch, Solution, SolveError, SolveOptions, SolveStatus,
    TangentKind, RHO_MIN,
};
use proptest::prelude::*;
use std::f64::consts::{PI, TAU};

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

fn unknowns(s: &Sketch<u32>) -> usize {
    s.entities
        .iter()
        .map(|e| match e.geometry {
            Geometry::Point { .. } => 2,
            Geometry::Circle { .. } | Geometry::Conic { .. } => 1,
            _ => 0,
        })
        .sum()
}

/// Debug text and bits: equal for identical sketches, NaNs included.
fn same(a: &Sketch<u32>, b: &Sketch<u32>) -> bool {
    format!("{a:?}") == format!("{b:?}") && bits(a) == bits(b)
}

/// Checks an early failure: the input back, nothing solved or diagnosed.
fn check_early(input: &Sketch<u32>, s: &Solution<u32>) {
    assert!(same(&s.sketch, input), "geometry unchanged");
    assert_eq!(s.dof, unknowns(input));
    assert!(s.entities.values().all(|&e| e == EntityState::Free));
    assert!(s.redundant.is_empty() && s.conflicting.is_empty());
    assert_eq!(s.max_residual, 0.0);
    assert_eq!(s.iterations, 0);
}

fn expect_invalid(s: &Sketch<u32>, what: &str) {
    let out = solve(s, &opts());
    match &out.status {
        SolveStatus::Failed(SolveError::Invalid(msg)) => {
            assert!(msg.contains(what), "expected {what:?} in {msg:?}")
        }
        other => panic!("expected Invalid ({what}), got {other:?}"),
    }
    check_early(s, &out);
}

fn expect_degenerate(s: &Sketch<u32>, what: &str) {
    let out = solve(s, &opts());
    match &out.status {
        SolveStatus::Failed(SolveError::Degenerate(msg)) => {
            assert!(msg.contains(what), "expected {what:?} in {msg:?}")
        }
        other => panic!("expected Degenerate ({what}), got {other:?}"),
    }
    check_early(s, &out);
}

/// Handles to the entities of [`base`].
#[allow(dead_code)]
struct Base {
    b: B,
    line: u32,
    l0: u32,
    l1: u32,
    line2: u32,
    circle: u32,
    cc: u32,
    arc: u32,
    ac: u32,
    a0: u32,
    a1: u32,
    spline: u32,
    ctrl: Vec<u32>,
    p: u32,
}

/// A valid sketch with one of each kind of entity.
fn base() -> Base {
    let mut b = B::new();
    let (line, l0, l1) = b.line(0.0, 0.0, 1.0, 0.0);
    let (line2, _, _) = b.line(0.0, 1.0, 1.0, 1.5);
    let (circle, cc) = b.circle(3.0, 0.0, 1.0);
    let (arc, ac, a0, a1) = b.arc(DVec2::new(0.0, 3.0), 1.0, 0.3, 2.0);
    let ctrl: Vec<u32> = [(4.0, 4.0), (5.0, 5.0), (6.0, 4.0), (7.0, 5.0)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    let spline = b.entity(Geometry::ControlSpline {
        control: ctrl.clone(),
        degree: 3,
    });
    let p = b.point(5.0, -5.0);
    Base {
        b,
        line,
        l0,
        l1,
        line2,
        circle,
        cc,
        arc,
        ac,
        a0,
        a1,
        spline,
        ctrl,
        p,
    }
}

#[test]
fn the_base_sketch_is_valid() {
    let s = solve(&base().b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
}

/// The base sketch with entity `id`'s geometry replaced.
fn with_entity(f: impl FnOnce(&Base) -> (u32, Geometry<u32>)) -> Sketch<u32> {
    let bs = base();
    let (id, g) = f(&bs);
    let mut s = bs.b.build();
    s.entities.iter_mut().find(|e| e.id == id).unwrap().geometry = g;
    s
}

#[test]
fn bad_entities_are_invalid() {
    // Duplicate key.
    let bs = base();
    let mut s = bs.b.build();
    s.entities.push(Entity {
        id: bs.p,
        geometry: Geometry::Point { at: DVec2::ZERO },
        construction: false,
    });
    expect_invalid(&s, "duplicate entity key");

    // Point coordinates.
    for at in [
        DVec2::new(f64::NAN, 0.0),
        DVec2::new(0.0, f64::INFINITY),
        DVec2::new(f64::NEG_INFINITY, 0.0),
        DVec2::new(1.000_001e6, 0.0),
        DVec2::new(0.0, -2e6),
    ] {
        let s = with_entity(|b| (b.p, Geometry::Point { at }));
        expect_invalid(&s, "point coordinates must be finite");
    }
    // The limit itself is allowed.
    let s = with_entity(|b| {
        (
            b.p,
            Geometry::Point {
                at: DVec2::new(1e6, -1e6),
            },
        )
    });
    assert_eq!(solve(&s, &opts()).status, SolveStatus::Solved);

    // Lines.
    let s = with_entity(|b| {
        (
            b.line,
            Geometry::Line {
                start: b.l0,
                end: 999,
            },
        )
    });
    expect_invalid(&s, "end 999 does not exist");
    let s = with_entity(|b| {
        (
            b.line,
            Geometry::Line {
                start: b.line2,
                end: b.l1,
            },
        )
    });
    expect_invalid(&s, "is not a point");
    let s = with_entity(|b| {
        (
            b.line,
            Geometry::Line {
                start: b.l0,
                end: b.l0,
            },
        )
    });
    expect_invalid(&s, "same point");
    let s = with_entity(|b| {
        (
            b.line,
            Geometry::Line {
                start: b.line,
                end: b.l1,
            },
        )
    });
    expect_invalid(&s, "is not a point");

    // Circles.
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY, 2e6, -0.0] {
        let s = with_entity(|b| {
            (
                b.circle,
                Geometry::Circle {
                    center: b.cc,
                    radius,
                },
            )
        });
        expect_invalid(&s, "circle radius must be positive");
    }
    let s = with_entity(|b| {
        (
            b.circle,
            Geometry::Circle {
                center: b.line,
                radius: 1.0,
            },
        )
    });
    expect_invalid(&s, "centre");

    // Arcs.
    let s = with_entity(|b| {
        (
            b.arc,
            Geometry::Arc {
                center: b.ac,
                start: b.ac,
                end: b.a1,
            },
        )
    });
    expect_invalid(&s, "distinct");
    let s = with_entity(|b| {
        (
            b.arc,
            Geometry::Arc {
                center: b.ac,
                start: b.a0,
                end: b.a0,
            },
        )
    });
    expect_invalid(&s, "distinct");

    // Control splines.
    for degree in [0, 6, u32::MAX] {
        let s = with_entity(|b| {
            (
                b.spline,
                Geometry::ControlSpline {
                    control: b.ctrl.clone(),
                    degree,
                },
            )
        });
        expect_invalid(&s, "outside 1..=5");
    }
    let s = with_entity(|b| {
        (
            b.spline,
            Geometry::ControlSpline {
                control: b.ctrl[..3].to_vec(),
                degree: 3,
            },
        )
    });
    expect_invalid(&s, "fewer than degree + 1");
    let s = with_entity(|b| {
        (
            b.spline,
            Geometry::ControlSpline {
                control: vec![b.ctrl[0], b.circle, b.ctrl[2], b.ctrl[3]],
                degree: 3,
            },
        )
    });
    expect_invalid(&s, "control point");

    // Interpolated splines.
    let interp = |through: Vec<u32>, params: Vec<f64>| {
        with_entity(move |b| {
            let through = through.into_iter().map(|i| b.ctrl[i as usize]).collect();
            (b.spline, Geometry::InterpolatedSpline { through, params })
        })
    };
    expect_invalid(&interp(vec![0], vec![0.0]), "at least 2 points");
    expect_invalid(&interp(vec![0, 1, 2], vec![0.0, 1.0]), "2 parameters for 3");
    expect_invalid(&interp(vec![0, 1, 2], vec![0.1, 0.5, 1.0]), "from 0 to 1");
    expect_invalid(&interp(vec![0, 1, 2], vec![0.0, 0.5, 0.9]), "from 0 to 1");
    expect_invalid(
        &interp(vec![0, 1, 2, 3], vec![0.0, 0.5, 0.5, 1.0]),
        "strictly increasing",
    );
    expect_invalid(
        &interp(vec![0, 1, 2, 3], vec![0.0, 0.6, 0.4, 1.0]),
        "strictly increasing",
    );
    expect_invalid(&interp(vec![0, 1, 2], vec![0.0, f64::NAN, 1.0]), "finite");
    expect_invalid(
        &interp(vec![0, 1, 2], vec![0.0, f64::INFINITY, 1.0]),
        "finite",
    );
}

/// The base sketch plus one constraint.
fn with(f: impl FnOnce(&Base) -> Constraint<u32>) -> Sketch<u32> {
    let mut bs = base();
    let c = f(&bs);
    bs.b.c(c);
    bs.b.build()
}

#[test]
fn bad_references_are_invalid() {
    use Constraint as C;
    let cases: Vec<(Sketch<u32>, &str)> = vec![
        (
            with(|_| C::Coincident { a: 999, b: 1 }),
            "999 does not exist",
        ),
        (with(|b| C::Coincident { a: b.p, b: b.p }), "same entity"),
        (
            with(|b| C::Coincident { a: b.p, b: b.line }),
            "must be a point",
        ),
        (
            with(|b| C::PointOnCurve {
                point: b.p,
                curve: b.l0,
            }),
            "must be a curve",
        ),
        (
            with(|b| C::PointOnCurve {
                point: b.line,
                curve: b.circle,
            }),
            "must be a point",
        ),
        (
            with(|b| C::PointOnCurve {
                point: b.cc,
                curve: b.circle,
            }),
            "cannot lie on it",
        ),
        (with(|b| C::Horizontal { line: b.p }), "must be a line"),
        (with(|b| C::Vertical { line: b.circle }), "must be a line"),
        (
            with(|b| C::HorizontalPoints { a: b.l0, b: b.l0 }),
            "same entity",
        ),
        (
            with(|b| C::VerticalPoints { a: b.l0, b: b.arc }),
            "must be a point",
        ),
        (
            with(|b| C::Parallel {
                a: b.line,
                b: b.line,
            }),
            "same entity",
        ),
        (
            with(|b| C::Perpendicular {
                a: b.line,
                b: b.arc,
            }),
            "must be a line",
        ),
        (
            with(|b| C::Equal {
                a: b.line,
                b: b.circle,
            }),
            "equal needs",
        ),
        (
            with(|b| C::Equal {
                a: b.spline,
                b: b.line,
            }),
            "equal needs",
        ),
        (
            with(|b| C::Midpoint {
                point: b.l1,
                line: b.line,
            }),
            "cannot be its midpoint",
        ),
        (
            with(|b| C::Concentric {
                a: b.circle,
                b: b.line,
            }),
            "circle or an arc",
        ),
        (
            with(|b| C::Symmetric {
                a: b.p,
                b: b.cc,
                axis: b.circle,
            }),
            "must be a line",
        ),
        (
            with(|b| C::Symmetric {
                a: b.p,
                b: b.p,
                axis: b.line,
            }),
            "same entity",
        ),
        (
            with(|b| C::Tangent {
                a: b.line,
                b: b.arc,
                kind: TangentKind::AtPoint { point: b.p },
            }),
            "not an end point",
        ),
        (
            with(|b| C::Tangent {
                a: b.line,
                b: b.circle,
                kind: TangentKind::AtPoint { point: b.l0 },
            }),
            "has no end points",
        ),
        (
            with(|b| C::Tangent {
                a: b.line,
                b: b.line2,
                kind: TangentKind::LineCircle { center_left: true },
            }),
            "line-circle tangency",
        ),
        (
            with(|b| C::Tangent {
                a: b.circle,
                b: b.line,
                kind: TangentKind::CircleCircle { internal: false },
            }),
            "circle-circle tangency",
        ),
        (
            with(|b| C::Tangent {
                a: b.arc,
                b: b.arc,
                kind: TangentKind::CircleCircle { internal: true },
            }),
            "same entity",
        ),
        (
            with(|b| C::Tangent {
                a: b.p,
                b: b.line,
                kind: TangentKind::AtPoint { point: b.l0 },
            }),
            "must be a curve",
        ),
        (
            with(|b| C::Normal {
                line: b.circle,
                curve: b.circle,
            }),
            "must be a line",
        ),
        (
            with(|b| C::Normal {
                line: b.line,
                curve: b.spline,
            }),
            "circle or an arc",
        ),
        (
            with(|b| C::Fix {
                point: b.line,
                at: DVec2::ZERO,
            }),
            "must be a point",
        ),
        (
            with(|b| C::Distance {
                a: b.p,
                b: b.circle,
                kind: DistanceKind::Aligned,
                value: 1.0,
            }),
            "a distance needs",
        ),
        (
            with(|b| C::Distance {
                a: b.p,
                b: b.line,
                kind: DistanceKind::Horizontal,
                value: 1.0,
            }),
            "need two points",
        ),
        (
            with(|b| C::Length {
                line: b.arc,
                value: 1.0,
            }),
            "must be a line",
        ),
        (
            with(|b| C::Angle {
                a: b.line,
                b: b.p,
                value: 1.0,
            }),
            "must be a line",
        ),
        (
            with(|b| C::Radius {
                curve: b.line,
                value: 1.0,
            }),
            "circle or an arc",
        ),
        (
            with(|b| C::Diameter {
                curve: b.spline,
                value: 1.0,
            }),
            "circle or an arc",
        ),
    ];
    for (s, what) in cases {
        expect_invalid(&s, what);
    }

    // Tangency at a spline end shared by both ends (a closed spline).
    let mut bs = base();
    let pts: Vec<u32> = [(0.0, -2.0), (1.0, -3.0), (2.0, -2.0)]
        .iter()
        .map(|&(x, y)| bs.b.point(x, y))
        .collect();
    let closed = bs.b.entity(Geometry::ControlSpline {
        control: vec![pts[0], pts[1], pts[2], pts[0]],
        degree: 2,
    });
    let (l, ls, _) = bs.b.line(0.0, -2.0, -1.0, -2.0);
    bs.b.coincident(ls, pts[0]);
    bs.b.c(Constraint::Tangent {
        a: l,
        b: closed,
        kind: TangentKind::AtPoint { point: ls },
    });
    expect_invalid(&bs.b.build(), "ambiguous");

    // Duplicate constraint key.
    let mut s = with(|b| Constraint::Horizontal { line: b.line });
    let dup = s.constraints[0].clone();
    s.constraints.push(ConstraintEntry {
        id: dup.id,
        constraint: Constraint::Vertical { line: 3 },
    });
    expect_invalid(&s, "duplicate constraint key");
}

#[test]
fn bad_values_are_invalid() {
    use Constraint as C;
    // Each makes a dimension of one kind with the given value.
    type Dim = fn(&Base, f64) -> Constraint<u32>;
    let dims: [Dim; 6] = [
        |b, value| C::Distance {
            a: b.p,
            b: b.cc,
            kind: DistanceKind::Aligned,
            value,
        },
        |b, value| C::Distance {
            a: b.p,
            b: b.cc,
            kind: DistanceKind::Vertical,
            value,
        },
        |b, value| C::Distance {
            a: b.p,
            b: b.line,
            kind: DistanceKind::Aligned,
            value,
        },
        |b, value| C::Length {
            line: b.line,
            value,
        },
        |b, value| C::Radius {
            curve: b.circle,
            value,
        },
        |b, value| C::Diameter {
            curve: b.arc,
            value,
        },
    ];
    for d in &dims {
        // A zero dimension points to the coincident constraint.
        let s = with(|b| d(b, 0.0));
        expect_invalid(&s, "use a coincident constraint instead");
        expect_invalid(&with(|b| d(b, -0.0)), "coincident");
        expect_invalid(&with(|b| d(b, -1.0)), "must be positive");
        expect_invalid(&with(|b| d(b, f64::NAN)), "must be finite");
        expect_invalid(&with(|b| d(b, f64::INFINITY)), "must be finite");
        expect_invalid(&with(|b| d(b, 1.5e6)), "exceeds");
    }
    for value in [f64::NAN, f64::INFINITY, -PI, PI + 1e-12, 4.0, -TAU] {
        let s = with(|b| C::Angle {
            a: b.line,
            b: b.line2,
            value,
        });
        expect_invalid(&s, "(−π, π]");
    }
    // π itself is in range.
    let s = with(|b| C::Angle {
        a: b.line,
        b: b.line2,
        value: PI,
    });
    assert!(!matches!(
        solve(&s, &opts()).status,
        SolveStatus::Failed(SolveError::Invalid(_))
    ));
    for at in [
        DVec2::new(f64::NAN, 0.0),
        DVec2::new(0.0, f64::NEG_INFINITY),
        DVec2::new(3e6, 0.0),
    ] {
        let s = with(|b| C::Fix { point: b.p, at });
        expect_invalid(&s, "fixed position must be finite");
    }
}

#[test]
fn degenerate_geometry_is_reported() {
    // A zero-length line (two point entities at one place) cannot be
    // parallel to anything.
    let mut b = B::new();
    let (l, ..) = b.line(1.0, 1.0, 1.0, 1.0);
    let (m, ..) = b.line(0.0, 0.0, 1.0, 0.0);
    let c = b.c(Constraint::Parallel { a: l, b: m });
    let msg = format!("constraint {c:?} needs a direction that has zero length");
    expect_degenerate(&b.build(), &msg);

    let mut b = B::new();
    let (l, ..) = b.line(0.0, 0.0, 1.0, 0.0);
    let (m, ..) = b.line(2.0, 2.0, 2.0, 2.0);
    b.c(Constraint::Angle {
        a: l,
        b: m,
        value: 0.5,
    });
    expect_degenerate(&b.build(), "zero length");

    // A line tangent to an arc whose start sits on its centre.
    let mut b = B::new();
    let cc = b.point(0.0, 0.0);
    let s0 = b.point(0.0, 0.0);
    let e0 = b.point(1.0, 0.0);
    let arc = b.entity(Geometry::Arc {
        center: cc,
        start: s0,
        end: e0,
    });
    let (l, ls, _) = b.line(0.0, 0.0, -1.0, 0.0);
    b.coincident(ls, s0);
    b.c(Constraint::Tangent {
        a: l,
        b: arc,
        kind: TangentKind::AtPoint { point: ls },
    });
    let out = solve(&b.build(), &opts());
    assert!(
        matches!(out.status, SolveStatus::Failed(SolveError::Degenerate(_))),
        "{:?}",
        out.status
    );
    check_early(&b.build(), &out);
}

#[test]
fn resolve_and_measure_reject_bad_input() {
    let s = with_entity(|b| {
        (
            b.line,
            Geometry::Line {
                start: b.l0,
                end: 999,
            },
        )
    });
    assert!(matches!(resolve_curves(&s), Err(SolveError::Invalid(_))));
    let s = base().b.build();
    let bs = base();
    let bad: Vec<Constraint<u32>> = vec![
        Constraint::Horizontal { line: bs.line },
        Constraint::Length {
            line: 999,
            value: 1.0,
        },
        Constraint::Length {
            line: bs.p,
            value: 1.0,
        },
        Constraint::Radius {
            curve: bs.line,
            value: 1.0,
        },
        Constraint::Angle {
            a: bs.line,
            b: bs.line,
            value: 0.0,
        },
    ];
    for c in bad {
        assert!(
            matches!(measure(&s, &c), Err(SolveError::Invalid(_))),
            "{c:?}"
        );
    }
}

// ---------------------------------------------------------------------
// Garbage fuzzing.

/// Random sketch content: references and values are well formed except
/// with probability about `h` each (so `h = 0` gives valid sketches with
/// random, mostly unsatisfiable constraints, and larger `h` mostly
/// malformed ones).
struct Junk {
    r: Rng,
    h: f64,
}

impl Junk {
    /// A value in `lo..hi`, or (rarely) a hostile one.
    fn value(&mut self, lo: f64, hi: f64) -> f64 {
        if !self.r.chance(self.h) {
            return self.r.range(lo, hi);
        }
        const NASTY: [f64; 16] = [
            0.0,
            -0.0,
            -1.0,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            1e-300,
            1e-15,
            1e6,
            1.000_000_1e6,
            2e6,
            f64::MAX,
            f64::MIN_POSITIVE,
            PI,
            -PI,
            1e-7,
        ];
        NASTY[self.r.int(NASTY.len())]
    }

    /// A key from `pool`, or (rarely) any key at all, existing or not.
    fn pick(&mut self, pool: &[u32], n_keys: u32) -> u32 {
        if !pool.is_empty() && !self.r.chance(self.h) {
            pool[self.r.int(pool.len())]
        } else {
            self.r.int(n_keys as usize + 2) as u32
        }
    }

    /// Two keys from `pool`, distinct unless (rarely) not.
    fn pair(&mut self, pool: &[u32], n_keys: u32) -> (u32, u32) {
        let a = self.pick(pool, n_keys);
        let mut b = self.pick(pool, n_keys);
        if pool.iter().any(|&k| k != a) && !self.r.chance(self.h) {
            while b == a {
                b = pool[self.r.int(pool.len())];
            }
        }
        (a, b)
    }
}

fn garbage(seed: u64) -> Sketch<u32> {
    let mut r = Rng(seed);
    let h = [0.0, 0.0, 0.01, 0.05, 0.2][r.int(5)];
    let mut j = Junk { r, h };
    let mut s = Sketch::<u32>::default();
    let mut next = 1u32;
    let key = |j: &mut Junk, next: &mut u32| -> u32 {
        // Rarely, reuse a key.
        if *next > 2 && j.r.chance(j.h * 0.1) {
            1 + j.r.int(*next as usize - 1) as u32
        } else {
            *next += 1;
            *next - 1
        }
    };
    let mut points: Vec<u32> = Vec::new();
    let mut positions: Vec<DVec2> = Vec::new();
    for _ in 0..3 + j.r.int(8) {
        let at = match j.r.int(20) {
            0 if !positions.is_empty() => positions[j.r.int(positions.len())],
            1 if !positions.is_empty() => {
                positions[j.r.int(positions.len())] + DVec2::new(1e-13, -1e-13)
            }
            _ => DVec2::new(j.value(-5.0, 5.0), j.value(-5.0, 5.0)),
        };
        let id = key(&mut j, &mut next);
        s.entities.push(Entity {
            id,
            geometry: Geometry::Point { at },
            construction: false,
        });
        points.push(id);
        positions.push(at);
    }
    let n_pts = next;
    let mut lines = Vec::new();
    let mut circulars = Vec::new();
    let mut splines = Vec::new();
    let mut ellipses = Vec::new();
    let mut conics = Vec::new();
    let mut ends: Vec<(u32, Vec<u32>)> = Vec::new();
    for i in 0..3 + j.r.int(6) {
        // At least two lines and a circle or arc, so that most constraints
        // have something to refer to.
        let pick = match i {
            0 | 1 => 0,
            2 => 4 + 2 * j.r.int(2),
            _ => j.r.int(13),
        };
        let (geometry, kind, curve_ends) = match pick {
            0..=3 => {
                let (start, end) = j.pair(&points, n_pts);
                (Geometry::Line { start, end }, 0, vec![start, end])
            }
            4 | 5 => (
                Geometry::Circle {
                    center: j.pick(&points, n_pts),
                    radius: j.value(0.1, 3.0),
                },
                1,
                vec![],
            ),
            6 | 7 => {
                let center = j.pick(&points, n_pts);
                let rest: Vec<u32> = points.iter().copied().filter(|&p| p != center).collect();
                let (start, end) = j.pair(&rest, n_pts);
                (Geometry::Arc { center, start, end }, 1, vec![start, end])
            }
            8 => {
                let n = 2 + j.r.int(5);
                let degree = if !j.r.chance(j.h) {
                    1 + j.r.int(n.min(6) - 1) as u32
                } else {
                    [0, 6, 7, u32::MAX][j.r.int(4)]
                };
                let control: Vec<u32> = (0..n).map(|_| j.pick(&points, n_pts)).collect();
                let e = vec![control[0], control[n - 1]];
                (Geometry::ControlSpline { control, degree }, 2, e)
            }
            10 => {
                let (center, major) = j.pair(&points, n_pts);
                let minor = j.pick(&points, n_pts);
                (
                    Geometry::Ellipse {
                        center,
                        major,
                        minor,
                    },
                    3,
                    vec![],
                )
            }
            11 => {
                let (center, major) = j.pair(&points, n_pts);
                let minor = j.pick(&points, n_pts);
                let (start, end) = j.pair(&points, n_pts);
                (
                    Geometry::EllipticalArc {
                        center,
                        major,
                        minor,
                        start,
                        end,
                    },
                    3,
                    vec![start, end],
                )
            }
            12 => {
                let (start, end) = j.pair(&points, n_pts);
                let control = j.pick(&points, n_pts);
                (
                    Geometry::Conic {
                        start,
                        control,
                        end,
                        rho: j.value(0.0, 1.0),
                    },
                    4,
                    vec![start, end],
                )
            }
            _ => {
                let n = 2 + j.r.int(5);
                let mut params: Vec<f64> = (1..n - 1).map(|_| j.r.f()).collect();
                params.push(0.0);
                params.push(1.0);
                params.sort_by(f64::total_cmp);
                if j.r.chance(j.h) {
                    let k = j.r.int(n);
                    params[k] = j.value(0.0, 1.0);
                }
                let m = if j.r.chance(j.h) { n + 1 } else { n };
                let through: Vec<u32> = (0..m).map(|_| j.pick(&points, n_pts)).collect();
                let e = vec![through[0], through[m - 1]];
                (Geometry::InterpolatedSpline { through, params }, 2, e)
            }
        };
        let id = key(&mut j, &mut next);
        ends.push((id, curve_ends));
        match kind {
            0 => lines.push(id),
            1 => circulars.push(id),
            3 => ellipses.push(id),
            4 => conics.push(id),
            _ => splines.push(id),
        }
        let construction = j.r.chance(0.2);
        s.entities.push(Entity {
            id,
            geometry,
            construction,
        });
    }
    let n_keys = next;
    let curves: Vec<u32> = lines
        .iter()
        .chain(&circulars)
        .chain(&splines)
        .chain(&ellipses)
        .chain(&conics)
        .copied()
        .collect();
    // Circles and arcs, and sometimes ellipses where a circular is wanted.
    let roundish: Vec<u32> = circulars.iter().chain(&ellipses).copied().collect();
    let with_ends: Vec<u32> = ends
        .iter()
        .filter(|(_, e)| !e.is_empty())
        .map(|(id, _)| *id)
        .collect();
    for _ in 0..j.r.int(9) {
        use Constraint as C;
        let mut extra: Option<Constraint<u32>> = None;
        let constraint = match j.r.int(21) {
            0 => {
                let (a, b) = j.pair(&points, n_keys);
                C::Coincident { a, b }
            }
            1 => C::PointOnCurve {
                point: j.pick(&points, n_keys),
                curve: j.pick(&curves, n_keys),
            },
            2 => C::Horizontal {
                line: j.pick(&lines, n_keys),
            },
            3 => C::Vertical {
                line: j.pick(&lines, n_keys),
            },
            4 => {
                let (a, b) = j.pair(&points, n_keys);
                C::HorizontalPoints { a, b }
            }
            5 => {
                let (a, b) = j.pair(&points, n_keys);
                C::VerticalPoints { a, b }
            }
            6 => {
                let (a, b) = j.pair(&lines, n_keys);
                C::Parallel { a, b }
            }
            7 => {
                let (a, b) = j.pair(&lines, n_keys);
                C::Perpendicular { a, b }
            }
            8 => {
                let (a, b) = if j.r.chance(0.5) {
                    j.pair(&lines, n_keys)
                } else {
                    j.pair(&circulars, n_keys)
                };
                C::Equal { a, b }
            }
            9 => C::Midpoint {
                point: j.pick(&points, n_keys),
                line: j.pick(&lines, n_keys),
            },
            10 => {
                let (a, b) = j.pair(&circulars, n_keys);
                C::Concentric { a, b }
            }
            11 => {
                let (a, b) = j.pair(&points, n_keys);
                C::Symmetric {
                    a,
                    b,
                    axis: j.pick(&lines, n_keys),
                }
            }
            12 => match j.r.int(3) {
                0 => {
                    // Joined at an end of each, made coincident first.
                    let (a, b) = j.pair(&with_ends, n_keys);
                    let end_of = |k: u32| ends.iter().find(|(id, _)| *id == k).map(|e| &e.1);
                    let (ea, eb) = (end_of(a), end_of(b));
                    let point = match (ea, eb) {
                        (Some(ea), Some(eb)) if !ea.is_empty() && !eb.is_empty() => {
                            let pa = ea[j.r.int(ea.len())];
                            let pb = eb[j.r.int(eb.len())];
                            if pa != pb {
                                extra = Some(C::Coincident { a: pa, b: pb });
                            }
                            pa
                        }
                        _ => j.pick(&points, n_keys),
                    };
                    C::Tangent {
                        a,
                        b,
                        kind: TangentKind::AtPoint { point },
                    }
                }
                1 => C::Tangent {
                    a: j.pick(&lines, n_keys),
                    b: j.pick(&roundish, n_keys),
                    kind: TangentKind::LineCircle {
                        center_left: j.r.chance(0.5),
                    },
                },
                _ => {
                    let (a, b) = j.pair(&roundish, n_keys);
                    C::Tangent {
                        a,
                        b,
                        kind: TangentKind::CircleCircle {
                            internal: j.r.chance(0.5),
                        },
                    }
                }
            },
            13 => C::Normal {
                line: j.pick(&lines, n_keys),
                curve: j.pick(&circulars, n_keys),
            },
            14 => C::Fix {
                point: j.pick(&points, n_keys),
                at: DVec2::new(j.value(-5.0, 5.0), j.value(-5.0, 5.0)),
            },
            15 => {
                let kind = [
                    DistanceKind::Aligned,
                    DistanceKind::Horizontal,
                    DistanceKind::Vertical,
                ][j.r.int(3)];
                let (a, b) = match (j.r.int(3), kind) {
                    (0, _) | (_, DistanceKind::Horizontal | DistanceKind::Vertical) => {
                        j.pair(&points, n_keys)
                    }
                    (1, _) => (j.pick(&points, n_keys), j.pick(&lines, n_keys)),
                    _ => j.pair(&lines, n_keys),
                };
                C::Distance {
                    a,
                    b,
                    kind,
                    value: j.value(0.01, 5.0),
                }
            }
            16 => C::Length {
                line: j.pick(&lines, n_keys),
                value: j.value(0.01, 5.0),
            },
            17 => {
                let (a, b) = j.pair(&lines, n_keys);
                C::Angle {
                    a,
                    b,
                    value: j.value(-PI, PI),
                }
            }
            18 => C::Radius {
                curve: j.pick(&roundish, n_keys),
                value: j.value(0.01, 5.0),
            },
            19 => C::Rho {
                curve: j.pick(&conics, n_keys),
                value: j.value(0.0, 1.0),
            },
            _ => C::Diameter {
                curve: j.pick(&circulars, n_keys),
                value: j.value(0.01, 5.0),
            },
        };
        for c in extra.into_iter().chain([constraint]) {
            let id = key(&mut j, &mut next);
            s.constraints.push(ConstraintEntry { id, constraint: c });
        }
    }
    s
}

fn geometry_ok(s: &Sketch<u32>) -> bool {
    s.entities.iter().all(|e| match e.geometry {
        Geometry::Point { at } => at.is_finite() && at.abs().max_element() <= 1e6,
        Geometry::Circle { radius, .. } => radius.is_finite() && radius > 0.0 && radius <= 1e6,
        Geometry::Conic { rho, .. } => (RHO_MIN..=1.0 - RHO_MIN).contains(&rho),
        _ => true,
    })
}

/// No line, circle or arc has collapsed: lengths and radii above 1e-10 of
/// the sketch size (the solver accepts nothing within 8e-10 of its scale,
/// which is at least that size).
fn not_collapsed(s: &Sketch<u32>) -> bool {
    let min = 1e-10 * common::size(s);
    s.entities.iter().all(|e| match e.geometry {
        Geometry::Line { .. } => {
            let (a, b) = common::line(s, e.id);
            (b - a).length() > min
        }
        Geometry::Circle { .. } | Geometry::Arc { .. } => common::circ(s, e.id).1 > min,
        Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } => {
            let (_, a, b) = common::ell(s, e.id);
            a.length() > min && b.length() > min
        }
        Geometry::Conic { .. } => {
            let ([p0, p1, p2], _) = common::cnc(s, e.id);
            let chord = p2 - p0;
            chord.length() > min && chord.perp_dot(p1 - p0).abs() / chord.length() > min
        }
        _ => true,
    })
}

fn resolved_ok(c: &ResolvedCurve) -> bool {
    match c {
        ResolvedCurve::Line { start, end } => start.is_finite() && end.is_finite() && start != end,
        ResolvedCurve::Circle { center, radius } => center.is_finite() && *radius > 0.0,
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            center.is_finite()
                && radius.is_finite()
                && *radius > 0.0
                && *start_angle > -PI
                && *start_angle <= PI
                && *sweep_angle > 0.0
                && *sweep_angle < TAU
        }
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        } => {
            center.is_finite()
                && major.is_finite()
                && minor.is_finite()
                && major.perp_dot(*minor) != 0.0
        }
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            start_param,
            sweep_param,
        } => {
            center.is_finite()
                && major.is_finite()
                && minor.is_finite()
                && major.perp_dot(*minor) != 0.0
                && start_param.is_finite()
                && sweep_param.abs() > 0.0
                && sweep_param.abs() < TAU
                // Counter-clockwise in the plane.
                && (*sweep_param > 0.0) == (major.perp_dot(*minor) > 0.0)
        }
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => {
            start.is_finite()
                && control.is_finite()
                && end.is_finite()
                && weight.is_finite()
                && *weight > 0.0
                && (*end - *start).perp_dot(*control - *start) != 0.0
        }
        ResolvedCurve::Spline {
            degree,
            knots,
            points,
        } => {
            (1..=5).contains(degree)
                && knots.len() == points.len() + degree + 1
                && knots.iter().all(|k| k.is_finite())
                && knots.windows(2).all(|w| w[0] <= w[1])
                && points.iter().all(|p| p.is_finite())
        }
    }
}

/// `resolve_curves` and `measure` on any sketch: no panic, no NaN.
fn check_queries(s: &Sketch<u32>) {
    if let Ok(curves) = resolve_curves(s) {
        for (_, _, c) in &curves {
            assert!(resolved_ok(c), "{c:?}");
        }
    }
    for c in &s.constraints {
        if let Ok(v) = measure(s, &c.constraint) {
            assert!(v.is_finite() && v >= -PI, "measure {v} for {c:?}");
            if let Constraint::Angle { .. } = c.constraint {
                assert!(v > -PI && v <= PI);
            } else {
                assert!(v >= 0.0);
            }
        }
    }
}

/// Tolerances for a claimed success: the model tolerance, or what f64
/// rounding of the coordinates allows where that is larger (directions of
/// length `short` among coordinates of magnitude `m`: angles to about
/// ε·m/short, lengths measured across a lever of up to `m` to about
/// ε·m·m/short). Returns (linear, angular).
fn success_tolerances(s: &Sketch<u32>) -> (f64, f64) {
    let mut m = 1.0f64;
    let mut short = f64::INFINITY;
    for e in &s.entities {
        match &e.geometry {
            Geometry::Point { at } => m = m.max(at.abs().max_element()),
            Geometry::Circle { radius, .. } => {
                m = m.max(*radius);
                short = short.min(*radius);
            }
            Geometry::Line { .. } => {
                let (a, b) = common::line(s, e.id);
                let d = (b - a).length();
                if d > 0.0 {
                    short = short.min(d);
                }
            }
            Geometry::Arc { .. } => short = short.min(common::circ(s, e.id).1),
            Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } => {
                let (_, a, b) = common::ell(s, e.id);
                short = short.min(a.length()).min(b.length());
            }
            Geometry::Conic { .. } => {
                let ([p0, p1, p2], _) = common::cnc(s, e.id);
                let chord = p2 - p0;
                let height = chord.perp_dot(p1 - p0).abs() / chord.length();
                short = short.min(chord.length()).min(height);
            }
            Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. } => {
                if let ResolvedCurve::Spline { points, .. } = common::resolved(s, e.id) {
                    for w in points.windows(2) {
                        let d = (w[1] - w[0]).length();
                        if d > 0.0 {
                            short = short.min(d);
                        }
                    }
                }
            }
        }
    }
    let q = (m / short.min(m)).max(1.0);
    let eps = f64::EPSILON;
    (1.5e-7 + 64.0 * eps * m * q, 1.5e-9 + 64.0 * eps * q)
}

/// A success claimed on `s` holds by the independent checker.
fn check_success(s: &Sketch<u32>) {
    let (lin, ang) = success_tolerances(s);
    for c in &s.constraints {
        let (v, angular) = violation(s, &c.constraint);
        let tol = if angular { ang } else { lin };
        assert!(v <= tol, "{c:?}: violation {v:e} > {tol:e}");
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 2000,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Random, often malformed sketches and drags: no panic, no NaN out,
    /// consistent reports, and any success holds.
    #[test]
    fn garbage_never_panics_or_lies(seed in any::<u64>(), drag_seed in any::<u64>()) {
        let input = garbage(seed);
        let mut j = Junk { r: Rng(drag_seed), h: 0.05 };
        let mut options = SolveOptions::default();
        if j.r.chance(0.3) {
            for _ in 0..1 + j.r.int(2) {
                let k = j.r.int(input.entities.len() + 2) as u32;
                let at = DVec2::new(j.value(-5.0, 5.0), j.value(-5.0, 5.0));
                options.drag.push((k, at));
            }
        }
        let out = solve(&input, &options);
        prop_assert!(out.max_residual.is_finite() && out.max_residual >= 0.0);
        prop_assert!(out.dof <= unknowns(&input));
        // (Debug text: NaN values compare equal there.)
        prop_assert_eq!(
            format!("{:?}", out.sketch.constraints),
            format!("{:?}", input.constraints)
        );
        let ids: Vec<u32> = input.constraints.iter().map(|c| c.id).collect();
        prop_assert!(out.redundant.iter().chain(&out.conflicting).all(|k| ids.contains(k)));
        prop_assert!(out.redundant.iter().all(|k| !out.conflicting.contains(k)));
        match &out.status {
            SolveStatus::Solved | SolveStatus::Redundant => {
                prop_assert!(geometry_ok(&out.sketch), "{:?}", out.sketch);
                prop_assert!(not_collapsed(&out.sketch), "{:?}", out.sketch);
                prop_assert_eq!(out.redundant.is_empty(), out.status == SolveStatus::Solved);
                prop_assert!(out.conflicting.is_empty());
                let (lin, ang) = success_tolerances(&out.sketch);
                prop_assert!(out.max_residual <= lin.max(ang), "{}", out.max_residual);
                check_success(&out.sketch);
                check_queries(&out.sketch);
                // Solving the result again changes nothing (beyond
                // rounding) and still succeeds.
                let again = solve(&out.sketch, &SolveOptions::default());
                prop_assert!(
                    !matches!(again.status, SolveStatus::Failed(_)),
                    "re-solve: {:?}", again.status
                );
                check_success(&again.sketch);
            }
            SolveStatus::Failed(e) => {
                prop_assert!(same(&out.sketch, &input));
                match e {
                    SolveError::Invalid(_) | SolveError::Degenerate(_) => {
                        check_early(&input, &out)
                    }
                    SolveError::Conflicting => prop_assert!(!out.conflicting.is_empty()),
                    SolveError::NotConverged => prop_assert!(out.conflicting.is_empty()),
                }
            }
        }
        check_queries(&input);
    }
}
