//! Determinism (design note 0005 §8): the same input gives bit-identical
//! output, every time, on any thread, whatever the key values.

mod common;

use common::gen::{generate, perturb};
use common::{bits, B};
use openshape_sketch_solver::{
    solve, Constraint, ConstraintEntry, DVec2, Entity, Geometry, Sketch, SketchKey, Solution,
    SolveOptions, TangentKind,
};
use proptest::prelude::*;

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// Asserts two solutions are identical, floating-point values to the bit.
fn identical<K: SketchKey>(a: &Solution<K>, b: &Solution<K>) {
    assert_eq!(a.status, b.status);
    assert_eq!(a.dof, b.dof);
    assert_eq!(a.entities, b.entities);
    assert_eq!(a.redundant, b.redundant);
    assert_eq!(a.conflicting, b.conflicting);
    assert_eq!(a.iterations, b.iterations);
    assert_eq!(a.max_residual.to_bits(), b.max_residual.to_bits());
    assert_eq!(geometry_bits(&a.sketch), geometry_bits(&b.sketch));
    assert_eq!(a.sketch.constraints, b.sketch.constraints);
}

fn geometry_bits<K>(s: &Sketch<K>) -> Vec<u64> {
    let mut out = Vec::new();
    for e in &s.entities {
        match e.geometry {
            Geometry::Point { at } => out.extend([at.x.to_bits(), at.y.to_bits()]),
            Geometry::Circle { radius, .. } => out.push(radius.to_bits()),
            _ => {}
        }
    }
    out
}

/// The same sketch with every key `k` replaced by `f(k)`.
fn relabel<F: Fn(u32) -> u64>(s: &Sketch<u32>, f: F) -> Sketch<u64> {
    let geometry = |g: &Geometry<u32>| match g.clone() {
        Geometry::Point { at } => Geometry::Point { at },
        Geometry::Line { start, end } => Geometry::Line {
            start: f(start),
            end: f(end),
        },
        Geometry::Circle { center, radius } => Geometry::Circle {
            center: f(center),
            radius,
        },
        Geometry::Arc { center, start, end } => Geometry::Arc {
            center: f(center),
            start: f(start),
            end: f(end),
        },
        Geometry::ControlSpline { control, degree } => Geometry::ControlSpline {
            control: control.into_iter().map(&f).collect(),
            degree,
        },
        Geometry::InterpolatedSpline { through, params } => Geometry::InterpolatedSpline {
            through: through.into_iter().map(&f).collect(),
            params,
        },
    };
    use Constraint as C;
    let constraint = |c: &Constraint<u32>| match c.clone() {
        C::Coincident { a, b } => C::Coincident { a: f(a), b: f(b) },
        C::PointOnCurve { point, curve } => C::PointOnCurve {
            point: f(point),
            curve: f(curve),
        },
        C::Horizontal { line } => C::Horizontal { line: f(line) },
        C::Vertical { line } => C::Vertical { line: f(line) },
        C::HorizontalPoints { a, b } => C::HorizontalPoints { a: f(a), b: f(b) },
        C::VerticalPoints { a, b } => C::VerticalPoints { a: f(a), b: f(b) },
        C::Parallel { a, b } => C::Parallel { a: f(a), b: f(b) },
        C::Perpendicular { a, b } => C::Perpendicular { a: f(a), b: f(b) },
        C::Equal { a, b } => C::Equal { a: f(a), b: f(b) },
        C::Midpoint { point, line } => C::Midpoint {
            point: f(point),
            line: f(line),
        },
        C::Concentric { a, b } => C::Concentric { a: f(a), b: f(b) },
        C::Symmetric { a, b, axis } => C::Symmetric {
            a: f(a),
            b: f(b),
            axis: f(axis),
        },
        C::Tangent { a, b, kind } => C::Tangent {
            a: f(a),
            b: f(b),
            kind: match kind {
                TangentKind::AtPoint { point } => TangentKind::AtPoint { point: f(point) },
                TangentKind::LineCircle { center_left } => TangentKind::LineCircle { center_left },
                TangentKind::CircleCircle { internal } => TangentKind::CircleCircle { internal },
            },
        },
        C::Normal { line, curve } => C::Normal {
            line: f(line),
            curve: f(curve),
        },
        C::Fix { point, at } => C::Fix {
            point: f(point),
            at,
        },
        C::Distance { a, b, kind, value } => C::Distance {
            a: f(a),
            b: f(b),
            kind,
            value,
        },
        C::Length { line, value } => C::Length {
            line: f(line),
            value,
        },
        C::Angle { a, b, value } => C::Angle {
            a: f(a),
            b: f(b),
            value,
        },
        C::Radius { curve, value } => C::Radius {
            curve: f(curve),
            value,
        },
        C::Diameter { curve, value } => C::Diameter {
            curve: f(curve),
            value,
        },
        C::Mirror { axis, pairs } => C::Mirror {
            axis: f(axis),
            pairs: pairs.into_iter().map(|(s, i)| (f(s), f(i))).collect(),
        },
        C::Offset {
            distance,
            pieces,
            smooth,
        } => C::Offset {
            distance,
            pieces: pieces
                .into_iter()
                .map(|(s, i, l)| (f(s), f(i), l))
                .collect(),
            smooth: smooth.into_iter().map(&f).collect(),
        },
        C::LinearPattern { direction, copies } => C::LinearPattern {
            direction: f(direction),
            copies: copies
                .into_iter()
                .map(|(s, i, k)| (f(s), f(i), k))
                .collect(),
        },
        C::CircularPattern {
            centre,
            angle,
            copies,
        } => C::CircularPattern {
            centre: f(centre),
            angle,
            copies: copies
                .into_iter()
                .map(|(s, i, k)| (f(s), f(i), k))
                .collect(),
        },
    };
    Sketch {
        entities: s
            .entities
            .iter()
            .map(|e| Entity {
                id: f(e.id),
                geometry: geometry(&e.geometry),
                construction: e.construction,
            })
            .collect(),
        constraints: s
            .constraints
            .iter()
            .map(|c| ConstraintEntry {
                id: f(c.id),
                constraint: constraint(&c.constraint),
            })
            .collect(),
    }
}

/// A sketch with a bit of everything, off its solution, with a redundant
/// constraint.
fn mixed() -> Sketch<u32> {
    let mut b = B::new();
    let (l0, a0, a1) = b.line(0.0, 0.0, 3.0, 0.1);
    let (l1, b0, b1) = b.line(3.1, 0.0, 3.0, 2.0);
    let (arc, _, s, _) = b.arc(DVec2::new(1.0, 2.0), 1.0, 0.2, 2.9);
    let (circle, cc) = b.circle(-2.0, 1.0, 0.7);
    let ctrl: Vec<u32> = [(-3.0, -1.0), (-2.0, -2.0), (-1.0, -0.5), (0.0, -1.5)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    let sp = b.entity(Geometry::ControlSpline {
        control: ctrl.clone(),
        degree: 3,
    });
    b.coincident(a1, b0);
    b.coincident(b1, s);
    b.c(Constraint::Horizontal { line: l0 });
    b.c(Constraint::Perpendicular { a: l0, b: l1 });
    b.c(Constraint::Vertical { line: l1 });
    b.c(Constraint::Length {
        line: l0,
        value: 3.2,
    });
    b.c(Constraint::Tangent {
        a: arc,
        b: circle,
        kind: TangentKind::CircleCircle { internal: false },
    });
    b.c(Constraint::PointOnCurve {
        point: cc,
        curve: sp,
    });
    b.coincident(ctrl[3], a0);
    b.fix(ctrl[0]);
    b.c(Constraint::Radius {
        curve: arc,
        value: 1.1,
    });
    b.build()
}

#[test]
fn solving_twice_gives_identical_bits() {
    let s = mixed();
    let a = solve(&s, &opts());
    assert!(a.iterations > 0);
    let b = solve(&s.clone(), &opts());
    identical(&a, &b);
    assert_eq!(a, b);
}

#[test]
fn solving_on_other_threads_gives_identical_bits() {
    let s = mixed();
    let here = solve(&s, &opts());
    let handles: Vec<_> = (0..4)
        .map(|_| {
            let s = s.clone();
            std::thread::spawn(move || solve(&s, &SolveOptions::default()))
        })
        .collect();
    for h in handles {
        identical(&here, &h.join().unwrap());
    }
}

#[test]
fn key_values_do_not_change_the_result() {
    let s = mixed();
    let a = solve(&s, &opts());
    // An order-reversing relabelling: entity lookups change, entity and
    // constraint order do not.
    let f = |k: u32| 1_000_000 - u64::from(k);
    let t = relabel(&s, f);
    let b = solve(&t, &SolveOptions::default());
    assert_eq!(a.status, b.status);
    assert_eq!(a.dof, b.dof);
    assert_eq!(a.iterations, b.iterations);
    assert_eq!(a.max_residual.to_bits(), b.max_residual.to_bits());
    assert_eq!(geometry_bits(&a.sketch), geometry_bits(&b.sketch));
    let mapped: Vec<u64> = a.redundant.iter().map(|&k| f(k)).collect();
    assert_eq!(mapped, b.redundant);
    for (k, st) in &a.entities {
        assert_eq!(b.entities[&f(*k)], *st);
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Determinism over generated sketches, including failures (large
    /// perturbations).
    #[test]
    fn generated_solves_are_deterministic(
        seed in any::<u64>(),
        motifs in 1usize..8,
        e in -4i32..=3,
        amount in 0.0f64..0.4,
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::new(0.3, -0.2));
        let s = perturb(&g.sketch, seed ^ 0x9e37, amount, scale);
        let a = solve(&s, &opts());
        let b = solve(&s, &opts());
        identical(&a, &b);
        prop_assert_eq!(bits(&a.sketch), bits(&b.sketch));
        let f = |k: u32| u64::from(k) * 7 + 3;
        let c = solve(&relabel(&s, f), &SolveOptions::default());
        prop_assert_eq!(geometry_bits(&a.sketch), geometry_bits(&c.sketch));
        prop_assert_eq!(a.max_residual.to_bits(), c.max_residual.to_bits());
    }
}
