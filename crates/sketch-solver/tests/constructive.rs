//! Constructive property tests (design note 0005 §8): random geometry with
//! constraints that hold on it and dimensions measured from it, perturbed,
//! then solved. Every constraint must be satisfied afterwards, at sketch
//! scales from 1e-4 to 1e3.

mod common;

use common::gen::{generate, perturb};
use common::{assert_satisfied, bits, size, worst};
use openshape_sketch_solver::{
    solve, DVec2, EntityState, Geometry, Sketch, SolveOptions, SolveStatus,
};
use proptest::prelude::*;

fn finite(s: &Sketch<u32>) -> bool {
    s.entities.iter().all(|e| match e.geometry {
        Geometry::Point { at } => at.is_finite(),
        Geometry::Circle { radius, .. } => radius.is_finite() && radius > 0.0,
        _ => true,
    })
}

fn largest_move(a: &Sketch<u32>, b: &Sketch<u32>) -> f64 {
    let mut m = 0.0f64;
    for (x, y) in a.entities.iter().zip(&b.entities) {
        match (&x.geometry, &y.geometry) {
            (Geometry::Point { at: p }, Geometry::Point { at: q }) => m = m.max((*p - *q).length()),
            (Geometry::Circle { radius: r, .. }, Geometry::Circle { radius: t, .. }) => {
                m = m.max((r - t).abs())
            }
            _ => {}
        }
    }
    m
}

/// Checks the parts of a solution that hold whatever the outcome.
fn check_common(input: &Sketch<u32>, s: &openshape_sketch_solver::Solution<u32>) {
    assert!(finite(&s.sketch), "non-finite geometry");
    assert!(s.max_residual.is_finite(), "max_residual {}", s.max_residual);
    assert_eq!(s.entities.len(), input.entities.len(), "entity states");
    for e in &input.entities {
        assert!(s.entities.contains_key(&e.id));
    }
    assert_eq!(s.sketch.constraints, input.constraints, "constraints are not changed");
    let unknowns: usize = input
        .entities
        .iter()
        .map(|e| match e.geometry {
            Geometry::Point { .. } => 2,
            Geometry::Circle { .. } => 1,
            _ => 0,
        })
        .sum();
    assert!(s.dof <= unknowns);
}

proptest! {
    #![proptest_config(ProptestConfig {
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// The generator itself: its sketches satisfy their constraints.
    #[test]
    fn generated_sketches_hold(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::new(ox, oy));
        let (w, what) = worst(&g.sketch, size(&g.sketch));
        prop_assert!(w <= 1e-11, "generator: violation {:e} at {}", w, what);
    }

    /// Solving an already satisfied sketch leaves it (essentially) where
    /// it is.
    #[test]
    fn satisfied_sketches_stay_put(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::new(ox, oy));
        let s = solve(&g.sketch, &SolveOptions::default());
        check_common(&g.sketch, &s);
        prop_assert!(
            matches!(s.status, SolveStatus::Solved | SolveStatus::Redundant),
            "status {:?}, conflicting {:?}, residual {:e}",
            s.status, s.conflicting, s.max_residual
        );
        let moved = largest_move(&g.sketch, &s.sketch) / size(&g.sketch).max(1e-300);
        prop_assert!(moved <= 1e-9, "moved by {:e} (relative)", moved);
        assert_satisfied(&s.sketch);
    }

    /// Perturbed by up to 2 % of the sketch size, every sketch solves and
    /// every constraint holds.
    #[test]
    fn perturbed_sketches_solve(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        ox in -3.0f64..3.0,
        oy in -3.0f64..3.0,
        amount in 0.0f64..0.02,
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::new(ox, oy));
        let input = perturb(&g.sketch, seed, amount, scale);
        let s = solve(&input, &SolveOptions::default());
        check_common(&input, &s);
        prop_assert!(
            matches!(s.status, SolveStatus::Solved | SolveStatus::Redundant),
            "status {:?}, conflicting {:?}, residual {:e}",
            s.status, s.conflicting, s.max_residual
        );
        prop_assert_eq!(
            s.status == SolveStatus::Redundant,
            !s.redundant.is_empty()
        );
        prop_assert!(s.conflicting.is_empty());
        prop_assert!(s.max_residual <= 1e-7, "max_residual {:e}", s.max_residual);
        assert_satisfied(&s.sketch);
        if s.dof == 0 {
            prop_assert!(s.entities.values().all(|&st| st == EntityState::Constrained));
        }
    }

    /// Larger perturbations may fail, but a success always satisfies every
    /// constraint and a failure returns the input unchanged.
    #[test]
    fn large_perturbations_never_lie(
        seed in any::<u64>(),
        motifs in 1usize..9,
        e in -4i32..=3,
        amount in 0.02f64..0.4,
    ) {
        let scale = 10f64.powi(e);
        let g = generate(seed, motifs, scale, DVec2::ZERO);
        let input = perturb(&g.sketch, seed, amount, scale);
        let s = solve(&input, &SolveOptions::default());
        check_common(&input, &s);
        match &s.status {
            SolveStatus::Solved | SolveStatus::Redundant => assert_satisfied(&s.sketch),
            SolveStatus::Failed(_) => prop_assert_eq!(bits(&s.sketch), bits(&input)),
        }
    }
}
