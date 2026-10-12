//! Timing smoke tests (design note 0005 §8): solves of 50- and 300-entity
//! sketches are timed and reported (run with `--nocapture` to see the
//! numbers). Limits are generous in debug builds and tight in release
//! builds, where they guard against complexity regressions.

mod common;

use std::time::{Duration, Instant};

use common::gen::{generate, perturb};
use common::{assert_satisfied, B};
use opendatum_sketch_solver::{
    solve, Constraint, DVec2, Sketch, Solution, SolveOptions, SolveStatus,
};

/// Generated sketches (many independent motifs, many components) with at
/// least `target` entities.
fn generated(target: usize, seed: u64) -> Sketch<u32> {
    let mut motifs = 1;
    loop {
        let g = generate(seed, motifs, 1.0, DVec2::new(0.3, -0.2));
        if g.sketch.entities.len() >= target {
            return g.sketch;
        }
        motifs += 1;
    }
}

/// One large component: a chain of `n` lines (3n entities) joined end to
/// end, each with a length and an angle to the previous one, the first
/// fixed and horizontal: fully constrained.
fn chain(n: usize) -> Sketch<u32> {
    chain_with(n, true)
}

/// [`chain`], with or without the angles (without: `n − 1` freedoms).
fn chain_with(n: usize, angles: bool) -> Sketch<u32> {
    let mut b = B::new();
    let mut at = DVec2::ZERO;
    let mut heading = 0.0f64;
    let mut prev: Option<(u32, u32)> = None;
    for i in 0..n {
        let len = 1.0 + 0.5 * ((i * 7 % 5) as f64);
        let turn = if i == 0 {
            0.0
        } else {
            0.3 + 0.9 * ((i * 3 % 4) as f64) / 3.0
        } * if i % 2 == 0 { 1.0 } else { -1.0 };
        heading += turn;
        let end = at + DVec2::new(heading.cos(), heading.sin()) * len;
        let (l, s, e) = b.line(at.x, at.y, end.x, end.y);
        b.c(Constraint::Length {
            line: l,
            value: len,
        });
        match prev {
            None => {
                b.fix(s);
                b.c(Constraint::Horizontal { line: l });
            }
            Some((pl, pe)) => {
                b.coincident(pe, s);
                if angles {
                    b.c(Constraint::Angle {
                        a: pl,
                        b: l,
                        value: turn,
                    });
                }
            }
        }
        prev = Some((l, e));
        at = end;
    }
    b.build()
}

/// Solves `s` (perturbed from its solution) a few times; returns the result
/// and the fastest time.
fn timed(s: &Sketch<u32>) -> (Solution<u32>, Duration) {
    let input = perturb(s, 7, 0.01, 1.0);
    let mut best = Duration::MAX;
    let mut out = None;
    for _ in 0..3 {
        let t0 = Instant::now();
        let r = solve(&input, &SolveOptions::default());
        best = best.min(t0.elapsed());
        out = Some(r);
    }
    (out.expect("solved"), best)
}

fn check(label: &str, s: &Sketch<u32>, release_limit: Duration) {
    let (r, t) = timed(s);
    println!(
        "{label}: {} entities, {} constraints: {:?} in {:.3} ms, {} iterations, dof {}",
        s.entities.len(),
        s.constraints.len(),
        r.status,
        t.as_secs_f64() * 1e3,
        r.iterations,
        r.dof,
    );
    assert!(
        matches!(r.status, SolveStatus::Solved | SolveStatus::Redundant),
        "{label}: {:?}",
        r.status
    );
    assert_satisfied(&r.sketch);
    // Debug builds: only a sanity bound. Release builds: the real limit.
    let limit = if cfg!(debug_assertions) {
        Duration::from_secs(60)
    } else {
        release_limit
    };
    assert!(t <= limit, "{label}: {t:?} exceeds {limit:?}");
}

#[test]
fn fifty_entities() {
    let g = generated(50, 11);
    assert!(g.entities.len() >= 50);
    // Design note §1: well within a 60 Hz frame.
    check("generated ~50", &g, Duration::from_millis(8));
    let c = chain(17);
    assert_eq!(c.entities.len(), 51);
    check("chain of 17 lines", &c, Duration::from_millis(8));
}

#[test]
fn three_hundred_entities() {
    let g = generated(300, 23);
    assert!(g.entities.len() >= 300);
    // Design note §1: within 50 ms.
    check("generated ~300", &g, Duration::from_millis(50));
    let c = chain(100);
    assert_eq!(c.entities.len(), 300);
    let (r, _) = timed(&c);
    assert_eq!(r.dof, 0, "the chain is fully constrained");
    check("chain of 100 lines", &c, Duration::from_millis(50));
    // Under-constrained (99 freedoms), and dragged by its free end as in
    // the UI: one solve per pointer move from the last result.
    let loose = chain_with(100, false);
    check(
        "chain of 100 lines, no angles",
        &loose,
        Duration::from_millis(50),
    );
    let end = loose.entities[loose.entities.len() - 2].id;
    let mut s = loose;
    let mut worst = Duration::ZERO;
    let start = s.point(end).unwrap();
    for i in 1..=20 {
        let target = start + DVec2::new(-0.5, 0.8) * i as f64;
        let options = SolveOptions {
            drag: vec![(end, target)],
            ..SolveOptions::default()
        };
        let t0 = Instant::now();
        let r = solve(&s, &options);
        worst = worst.max(t0.elapsed());
        assert!(
            matches!(r.status, SolveStatus::Solved | SolveStatus::Redundant),
            "drag step {i}: {:?}",
            r.status
        );
        assert_satisfied(&r.sketch);
        s = r.sketch;
    }
    println!(
        "dragging the loose chain: slowest of 20 steps {:.3} ms",
        worst.as_secs_f64() * 1e3
    );
    let limit = if cfg!(debug_assertions) {
        Duration::from_secs(60)
    } else {
        Duration::from_millis(50)
    };
    assert!(worst <= limit, "drag: {worst:?}");
}
