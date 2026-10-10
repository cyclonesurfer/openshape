mod common;
use common::gen::{generate, perturb};
use openshape_sketch_solver::{solve, DVec2, SolveOptions, SolveStatus};
use std::collections::BTreeMap;

#[test]
fn stats() {
    let mut kinds: BTreeMap<String, usize> = BTreeMap::new();
    let mut status: BTreeMap<String, usize> = BTreeMap::new();
    let (mut ne, mut nc, mut it, mut dof0, mut red) = (0, 0, 0, 0, 0);
    let n = 2000;
    let t = std::time::Instant::now();
    for seed in 0..n {
        let g = generate(seed * 7919 + 13, 1 + (seed as usize % 8), 1.0, DVec2::ZERO);
        ne += g.sketch.entities.len();
        nc += g.sketch.constraints.len();
        for c in &g.sketch.constraints {
            let k = format!("{:?}", c.constraint);
            let k = k.split(|ch: char| !ch.is_alphanumeric()).next().unwrap().to_string();
            *kinds.entry(k).or_default() += 1;
        }
        let amt: f64 = std::env::var("AMT").ok().and_then(|v| v.parse().ok()).unwrap_or(0.02);
        let input = perturb(&g.sketch, seed, amt, 1.0);
        let s = solve(&input, &SolveOptions::default());
        *status.entry(format!("{:?}", s.status)).or_default() += 1;
        if matches!(s.status, SolveStatus::Failed(_)) && std::env::var("SHOW").is_ok() {
            println!("seed {seed} motifs {} status {:?} it {} res {:e} conflicting {:?}", 1 + (seed as usize % 8), s.status, s.iterations, s.max_residual, s.conflicting);
            for c in &g.sketch.constraints { if s.conflicting.contains(&c.id) { println!("   {:?}", c); } }
        }
        it += s.iterations;
        if s.dof == 0 { dof0 += 1; }
        if !s.redundant.is_empty() { red += 1; }
    }
    println!("time {:?}", t.elapsed());
    println!("avg entities {} constraints {} iterations {}", ne as f64 / n as f64, nc as f64 / n as f64, it as f64 / n as f64);
    println!("dof0 {dof0} redundant {red}");
    println!("{kinds:#?}\n{status:#?}");
    let _ = SolveStatus::Solved;
}
