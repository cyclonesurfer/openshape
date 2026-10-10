mod common;
use common::gen::{generate, perturb};
use common::violation;
use openshape_sketch_solver::{solve, DVec2, SolveOptions, SolveStatus, Sketch};

fn fails(s: &Sketch<u32>) -> bool {
    matches!(solve(s, &SolveOptions::default()).status, SolveStatus::Failed(openshape_sketch_solver::SolveError::NotConverged | openshape_sketch_solver::SolveError::Conflicting))
}

#[test]
fn one() {
    let seed: u64 = std::env::var("SEED").unwrap().parse().unwrap();
    let motifs = 1 + (seed as usize % 8);
    let amt: f64 = std::env::var("AMT").ok().and_then(|v| v.parse().ok()).unwrap_or(0.02);
    let g = generate(seed * 7919 + 13, motifs, 1.0, DVec2::ZERO);
    let mut input = perturb(&g.sketch, seed, amt, 1.0);
    assert!(fails(&input));
    // Minimise constraints.
    let mut i = 0;
    while i < input.constraints.len() {
        let mut t = input.clone();
        t.constraints.remove(i);
        if fails(&t) { input = t; } else { i += 1; }
    }
    // Remove unreferenced entities (points/curves not used).
    loop {
        let mut changed = false;
        for k in 0..input.entities.len() {
            let id = input.entities[k].id;
            let used = input.constraints.iter().any(|c| c.constraint.referenced().contains(&id))
                || input.entities.iter().any(|e| e.geometry.referenced_points().contains(&id));
            if !used {
                let mut t = input.clone();
                t.entities.remove(k);
                if fails(&t) { input = t; changed = true; break; }
            }
        }
        if !changed { break; }
    }
    for e in &input.entities { println!("{:?}", e); }
    for c in &input.constraints { println!("{:?}  viol {:?}", c, violation(&input, &c.constraint)); }
    let s = solve(&input, &SolveOptions::default());
    println!("{:?} it {} res {:e} dof {} red {:?}", s.status, s.iterations, s.max_residual, s.dof, s.redundant);
    println!("{}", serde_json::to_string(&input).unwrap());
}
