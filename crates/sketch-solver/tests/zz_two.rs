use openshape_sketch_solver::{solve, SolveOptions, Sketch};

#[test]
fn two() {
    let j = std::fs::read_to_string(std::env::var("SK").unwrap()).unwrap();
    let s: Sketch<u32> = serde_json::from_str(&j).unwrap();
    let r = solve(&s, &SolveOptions::default());
    println!("{:?} it {} res {:e} dof {}", r.status, r.iterations, r.max_residual, r.dof);
}
