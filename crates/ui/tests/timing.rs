//! Time guards for the application (docs/PERFORMANCE.md): regenerating a
//! document, which the app does after every edit. Release builds hold it to
//! a limit about ten times its measured time; a debug build only reports
//! it. CI runs this in release (`cargo test --release -p opendatum --test
//! timing`).

use std::time::{Duration, Instant};

use opendatum_app::screenshots::demo_model;

/// Release-build limit in milliseconds (measured: see docs/PERFORMANCE.md).
const REGENERATE_DEMO_MS: f64 = 400.0;

#[test]
fn regenerate_the_demo_document() {
    let (mut m, _, _) = demo_model().expect("the demo document");
    let mut fastest = Duration::MAX;
    for _ in 0..3 {
        let t0 = Instant::now();
        m.regenerate();
        fastest = fastest.min(t0.elapsed());
    }
    let ms = fastest.as_secs_f64() * 1e3;
    println!("regenerate the demo document: {ms:.2} ms (limit {REGENERATE_DEMO_MS} ms)");
    if !cfg!(debug_assertions) {
        assert!(
            ms <= REGENERATE_DEMO_MS,
            "regenerating the demo document took {ms:.2} ms, over its limit of \
             {REGENERATE_DEMO_MS} ms (docs/PERFORMANCE.md)"
        );
    }
}
