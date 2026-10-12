//! End-to-end: op log → document → features → kernel → tessellation → wgpu
//! offscreen render, with orbit, undo/redo, save/open and STL export checks.
//! Runs on any Vulkan/GL adapter, including Mesa's software lavapipe in CI.

#[test]
fn headless_demo_passes_all_checks() {
    let dir = std::env::temp_dir().join(format!("openshape-demo-{}", std::process::id()));
    let report = openshape_app::demo::run(&dir).expect("demo runs (needs a Vulkan or GL adapter)");
    for (name, ok) in &report.checks {
        println!("[{}] {name}", if *ok { "PASS" } else { "FAIL" });
    }
    assert!(
        report.ok(),
        "failed checks: {:?}",
        report.checks.iter().filter(|c| !c.1).collect::<Vec<_>>()
    );
    assert_eq!(report.images.len(), 6);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn model_undo_redo_and_display_meshes() {
    use openshape_app::dialogs::{FeatureDialog, Template};
    let mut m = openshape_app::Model::new();
    assert!(m.meshes().is_empty());
    m.add_feature("Box 1", FeatureDialog::new(Template::Box).build().unwrap());
    m.add_feature(
        "Ring 1",
        FeatureDialog::new(Template::Ring).build().unwrap(),
    );
    assert_eq!(m.meshes().len(), 2);
    assert!(m.meshes().iter().all(|d| !d.indices.is_empty()));
    // The box shows its 12 edges; a full torus has no creases, only the
    // seams where its faces are split, which are not drawn.
    assert_eq!(m.meshes()[0].edge_segment_count(), 12);
    assert!(m.meshes()[1].edge_lines.is_empty());
    assert!(m.undo());
    assert_eq!(m.meshes().len(), 1);
    assert!(m.redo());
    assert_eq!(m.meshes().len(), 2);
    assert!(!m.redo());
    assert_eq!(m.next_name("Revolve"), "Revolve 2");
}

#[test]
fn bench_times_every_frame() {
    let r = openshape_app::demo::bench(5, 320, 240).expect("bench runs (needs an adapter)");
    assert_eq!(r.frame_ms.len(), 5);
    assert!(r.frame_ms.iter().all(|t| t.is_finite() && *t > 0.0));
    assert!(r.triangles > 1000 && r.edge_segments > 0);
    assert!(r.quantile(0.0) <= r.quantile(0.5) && r.quantile(0.5) <= r.quantile(1.0));
    assert!(r.summary().contains("frame time ms"));
}
