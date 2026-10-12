//! Feature previews (ADR 0013): an open dialog's feature regenerates as if
//! applied, with nothing recorded until it is accepted, where a new feature
//! would go, and never into what is exported.

mod common;

use common::*;
use opendatum_app::dialogs::{FeatureDialog, Template};
use opendatum_app::{FeaturePreview, Model, PREVIEW_FEATURE};
use opendatum_core_model::{DVec3, Feature, InsertAt};
use opendatum_features::FeatureStatus;

fn cube(side: f64) -> Feature {
    Feature::Box {
        corner: DVec3::ZERO,
        size: DVec3::splat(side),
    }
}

fn bodies(m: &Model) -> usize {
    m.bodies().bodies.len()
}

#[test]
fn a_new_feature_is_shown_without_being_recorded() {
    let mut m = Model::new();
    let ops = m.session().log().len();
    assert!(m.set_preview(Some(FeaturePreview::New(cube(0.01)))));
    assert!(
        !m.set_preview(Some(FeaturePreview::New(cube(0.01)))),
        "unchanged"
    );
    assert_eq!(m.mesh_owners(), &[PREVIEW_FEATURE]);
    assert!(close(volume(&m, PREVIEW_FEATURE), 1e-6));
    assert_eq!(m.preview_status(), Some(&FeatureStatus::Ok));
    assert_eq!(m.session().log().len(), ops, "nothing recorded");
    assert!(m.studio().features.is_empty());
    // Changing the values regenerates; clearing shows the document again.
    assert!(m.set_preview(Some(FeaturePreview::New(cube(0.02)))));
    assert!(close(volume(&m, PREVIEW_FEATURE), 8e-6));
    assert!(m.set_preview(None));
    assert!(m.meshes().is_empty());
    assert!(!m.session().can_undo());
}

#[test]
fn an_edit_preview_changes_what_is_shown_not_the_document() {
    let mut m = Model::new();
    let f = m.add_feature("Box 1", cube(0.01));
    let before = m.studio().clone();
    let undo = m.session().undo_depth();
    m.set_view(Some(f), Some(FeaturePreview::Edit(f, cube(0.03))));
    assert!(close(volume(&m, f), 27e-6));
    assert_eq!(m.studio(), &before);
    assert_eq!(m.session().undo_depth(), undo);
    m.set_view(None, None);
    assert!(close(volume(&m, f), 1e-6));
}

#[test]
fn a_new_feature_is_previewed_at_the_rollback_bar() {
    let mut m = Model::new();
    let a = m.add_feature("A", cube(0.01));
    let b = m.add_feature("B", cube(0.02));
    m.set_rollback(InsertAt::After(a));
    assert_eq!(bodies(&m), 1);
    m.set_preview(Some(FeaturePreview::New(cube(0.005))));
    let owners: Vec<_> = m.bodies().bodies.iter().map(|x| x.feature).collect();
    assert_eq!(owners, vec![a, PREVIEW_FEATURE], "B stays rolled back");
    assert_eq!(m.status(b), Some(&FeatureStatus::RolledBack));
    // Accepting inserts it there: the same bodies, now from the document.
    let Some(FeaturePreview::New(feature)) = m.preview().cloned() else {
        unreachable!()
    };
    m.clear_preview();
    let c = m.add_feature("C", feature);
    let owners: Vec<_> = m.bodies().bodies.iter().map(|x| x.feature).collect();
    assert_eq!(owners, vec![a, c]);
    assert_eq!(m.studio().index_of(c), Some(1));
}

#[test]
fn a_preview_that_does_not_regenerate_reports_why() {
    let mut m = Model::new();
    // A ring whose tube crosses the axis cannot be revolved.
    let mut d = FeatureDialog::new(Template::Ring);
    d.fields[0].text = "5 mm".into();
    d.fields[1].text = "6 mm".into();
    m.set_preview(Some(FeaturePreview::New(d.build().unwrap())));
    match m.preview_status() {
        Some(FeatureStatus::Failed(e)) => assert!(e.contains("crosses the revolve axis"), "{e}"),
        s => panic!("expected a failure, got {s:?}"),
    }
    assert!(m.meshes().is_empty(), "no partial geometry");
}

#[test]
fn export_writes_the_document_not_the_preview() {
    let dir = std::env::temp_dir().join(format!("opendatum-preview-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut m = Model::new();
    let f = m.add_feature("Box 1", cube(0.01));
    let plain = dir.join("plain.stl");
    m.export_stl(&plain).unwrap();
    m.set_preview(Some(FeaturePreview::New(cube(0.02))));
    let new = dir.join("new.stl");
    m.export_stl(&new).unwrap();
    m.set_view(Some(f), Some(FeaturePreview::Edit(f, cube(0.03))));
    let edited = dir.join("edited.stl");
    m.export_stl(&edited).unwrap();
    let bytes = |p: &std::path::Path| std::fs::read(p).unwrap();
    assert_eq!(bytes(&new), bytes(&plain));
    assert_eq!(bytes(&edited), bytes(&plain));
    std::fs::remove_dir_all(&dir).ok();
}
