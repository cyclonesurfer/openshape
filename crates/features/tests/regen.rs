//! Regeneration through the op log: failures are isolated, suppression and
//! undo change the produced bodies.

use openshape_core_model::{
    CurveDef, DVec2, DVec3, Feature, FeatureId, InsertAt, ModelOp, PlaneRef, ProfileDef, ReplicaId,
    Session,
};
use openshape_features::{regenerate, FeatureStatus};
use openshape_kernel::{Kernel, OsKernel};

fn insert(s: &mut Session, name: &str, f: Feature) -> FeatureId {
    let el = s.first_studio().unwrap();
    FeatureId(s.apply(ModelOp::InsertFeature {
        element: el,
        at: InsertAt::End,
        name: name.into(),
        feature: f,
    }))
}

#[test]
fn failed_feature_does_not_stop_later_features() {
    let k = OsKernel::default();
    let mut s = Session::new_document(ReplicaId(1));
    let a = insert(
        &mut s,
        "Box",
        Feature::Box {
            corner: DVec3::ZERO,
            size: DVec3::new(0.01, 0.02, 0.03),
        },
    );
    let bad = insert(
        &mut s,
        "Open profile",
        Feature::Extrude {
            profile: ProfileDef {
                plane: PlaneRef::Top,
                curves: vec![CurveDef::Line {
                    start: DVec2::ZERO,
                    end: DVec2::X,
                }],
            },
            distance: 0.01,
            symmetric: false,
        },
    );
    let c = insert(
        &mut s,
        "Revolve",
        Feature::Revolve {
            profile: ProfileDef {
                plane: PlaneRef::Front,
                curves: vec![CurveDef::Circle {
                    center: DVec2::new(0.05, 0.0),
                    radius: 0.01,
                }],
            },
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: std::f64::consts::TAU,
        },
    );
    let studio = &s.document().elements[0].studio;
    let r = regenerate(&k, studio);
    assert_eq!(r.status(a), Some(&FeatureStatus::Ok));
    assert!(matches!(r.status(bad), Some(FeatureStatus::Failed(m)) if m.contains("open")));
    assert_eq!(r.status(c), Some(&FeatureStatus::Ok));
    assert_eq!(r.bodies.len(), 2);
    let v = k.mass_properties(&r.bodies[0].body, 1e-3).unwrap().volume;
    assert!((v - 0.01 * 0.02 * 0.03).abs() < 1e-12);

    // Suppress the box: one body left. Undo: two again.
    let el = s.first_studio().unwrap();
    s.apply(ModelOp::SetSuppressed {
        element: el,
        feature: a,
        suppressed: true,
    });
    let r = regenerate(&k, &s.document().elements[0].studio);
    assert_eq!(r.status(a), Some(&FeatureStatus::Suppressed));
    assert_eq!(r.bodies.len(), 1);
    assert!(s.undo());
    let r = regenerate(&k, &s.document().elements[0].studio);
    assert_eq!(r.bodies.len(), 2);
}
