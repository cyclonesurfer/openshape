//! Sketch edits through the op log (ADR 0011): id resolution, atomic
//! batches, cascading deletes, undo, merging, and serialization.

use opendatum_core_model::sketch::{Constraint, ConstraintEntry, Entity, Geometry};
use opendatum_core_model::{
    DVec2, Feature, FeatureId, InsertAt, ModelOp, PlaneRef, ReplicaId, Session, SketchData,
    SketchDef, SketchEdit, SketchId,
};

fn point(index: u32, x: f64, y: f64) -> SketchEdit {
    SketchEdit::AddEntity(Entity {
        id: SketchId::new_in_batch(index),
        geometry: Geometry::Point {
            at: DVec2::new(x, y),
        },
        construction: false,
    })
}

fn line(index: u32, a: SketchId, b: SketchId) -> SketchEdit {
    SketchEdit::AddEntity(Entity {
        id: SketchId::new_in_batch(index),
        geometry: Geometry::Line { start: a, end: b },
        construction: false,
    })
}

fn constraint(index: u32, c: Constraint<SketchId>) -> SketchEdit {
    SketchEdit::AddConstraint(ConstraintEntry {
        id: SketchId::new_in_batch(index),
        constraint: c,
    })
}

fn new_sketch(s: &mut Session) -> FeatureId {
    let element = s.first_studio().unwrap();
    FeatureId(s.apply(ModelOp::InsertFeature {
        element,
        at: InsertAt::End,
        name: "Sketch 1".into(),
        feature: Feature::Sketch(SketchDef::new(PlaneRef::Top, SketchData::default())),
    }))
}

fn edit(s: &mut Session, feature: FeatureId, edits: Vec<SketchEdit>) -> opendatum_core_model::OpId {
    let element = s.first_studio().unwrap();
    s.apply(ModelOp::EditSketch {
        element,
        feature,
        edits,
    })
}

fn sketch(s: &Session, feature: FeatureId) -> SketchDef {
    let el = s.document().element(s.first_studio().unwrap()).unwrap();
    let f = el.studio.features.iter().find(|f| f.id == feature).unwrap();
    match &f.feature {
        Feature::Sketch(d) => d.clone(),
        _ => panic!("not a sketch"),
    }
}

/// A line with a horizontal constraint and a length, in one batch.
fn draw_line(s: &mut Session, f: FeatureId) -> opendatum_core_model::OpId {
    let n = SketchId::new_in_batch;
    edit(
        s,
        f,
        vec![
            point(0, 0.0, 0.0),
            point(1, 0.01, 0.0),
            line(2, n(0), n(1)),
            constraint(3, Constraint::Horizontal { line: n(2) }),
            constraint(
                4,
                Constraint::Length {
                    line: n(2),
                    value: 0.01,
                },
            ),
        ],
    )
}

#[test]
fn batch_ids_resolve_to_the_applying_op() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_line(&mut s, f);
    let d = sketch(&s, f);
    assert!(s.document().conflicts.is_empty());
    assert_eq!(d.sketch.entities.len(), 3);
    assert_eq!(d.sketch.constraints.len(), 2);
    let id = |index| SketchId { op, index };
    assert_eq!(
        d.sketch.entity(id(2)).unwrap().geometry,
        Geometry::Line {
            start: id(0),
            end: id(1)
        }
    );
    assert_eq!(
        d.sketch.constraint(id(3)).unwrap().constraint,
        Constraint::Horizontal { line: id(2) }
    );
    assert!(d
        .sketch
        .entities
        .iter()
        .chain(std::iter::empty())
        .all(|e| !e.id.is_new_in_batch()));
}

#[test]
fn a_failing_edit_leaves_the_sketch_unchanged() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    draw_line(&mut s, f);
    let before = sketch(&s, f);
    let n = SketchId::new_in_batch;
    // The second edit references a point that doesn't exist.
    let bad = edit(
        &mut s,
        f,
        vec![
            point(0, 1.0, 1.0),
            line(
                1,
                n(0),
                SketchId {
                    op: opendatum_core_model::OpId::NULL,
                    index: 99,
                },
            ),
        ],
    );
    assert_eq!(sketch(&s, f), before);
    let conflicts = &s.document().conflicts;
    assert_eq!(conflicts.len(), 1);
    assert_eq!(conflicts[0].op, bad);
    assert!(conflicts[0].reason.contains("sketch edit 2"));
}

#[test]
fn invalid_values_are_rejected() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_line(&mut s, f);
    let len = SketchId { op, index: 4 };
    let p0 = SketchId { op, index: 0 };
    let ln = SketchId { op, index: 2 };
    for e in [
        SketchEdit::SetValue {
            constraint: len,
            value: f64::NAN,
        },
        SketchEdit::SetValue {
            constraint: SketchId { op, index: 3 },
            value: 1.0,
        },
        SketchEdit::MovePoints {
            points: vec![(p0, DVec2::new(f64::INFINITY, 0.0))],
        },
        SketchEdit::MovePoints {
            points: vec![(ln, DVec2::ZERO)],
        },
        SketchEdit::SetRadius {
            entity: ln,
            radius: 1.0,
        },
        point(0, f64::NAN, 0.0),
        line(0, p0, p0),
    ] {
        let before = sketch(&s, f);
        let n_conflicts = s.document().conflicts.len();
        edit(&mut s, f, vec![e.clone()]);
        assert_eq!(sketch(&s, f), before, "{e:?}");
        assert_eq!(s.document().conflicts.len(), n_conflicts + 1, "{e:?}");
    }
}

/// A conic on three new points with ρ `rho`, and its ρ dimension.
fn draw_conic(s: &mut Session, f: FeatureId, rho: f64) -> opendatum_core_model::OpId {
    let n = SketchId::new_in_batch;
    edit(
        s,
        f,
        vec![
            point(0, 0.0, 0.0),
            point(1, 0.003, 0.008),
            point(2, 0.01, 0.001),
            SketchEdit::AddEntity(Entity {
                id: n(3),
                geometry: Geometry::Conic {
                    start: n(0),
                    control: n(1),
                    end: n(2),
                    rho,
                },
                construction: false,
            }),
            constraint(
                4,
                Constraint::Rho {
                    curve: n(3),
                    value: rho,
                },
            ),
        ],
    )
}

#[test]
fn conics_keep_their_rho_in_range() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_conic(&mut s, f, 0.4);
    assert!(s.document().conflicts.is_empty());
    let id = |index| SketchId { op, index };
    edit(
        &mut s,
        f,
        vec![
            SketchEdit::SetRho {
                entity: id(3),
                rho: 0.7,
            },
            SketchEdit::SetValue {
                constraint: id(4),
                value: 0.7,
            },
        ],
    );
    assert!(s.document().conflicts.is_empty());
    let d = sketch(&s, f);
    assert!(matches!(
        d.sketch.entity(id(3)).unwrap().geometry,
        Geometry::Conic { rho, .. } if rho == 0.7
    ));
    assert_eq!(
        d.sketch.constraint(id(4)).unwrap().constraint.value(),
        Some(0.7)
    );

    let line_op = draw_line(&mut s, f);
    let before = sketch(&s, f);
    for e in [
        SketchEdit::SetRho {
            entity: id(3),
            rho: 1.0,
        },
        SketchEdit::SetRho {
            entity: id(3),
            rho: 0.0,
        },
        SketchEdit::SetRho {
            entity: id(3),
            rho: f64::NAN,
        },
        SketchEdit::SetRho {
            entity: SketchId {
                op: line_op,
                index: 2,
            },
            rho: 0.5,
        },
        SketchEdit::SetValue {
            constraint: id(4),
            value: 1.5,
        },
        SketchEdit::AddEntity(Entity {
            id: SketchId::new_in_batch(0),
            geometry: Geometry::Conic {
                start: id(0),
                control: id(1),
                end: id(2),
                rho: -0.2,
            },
            construction: false,
        }),
        SketchEdit::AddEntity(Entity {
            id: SketchId::new_in_batch(0),
            geometry: Geometry::Conic {
                start: id(0),
                control: id(0),
                end: id(2),
                rho: 0.5,
            },
            construction: false,
        }),
        constraint(
            0,
            Constraint::Rho {
                curve: id(3),
                value: 0.0,
            },
        ),
    ] {
        let n_conflicts = s.document().conflicts.len();
        edit(&mut s, f, vec![e.clone()]);
        assert_eq!(sketch(&s, f), before, "{e:?}");
        assert_eq!(s.document().conflicts.len(), n_conflicts + 1, "{e:?}");
    }
}

#[test]
fn deleting_a_point_cascades_to_curves_and_constraints() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_line(&mut s, f);
    edit(
        &mut s,
        f,
        vec![SketchEdit::Delete {
            items: vec![SketchId { op, index: 1 }],
        }],
    );
    let d = sketch(&s, f);
    // The line used the point; its constraints referenced the line.
    assert_eq!(d.sketch.entities.len(), 1);
    assert!(d.sketch.constraints.is_empty());
    // Deleting again (as a concurrent replica might) is harmless.
    edit(
        &mut s,
        f,
        vec![SketchEdit::Delete {
            items: vec![SketchId { op, index: 1 }],
        }],
    );
    assert!(s.document().conflicts.is_empty());
}

#[test]
fn ellipses_are_checked_and_deleted_with_their_points() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let ellipse = |index: u32, g: Geometry<SketchId>| {
        SketchEdit::AddEntity(Entity {
            id: SketchId::new_in_batch(index),
            geometry: g,
            construction: false,
        })
    };
    let b = SketchId::new_in_batch;
    let op = edit(
        &mut s,
        f,
        vec![
            point(0, 0.0, 0.0),
            point(1, 2.0, 0.0),
            point(2, 0.0, 1.0),
            point(3, 2.0, 0.0),
            point(4, 0.0, 1.0),
            ellipse(
                5,
                Geometry::Ellipse {
                    center: b(0),
                    major: b(1),
                    minor: b(2),
                },
            ),
            ellipse(
                6,
                Geometry::EllipticalArc {
                    center: b(0),
                    major: b(1),
                    minor: b(2),
                    start: b(3),
                    end: b(4),
                },
            ),
        ],
    );
    assert!(s.document().conflicts.is_empty());
    assert_eq!(sketch(&s, f).sketch.entities.len(), 7);
    let id = |index| SketchId { op, index };
    // Repeated points are refused, each edit on its own.
    for g in [
        Geometry::Ellipse {
            center: id(0),
            major: id(0),
            minor: id(2),
        },
        Geometry::Ellipse {
            center: id(0),
            major: id(1),
            minor: id(1),
        },
        Geometry::EllipticalArc {
            center: id(0),
            major: id(1),
            minor: id(2),
            start: id(3),
            end: id(3),
        },
        Geometry::EllipticalArc {
            center: id(0),
            major: id(1),
            minor: id(2),
            start: id(0),
            end: id(4),
        },
    ] {
        let before = sketch(&s, f);
        let n = s.document().conflicts.len();
        edit(&mut s, f, vec![ellipse(0, g.clone())]);
        assert_eq!(sketch(&s, f), before, "{g:?}");
        assert_eq!(s.document().conflicts.len(), n + 1, "{g:?}");
    }
    // Deleting the minor end deletes both curves that use it.
    edit(&mut s, f, vec![SketchEdit::Delete { items: vec![id(2)] }]);
    let d = sketch(&s, f);
    assert_eq!(d.sketch.entities.len(), 4);
    assert!(d
        .sketch
        .entities
        .iter()
        .all(|e| matches!(e.geometry, Geometry::Point { .. })));
}

#[test]
fn a_batch_is_one_undo_step() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let empty = sketch(&s, f);
    draw_line(&mut s, f);
    let drawn = sketch(&s, f);
    assert!(s.undo());
    assert_eq!(sketch(&s, f), empty);
    assert!(s.redo());
    assert_eq!(sketch(&s, f), drawn);
}

#[test]
fn concurrent_edits_to_different_items_merge() {
    let mut a = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut a);
    let op = draw_line(&mut a, f);
    let mut b = Session::open(&a.save().unwrap(), ReplicaId(2)).unwrap();
    let len = SketchId { op, index: 4 };
    let p1 = SketchId { op, index: 1 };
    // A changes the length; B makes the line construction and moves a point.
    edit(
        &mut a,
        f,
        vec![SketchEdit::SetValue {
            constraint: len,
            value: 0.02,
        }],
    );
    edit(
        &mut b,
        f,
        vec![
            SketchEdit::SetConstruction {
                entity: SketchId { op, index: 2 },
                construction: true,
            },
            SketchEdit::MovePoints {
                points: vec![(p1, DVec2::new(0.02, 0.0))],
            },
        ],
    );
    let to_b = a.log().missing_for(&b.log().heads());
    let to_a = b.log().missing_for(&a.log().heads());
    b.integrate(to_b).unwrap();
    a.integrate(to_a).unwrap();
    assert_eq!(sketch(&a, f), sketch(&b, f));
    let d = sketch(&a, f);
    assert_eq!(
        d.sketch.constraint(len).unwrap().constraint.value(),
        Some(0.02)
    );
    assert!(
        d.sketch
            .entity(SketchId { op, index: 2 })
            .unwrap()
            .construction
    );
    assert_eq!(d.sketch.point(p1), Some(DVec2::new(0.02, 0.0)));
}

#[test]
fn a_constraint_on_a_concurrently_deleted_entity_fails_visibly() {
    let mut a = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut a);
    let op = draw_line(&mut a, f);
    let mut b = Session::open(&a.save().unwrap(), ReplicaId(2)).unwrap();
    let ln = SketchId { op, index: 2 };
    let deleted = edit(&mut a, f, vec![SketchEdit::Delete { items: vec![ln] }]);
    let added = edit(
        &mut b,
        f,
        vec![constraint(0, Constraint::Vertical { line: ln })],
    );
    let to_b = a.log().missing_for(&b.log().heads());
    let to_a = b.log().missing_for(&a.log().heads());
    b.integrate(to_b).unwrap();
    a.integrate(to_a).unwrap();
    assert_eq!(sketch(&a, f), sketch(&b, f));
    // Whichever order the ops replay in, the line is gone and the
    // constraint is either rejected (recorded, if it replays after the
    // delete) or removed with the line; it never dangles.
    let d = sketch(&a, f);
    assert!(d.sketch.entity(ln).is_none());
    assert!(d.sketch.constraints.iter().all(|c| c
        .constraint
        .referenced()
        .iter()
        .all(|e| d.sketch.entity(*e).is_some())));
    let rejected = a.document().conflicts.iter().any(|c| c.op == added);
    let present = d
        .sketch
        .constraint(SketchId {
            op: added,
            index: 0,
        })
        .is_some();
    assert!(!present);
    assert_eq!(rejected, added > deleted);
}

#[test]
fn sketch_ops_round_trip_through_save_and_open() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    draw_line(&mut s, f);
    let reopened = Session::open(&s.save().unwrap(), ReplicaId(3)).unwrap();
    assert_eq!(sketch(&reopened, f), sketch(&s, f));
}

#[test]
fn dimension_labels_are_placed_moved_and_returned_to_automatic() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_line(&mut s, f);
    let len = SketchId { op, index: 4 };
    let place = |offset| {
        vec![SketchEdit::PlaceLabel {
            dimension: len,
            offset,
        }]
    };
    assert_eq!(sketch(&s, f).label(len), None);
    edit(&mut s, f, place(Some(DVec2::new(0.0, 0.004))));
    assert_eq!(sketch(&s, f).label(len), Some(DVec2::new(0.0, 0.004)));
    edit(&mut s, f, place(Some(DVec2::new(0.001, -0.003))));
    let d = sketch(&s, f);
    assert_eq!(d.label(len), Some(DVec2::new(0.001, -0.003)));
    assert_eq!(d.labels.len(), 1);
    // A second dimension's label keeps the list sorted by id.
    let n = SketchId::new_in_batch;
    let op2 = edit(
        &mut s,
        f,
        vec![
            point(0, 0.0, 0.01),
            constraint(
                1,
                Constraint::Distance {
                    a: SketchId { op, index: 0 },
                    b: n(0),
                    kind: opendatum_core_model::sketch::DistanceKind::Aligned,
                    value: 0.01,
                },
            ),
            SketchEdit::PlaceLabel {
                dimension: n(1),
                offset: Some(DVec2::new(-0.002, 0.0)),
            },
            SketchEdit::PlaceLabel {
                dimension: len,
                offset: Some(DVec2::new(0.0, 0.005)),
            },
        ],
    );
    let d = sketch(&s, f);
    let dist = SketchId { op: op2, index: 1 };
    assert_eq!(
        d.labels,
        vec![
            (len, DVec2::new(0.0, 0.005)),
            (dist, DVec2::new(-0.002, 0.0))
        ]
    );
    // Back to automatic placement; again is harmless.
    edit(&mut s, f, place(None));
    edit(&mut s, f, place(None));
    let d = sketch(&s, f);
    assert_eq!(d.label(len), None);
    assert_eq!(d.label(dist), Some(DVec2::new(-0.002, 0.0)));
    assert!(s.document().conflicts.is_empty());
    // Placing a label is one undo step.
    assert!(s.undo());
    assert!(s.undo());
    assert_eq!(sketch(&s, f).label(len), Some(DVec2::new(0.0, 0.005)));
}

#[test]
fn labels_only_go_on_dimensions_and_finite_offsets() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_line(&mut s, f);
    let id = |index| SketchId { op, index };
    for (dimension, offset) in [
        // The horizontal constraint, the line, and an id never used.
        (id(3), DVec2::ZERO),
        (id(2), DVec2::ZERO),
        (id(77), DVec2::ZERO),
        (id(4), DVec2::new(f64::NAN, 0.0)),
        (id(4), DVec2::new(0.0, f64::INFINITY)),
    ] {
        let before = sketch(&s, f);
        let n_conflicts = s.document().conflicts.len();
        edit(
            &mut s,
            f,
            vec![SketchEdit::PlaceLabel {
                dimension,
                offset: Some(offset),
            }],
        );
        assert_eq!(sketch(&s, f), before, "{dimension:?} {offset}");
        assert_eq!(s.document().conflicts.len(), n_conflicts + 1);
    }
}

#[test]
fn deleting_a_dimension_or_what_it_measures_drops_its_label() {
    let mut s = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut s);
    let op = draw_line(&mut s, f);
    let len = SketchId { op, index: 4 };
    let place = SketchEdit::PlaceLabel {
        dimension: len,
        offset: Some(DVec2::new(0.0, 0.004)),
    };
    edit(&mut s, f, vec![place.clone()]);
    edit(&mut s, f, vec![SketchEdit::Delete { items: vec![len] }]);
    assert!(sketch(&s, f).labels.is_empty());
    assert!(s.undo());
    assert_eq!(sketch(&s, f).labels.len(), 1);
    // Deleting an end point takes the line, its length, and the label.
    edit(
        &mut s,
        f,
        vec![SketchEdit::Delete {
            items: vec![SketchId { op, index: 0 }],
        }],
    );
    let d = sketch(&s, f);
    assert!(d.sketch.constraints.is_empty());
    assert!(d.labels.is_empty());
}

#[test]
fn concurrent_label_moves_converge_and_labels_survive_save() {
    let mut a = Session::new_document(ReplicaId(1));
    let f = new_sketch(&mut a);
    let op = draw_line(&mut a, f);
    let len = SketchId { op, index: 4 };
    let mut b = Session::open(&a.save().unwrap(), ReplicaId(2)).unwrap();
    let place = |s: &mut Session, y| {
        edit(
            s,
            f,
            vec![SketchEdit::PlaceLabel {
                dimension: len,
                offset: Some(DVec2::new(0.0, y)),
            }],
        )
    };
    let in_a = place(&mut a, 0.003);
    let in_b = place(&mut b, -0.003);
    let to_b = a.log().missing_for(&b.log().heads());
    let to_a = b.log().missing_for(&a.log().heads());
    b.integrate(to_b).unwrap();
    a.integrate(to_a).unwrap();
    assert_eq!(sketch(&a, f), sketch(&b, f));
    // The later op in replay order wins.
    let last = if in_a > in_b { 0.003 } else { -0.003 };
    assert_eq!(sketch(&a, f).label(len), Some(DVec2::new(0.0, last)));
    let reopened = Session::open(&a.save().unwrap(), ReplicaId(3)).unwrap();
    assert_eq!(sketch(&reopened, f), sketch(&a, f));
}

#[test]
fn a_sketch_saved_without_labels_reads_with_none() {
    let def = SketchDef::new(PlaneRef::Front, SketchData::default());
    let mut json: serde_json::Value = serde_json::to_value(&def).unwrap();
    assert!(json.as_object_mut().unwrap().remove("labels").is_some());
    let back: SketchDef = serde_json::from_value(json).unwrap();
    assert_eq!(back, def);
}
