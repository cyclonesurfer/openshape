//! The rollback bar and element deletion in the op log (ADR 0013).

use opendatum_core_model::{
    DVec3, ElementId, ElementKind, Feature, FeatureId, InsertAt, ModelOp, OpId, PartStudio,
    ReplicaId, Session,
};
use proptest::prelude::*;

fn box_feature(s: f64) -> Feature {
    Feature::Box {
        corner: DVec3::ZERO,
        size: DVec3::splat(s),
    }
}

fn studio(s: &Session) -> &PartStudio {
    &s.document().elements[0].studio
}

fn names(s: &Session) -> Vec<String> {
    studio(s).features.iter().map(|f| f.name.clone()).collect()
}

fn active(s: &Session) -> Vec<String> {
    studio(s)
        .active_features()
        .iter()
        .map(|f| f.name.clone())
        .collect()
}

fn insert(s: &mut Session, el: ElementId, at: InsertAt, name: &str) -> FeatureId {
    FeatureId(s.apply(ModelOp::InsertFeature {
        element: el,
        at,
        name: name.into(),
        feature: box_feature(1.0),
    }))
}

/// A document with features A B C D and the bar where `bar` says.
fn abcd(bar: Option<usize>) -> (Session, ElementId, Vec<FeatureId>) {
    let mut s = Session::new_document(ReplicaId(1));
    let el = s.first_studio().unwrap();
    let ids: Vec<FeatureId> = ["A", "B", "C", "D"]
        .iter()
        .map(|n| insert(&mut s, el, InsertAt::End, n))
        .collect();
    if let Some(i) = bar {
        let at = if i == 0 {
            InsertAt::Start
        } else {
            InsertAt::After(ids[i - 1])
        };
        s.apply(ModelOp::SetRollback { element: el, at });
    }
    (s, el, ids)
}

#[test]
fn features_after_the_bar_are_rolled_back() {
    let (s, _, ids) = abcd(Some(2));
    assert_eq!(names(&s), ["A", "B", "C", "D"]);
    assert_eq!(active(&s), ["A", "B"]);
    assert_eq!(studio(&s).rollback_index(), 2);
    assert!(!studio(&s).is_rolled_back(ids[1]));
    assert!(studio(&s).is_rolled_back(ids[2]));
    assert!(studio(&s).is_rolled_back(ids[3]));
    assert!(s.document().conflicts.is_empty());
}

#[test]
fn a_new_document_has_the_bar_at_the_end() {
    let (s, _, _) = abcd(None);
    assert_eq!(studio(&s).rollback, InsertAt::End);
    assert_eq!(active(&s), ["A", "B", "C", "D"]);
}

#[test]
fn inserting_at_the_bar_puts_the_feature_above_it() {
    let (mut s, el, ids) = abcd(Some(2));
    let e = insert(&mut s, el, InsertAt::After(ids[1]), "E");
    assert_eq!(names(&s), ["A", "B", "E", "C", "D"]);
    assert_eq!(active(&s), ["A", "B", "E"]);
    assert_eq!(studio(&s).rollback, InsertAt::After(e));
    // From the start too.
    s.apply(ModelOp::SetRollback {
        element: el,
        at: InsertAt::Start,
    });
    assert!(active(&s).is_empty());
    insert(&mut s, el, InsertAt::Start, "F");
    assert_eq!(active(&s), ["F"]);
}

#[test]
fn inserting_before_or_after_the_bar_keeps_it_in_place() {
    let (mut s, el, ids) = abcd(Some(2));
    insert(&mut s, el, InsertAt::After(ids[0]), "E");
    assert_eq!(active(&s), ["A", "E", "B"]);
    insert(&mut s, el, InsertAt::End, "F");
    assert_eq!(names(&s), ["A", "E", "B", "C", "D", "F"]);
    assert_eq!(active(&s), ["A", "E", "B"]);
}

#[test]
fn appending_with_the_bar_at_the_end_stays_active() {
    let (mut s, el, _) = abcd(None);
    insert(&mut s, el, InsertAt::End, "E");
    assert_eq!(studio(&s).rollback, InsertAt::End);
    assert_eq!(active(&s).len(), 5);
}

#[test]
fn the_bar_after_the_last_feature_is_the_end() {
    let (mut s, el, ids) = abcd(None);
    s.apply(ModelOp::SetRollback {
        element: el,
        at: InsertAt::After(ids[3]),
    });
    assert_eq!(studio(&s).rollback, InsertAt::End);
    // So a later append is active.
    insert(&mut s, el, InsertAt::End, "E");
    assert_eq!(active(&s).len(), 5);
}

#[test]
fn deleting_the_anchor_keeps_the_same_features_active() {
    let (mut s, el, ids) = abcd(Some(2));
    s.apply(ModelOp::DeleteFeature {
        element: el,
        feature: ids[1],
    });
    assert_eq!(active(&s), ["A"]);
    assert_eq!(studio(&s).rollback, InsertAt::After(ids[0]));
    s.apply(ModelOp::DeleteFeature {
        element: el,
        feature: ids[0],
    });
    assert!(active(&s).is_empty());
    assert_eq!(studio(&s).rollback, InsertAt::Start);
    // Deleting everything after the bar leaves it at the end.
    s.apply(ModelOp::DeleteFeature {
        element: el,
        feature: ids[2],
    });
    s.apply(ModelOp::DeleteFeature {
        element: el,
        feature: ids[3],
    });
    assert_eq!(studio(&s).rollback, InsertAt::End);
    assert!(s.document().conflicts.is_empty());
}

#[test]
fn moving_features_across_the_bar() {
    // A B | C D: moving the anchor below the bar rolls it back.
    let (mut s, el, ids) = abcd(Some(2));
    s.apply(ModelOp::MoveFeature {
        element: el,
        feature: ids[1],
        to: InsertAt::After(ids[3]),
    });
    assert_eq!(names(&s), ["A", "C", "D", "B"]);
    assert_eq!(active(&s), ["A"]);
    // A | C D B: moving D to the bar puts it above the bar.
    s.apply(ModelOp::MoveFeature {
        element: el,
        feature: ids[3],
        to: InsertAt::After(ids[0]),
    });
    assert_eq!(names(&s), ["A", "D", "C", "B"]);
    assert_eq!(active(&s), ["A", "D"]);
    // A D | C B: moving C above A activates it, the bar stays after D.
    s.apply(ModelOp::MoveFeature {
        element: el,
        feature: ids[2],
        to: InsertAt::Start,
    });
    assert_eq!(names(&s), ["C", "A", "D", "B"]);
    assert_eq!(active(&s), ["C", "A", "D"]);
    assert!(s.document().conflicts.is_empty());
}

#[test]
fn a_bar_after_a_missing_feature_is_a_conflict() {
    let (mut s, el, _) = abcd(Some(2));
    let id = s.apply(ModelOp::SetRollback {
        element: el,
        at: InsertAt::After(FeatureId(OpId {
            lamport: 999,
            replica: ReplicaId(9),
        })),
    });
    assert_eq!(s.document().conflicts.len(), 1);
    assert_eq!(s.document().conflicts[0].op, id);
    assert_eq!(active(&s), ["A", "B"], "the bar did not move");
}

#[test]
fn moving_the_bar_is_undoable_and_saved() {
    let (mut s, el, ids) = abcd(None);
    s.apply(ModelOp::SetRollback {
        element: el,
        at: InsertAt::After(ids[0]),
    });
    assert_eq!(active(&s), ["A"]);
    let bytes = s.save().unwrap();
    let back = Session::open(&bytes, ReplicaId(1)).unwrap();
    assert_eq!(back.document(), s.document());
    assert!(s.undo());
    assert_eq!(active(&s).len(), 4);
    assert!(s.redo());
    assert_eq!(active(&s), ["A"]);
}

#[test]
fn undo_depth_counts_undoable_actions() {
    let (mut s, el, _) = abcd(None);
    assert_eq!(s.undo_depth(), 4);
    insert(&mut s, el, InsertAt::End, "E");
    assert_eq!(s.undo_depth(), 5);
    while s.undo_depth() > 3 {
        assert!(s.undo());
    }
    assert_eq!(names(&s), ["A", "B", "C"]);
}

#[test]
fn elements_can_be_deleted_but_not_the_last() {
    let mut s = Session::new_document(ReplicaId(1));
    let first = s.first_studio().unwrap();
    let second = ElementId(s.apply(ModelOp::CreateElement {
        name: "Part Studio 2".into(),
        kind: ElementKind::PartStudio,
    }));
    insert(&mut s, second, InsertAt::End, "A");
    s.apply(ModelOp::DeleteElement { element: first });
    assert_eq!(s.document().elements.len(), 1);
    assert_eq!(s.document().elements[0].id, second);
    let refused = s.apply(ModelOp::DeleteElement { element: second });
    assert_eq!(s.document().elements.len(), 1);
    assert_eq!(s.document().conflicts.len(), 1);
    assert_eq!(s.document().conflicts[0].op, refused);
    // Undo brings the first element back with its contents.
    assert!(s.undo());
    assert!(s.undo());
    assert_eq!(s.document().elements.len(), 2);
    assert_eq!(s.document().elements[0].id, first);
}

#[test]
fn concurrent_deletes_of_the_two_last_elements_keep_one() {
    let mut a = Session::new_document(ReplicaId(1));
    let first = a.first_studio().unwrap();
    let second = ElementId(a.apply(ModelOp::CreateElement {
        name: "Part Studio 2".into(),
        kind: ElementKind::PartStudio,
    }));
    let mut b = Session::open(&a.save().unwrap(), ReplicaId(2)).unwrap();
    a.apply(ModelOp::DeleteElement { element: first });
    b.apply(ModelOp::DeleteElement { element: second });
    let to_b = a.log().missing_for(&b.log().heads());
    let to_a = b.log().missing_for(&a.log().heads());
    a.integrate(to_a).unwrap();
    b.integrate(to_b).unwrap();
    assert_eq!(a.document(), b.document());
    assert_eq!(a.document().elements.len(), 1);
    assert_eq!(a.document().conflicts.len(), 1);
}

/// An edit for the property tests; indices pick features modulo the count.
#[derive(Clone, Debug)]
enum Edit {
    Insert(Pos),
    Move(usize, Pos),
    Delete(usize),
    Bar(Pos),
    Suppress(usize),
    Undo,
    Redo,
}

#[derive(Clone, Copy, Debug)]
enum Pos {
    Start,
    End,
    After(usize),
    AtBar,
}

fn pos() -> impl Strategy<Value = Pos> {
    prop_oneof![
        Just(Pos::Start),
        Just(Pos::End),
        (0usize..8).prop_map(Pos::After),
        Just(Pos::AtBar),
    ]
}

fn edit() -> impl Strategy<Value = Edit> {
    prop_oneof![
        3 => pos().prop_map(Edit::Insert),
        2 => (0usize..8, pos()).prop_map(|(i, p)| Edit::Move(i, p)),
        1 => (0usize..8).prop_map(Edit::Delete),
        2 => pos().prop_map(Edit::Bar),
        1 => (0usize..8).prop_map(Edit::Suppress),
        1 => Just(Edit::Undo),
        1 => Just(Edit::Redo),
    ]
}

fn resolve(st: &PartStudio, p: Pos) -> InsertAt {
    match p {
        Pos::Start => InsertAt::Start,
        Pos::End => InsertAt::End,
        Pos::After(i) if !st.features.is_empty() => {
            InsertAt::After(st.features[i % st.features.len()].id)
        }
        Pos::After(_) => InsertAt::End,
        Pos::AtBar => st.rollback,
    }
}

/// The index a position resolves to in `ids`.
fn index_in(ids: &[FeatureId], at: InsertAt) -> usize {
    match at {
        InsertAt::Start => 0,
        InsertAt::End => ids.len(),
        InsertAt::After(f) => ids.iter().position(|x| *x == f).unwrap() + 1,
    }
}

/// The reference model: the feature order and the number of active
/// features, updated by the rules of ADR 0013 written as index arithmetic.
#[derive(Clone, Debug, PartialEq)]
struct Reference {
    ids: Vec<FeatureId>,
    bar: usize,
}

impl Reference {
    fn insert(&mut self, at: usize, id: FeatureId) {
        let end = self.bar == self.ids.len();
        self.ids.insert(at, id);
        if end || at <= self.bar {
            self.bar += 1;
        }
    }

    fn remove(&mut self, i: usize) -> FeatureId {
        if i < self.bar {
            self.bar -= 1;
        }
        self.ids.remove(i)
    }
}

fn check_normalized(st: &PartStudio) {
    match st.rollback {
        InsertAt::End => {}
        InsertAt::Start => assert!(!st.features.is_empty(), "Start on an empty list"),
        InsertAt::After(f) => {
            let i = st.index_of(f).expect("the anchor exists");
            assert!(i + 1 < st.features.len(), "After(last) instead of End");
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    /// The reducer follows the reference rules op by op, stays normalized,
    /// and incremental application equals replay (including undo and redo).
    #[test]
    fn rollback_rules_match_the_reference(edits in proptest::collection::vec(edit(), 1..40)) {
        let mut s = Session::new_document(ReplicaId(1));
        let el = s.first_studio().unwrap();
        let mut reference = Reference { ids: Vec::new(), bar: 0 };
        for e in edits {
            let st = studio(&s).clone();
            let n = st.features.len();
            let op = match e {
                Edit::Insert(p) => Some(ModelOp::InsertFeature {
                    element: el,
                    at: resolve(&st, p),
                    name: "F".into(),
                    feature: box_feature(1.0),
                }),
                Edit::Move(i, p) if n > 0 => Some(ModelOp::MoveFeature {
                    element: el,
                    feature: st.features[i % n].id,
                    to: resolve(&st, p),
                }),
                Edit::Delete(i) if n > 0 => Some(ModelOp::DeleteFeature {
                    element: el,
                    feature: st.features[i % n].id,
                }),
                Edit::Bar(p) => Some(ModelOp::SetRollback { element: el, at: resolve(&st, p) }),
                Edit::Suppress(i) if n > 0 => Some(ModelOp::SetSuppressed {
                    element: el,
                    feature: st.features[i % n].id,
                    suppressed: !st.features[i % n].suppressed,
                }),
                Edit::Undo => {
                    s.undo();
                    None
                }
                Edit::Redo => {
                    s.redo();
                    None
                }
                _ => None,
            };
            match op {
                Some(op) => {
                    let id = s.apply(op.clone());
                    prop_assert!(s.document().conflicts.iter().all(|c| c.op != id), "{op:?} conflicted");
                    match op {
                        ModelOp::InsertFeature { at, .. } => {
                            let i = index_in(&reference.ids, at);
                            reference.insert(i, FeatureId(id));
                        }
                        ModelOp::MoveFeature { feature, to, .. } => {
                            if to != InsertAt::After(feature) {
                                let from = reference.ids.iter().position(|x| *x == feature).unwrap();
                                reference.remove(from);
                                let i = index_in(&reference.ids, to);
                                reference.insert(i, feature);
                            }
                        }
                        ModelOp::DeleteFeature { feature, .. } => {
                            let i = reference.ids.iter().position(|x| *x == feature).unwrap();
                            reference.remove(i);
                        }
                        ModelOp::SetRollback { at, .. } => {
                            reference.bar = index_in(&reference.ids, at);
                        }
                        _ => {}
                    }
                }
                None => {
                    // Undo/redo: re-derive the reference from the document.
                    let st = studio(&s);
                    reference = Reference {
                        ids: st.features.iter().map(|f| f.id).collect(),
                        bar: st.rollback_index(),
                    };
                }
            }
            let st = studio(&s);
            check_normalized(st);
            let ids: Vec<FeatureId> = st.features.iter().map(|f| f.id).collect();
            prop_assert_eq!(&ids, &reference.ids);
            prop_assert_eq!(st.rollback_index(), reference.bar);
            let replayed = opendatum_core_model::Document::replay(s.log());
            prop_assert_eq!(&replayed, s.document());
        }
    }
}
