//! Behavioural and property tests for the operation log.

use openshape_oplog::{Integrated, OpKind, OpLog, OpLogError, ReplicaId};
use proptest::prelude::*;

fn payloads(log: &OpLog<u32>) -> Vec<u32> {
    log.effective().into_iter().map(|(_, p)| *p).collect()
}

#[test]
fn append_undo_redo() {
    let mut log = OpLog::new(ReplicaId(1));
    log.append(10);
    log.append(20);
    assert_eq!(payloads(&log), vec![10, 20]);
    log.undo();
    assert_eq!(payloads(&log), vec![10]);
    log.undo();
    assert_eq!(payloads(&log), Vec::<u32>::new());
    assert!(log.undo().is_none());
    log.redo();
    assert_eq!(payloads(&log), vec![10]);
    log.redo();
    assert_eq!(payloads(&log), vec![10, 20]);
    assert!(log.redo().is_none());
    // Undo after redo works (undo of the redo op).
    log.undo();
    assert_eq!(payloads(&log), vec![10]);
    // A new action clears redo.
    log.append(30);
    assert!(!log.can_redo());
    assert_eq!(payloads(&log), vec![10, 30]);
    // Nothing is ever removed from the log.
    assert_eq!(log.len(), 8);
}

#[test]
fn history_restore_with_effective_at() {
    let mut log = OpLog::new(ReplicaId(1));
    log.append(1);
    let v1 = log.heads();
    log.append(2);
    log.undo();
    log.append(3);
    assert_eq!(payloads(&log), vec![1, 3]);
    let at: Vec<u32> = log.effective_at(&v1).into_iter().map(|(_, p)| *p).collect();
    assert_eq!(at, vec![1]);
}

#[test]
fn encode_decode_roundtrip_preserves_state_and_undo() {
    let mut log = OpLog::new(ReplicaId(7));
    log.append(1);
    log.append(2);
    log.undo();
    let bytes = log.encode().unwrap();
    let mut back: OpLog<u32> = OpLog::decode(&bytes, ReplicaId(7)).unwrap();
    assert_eq!(payloads(&back), vec![1]);
    assert_eq!(back.heads(), log.heads());
    assert!(back.can_redo());
    back.redo();
    assert_eq!(payloads(&back), vec![1, 2]);
}

#[test]
fn decode_rejects_garbage() {
    assert!(matches!(
        OpLog::<u32>::decode(b"not cbor at all", ReplicaId(1)),
        Err(OpLogError::Decode(_))
    ));
}

#[test]
fn out_of_order_delivery_is_buffered() {
    let mut a = OpLog::new(ReplicaId(1));
    a.append(1);
    a.append(2);
    a.append(3);
    let ops: Vec<_> = a.ops().cloned().collect();
    let mut b: OpLog<u32> = OpLog::new(ReplicaId(2));
    assert_eq!(b.integrate(ops[2].clone()), Ok(Integrated::Pending));
    assert_eq!(b.integrate(ops[1].clone()), Ok(Integrated::Pending));
    assert_eq!(b.integrate(ops[0].clone()), Ok(Integrated::Added));
    assert_eq!(b.pending_len(), 0);
    assert_eq!(payloads(&b), vec![1, 2, 3]);
    assert_eq!(b.integrate(ops[0].clone()), Ok(Integrated::Duplicate));
}

#[test]
fn rejects_causality_violation_and_bad_undo() {
    use openshape_oplog::{Op, OpId};
    let mut a = OpLog::new(ReplicaId(1));
    let x = a.append(1);
    let mut b: OpLog<u32> = OpLog::new(ReplicaId(2));
    b.integrate(a.get(x).unwrap().clone()).unwrap();
    let bad_id = OpId {
        lamport: 1,
        replica: ReplicaId(9),
    };
    let bad = Op {
        id: bad_id,
        deps: vec![x],
        kind: OpKind::Apply(5),
    };
    assert_eq!(
        b.integrate(bad),
        Err(OpLogError::CausalityViolation(bad_id))
    );
    let undo = Op {
        id: OpId {
            lamport: 3,
            replica: ReplicaId(3),
        },
        deps: vec![],
        kind: OpKind::Undo { target: x },
    };
    assert!(matches!(
        b.integrate(undo),
        Err(OpLogError::BadUndoTarget(_))
    ));
    // Same id, different content.
    let clash = Op {
        id: x,
        deps: vec![],
        kind: OpKind::Apply(99),
    };
    assert_eq!(b.integrate(clash), Err(OpLogError::ConflictingDuplicate(x)));
}

/// A user action in the randomized tests.
#[derive(Clone, Debug)]
enum Action {
    Append(u32),
    Undo,
    Redo,
}

fn action() -> impl Strategy<Value = Action> {
    prop_oneof![
        3 => (0u32..1000).prop_map(Action::Append),
        2 => Just(Action::Undo),
        1 => Just(Action::Redo),
    ]
}

proptest! {
    /// Undo/redo on the log behaves exactly like a classic two-stack editor.
    #[test]
    fn undo_redo_matches_reference_model(actions in prop::collection::vec(action(), 0..60)) {
        let mut log = OpLog::new(ReplicaId(1));
        let mut done: Vec<u32> = Vec::new();
        let mut undone: Vec<u32> = Vec::new();
        for a in &actions {
            match a {
                Action::Append(v) => { log.append(*v); done.push(*v); undone.clear(); }
                Action::Undo => { log.undo(); if let Some(v) = done.pop() { undone.push(v); } }
                Action::Redo => { log.redo(); if let Some(v) = undone.pop() { done.push(v); } }
            }
            // The effective set equals the reference, in the same order.
            prop_assert_eq!(payloads(&log), done.clone());
        }
        // Rebuilding the stacks from the log gives the same behaviour.
        let mut rebuilt = OpLog::decode(&log.encode().unwrap(), ReplicaId(1)).unwrap();
        prop_assert_eq!(rebuilt.can_undo(), !done.is_empty());
        prop_assert_eq!(rebuilt.can_redo(), !undone.is_empty());
        while rebuilt.undo().is_some() {}
        prop_assert!(payloads(&rebuilt).is_empty());
    }

    /// Replicas that exchange all ops in any order converge to the same state.
    #[test]
    fn replicas_converge(
        a_actions in prop::collection::vec(action(), 0..25),
        b_actions in prop::collection::vec(action(), 0..25),
        sync_points in prop::collection::vec(any::<bool>(), 50),
        shuffle_seed in any::<u64>(),
    ) {
        let mut a = OpLog::new(ReplicaId(1));
        let mut b = OpLog::new(ReplicaId(2));
        let n = a_actions.len().max(b_actions.len());
        let apply = |log: &mut OpLog<u32>, act: &Action| match act {
            Action::Append(v) => { log.append(*v); }
            Action::Undo => { log.undo(); }
            Action::Redo => { log.redo(); }
        };
        for i in 0..n {
            if let Some(x) = a_actions.get(i) { apply(&mut a, x); }
            if let Some(x) = b_actions.get(i) { apply(&mut b, x); }
            if sync_points[i % sync_points.len()] {
                for op in a.missing_for(&b.heads()) { b.integrate(op).unwrap(); }
                for op in b.missing_for(&a.heads()) { a.integrate(op).unwrap(); }
            }
        }
        // Final exchange in a scrambled order to exercise buffering.
        let mut to_b = a.missing_for(&b.heads());
        let mut to_a = b.missing_for(&a.heads());
        let mut s = shuffle_seed | 1;
        let mut scramble = |v: &mut Vec<_>| {
            for i in (1..v.len()).rev() {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                v.swap(i, (s % (i as u64 + 1)) as usize);
            }
        };
        scramble(&mut to_b);
        scramble(&mut to_a);
        for op in to_b { b.integrate(op).unwrap(); }
        for op in to_a { a.integrate(op).unwrap(); }
        prop_assert_eq!(a.pending_len(), 0);
        prop_assert_eq!(b.pending_len(), 0);
        prop_assert_eq!(a.heads(), b.heads());
        prop_assert_eq!(payloads(&a), payloads(&b));
    }
}
