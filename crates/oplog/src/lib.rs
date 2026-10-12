//! The OpenDatum operation log.
//!
//! A document is an append-only log of operations. Everything else (the
//! current model, undo/redo, history, versions, branches, merge and sync) is
//! derived from it. See `docs/adr/0004-oplog-and-sync.md` for the design.
//!
//! Key properties:
//!
//! * **Stable IDs.** Every operation is identified by an [`OpId`], a Lamport
//!   timestamp paired with the [`ReplicaId`] that created it. IDs never change
//!   once assigned, and two replicas never mint the same ID.
//! * **Causal ordering.** Every operation records the heads of the log it was
//!   created against (`deps`). An operation is only integrated once all of its
//!   dependencies are present, so the log is always causally closed.
//! * **Deterministic replay.** [`OpId`] ordering (Lamport time, then replica)
//!   is a total order that extends the causal order. Replaying the effective
//!   operations in that order gives every replica the same state once they
//!   hold the same set of operations, whatever order they received them in.
//! * **Undo is an operation.** Undo appends an [`OpKind::Undo`] op that
//!   cancels its target; redo is an undo of that undo. Nothing is ever
//!   removed, so undo composes with history, versions and sync.
//!
//! The log is generic over the payload type `P`; the document model supplies
//! the payload and a deterministic reducer.

use std::collections::{BTreeMap, BTreeSet};

use serde::{de::DeserializeOwned, Deserialize, Serialize};

/// Identifies one writer of a log (one editing session on one device).
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Default,
)]
pub struct ReplicaId(pub u64);

impl ReplicaId {
    /// Creates a fresh, practically unique replica id from the system clock,
    /// the process id and a per-process counter (no external RNG dependency).
    pub fn generate() -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        use std::time::{SystemTime, UNIX_EPOCH};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(0);
        let pid = u64::from(std::process::id());
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        // SplitMix64 finaliser over the mixed inputs.
        let mut z = nanos ^ pid.rotate_left(32) ^ n.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ReplicaId(z ^ (z >> 31))
    }
}

/// Globally unique, totally ordered identifier of an operation.
///
/// Ordering is by Lamport time first, then replica. This total order is
/// consistent with causality: if `a` happened before `b`, then `a < b`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct OpId {
    /// Lamport timestamp.
    pub lamport: u64,
    /// Replica that created the operation.
    pub replica: ReplicaId,
}

impl OpId {
    /// A reserved id no real operation has (Lamport times start at 1). The
    /// document model uses it to mean "the operation being applied" inside
    /// an op's payload.
    pub const NULL: OpId = OpId {
        lamport: 0,
        replica: ReplicaId(0),
    };
}

impl std::fmt::Display for OpId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}@{:016x}", self.lamport, self.replica.0)
    }
}

/// What an operation does.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum OpKind<P> {
    /// Applies a model payload.
    Apply(P),
    /// Cancels the target operation (which may itself be an undo, giving redo).
    Undo {
        /// The operation being cancelled.
        target: OpId,
    },
}

/// One entry in the log.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Op<P> {
    /// Unique id.
    pub id: OpId,
    /// Causal dependencies: the heads of the log this op was created against.
    pub deps: Vec<OpId>,
    /// The operation itself.
    pub kind: OpKind<P>,
}

/// Errors produced by the log.
#[derive(Clone, Debug, thiserror::Error, PartialEq)]
pub enum OpLogError {
    /// An operation's Lamport time is not greater than one of its dependencies.
    #[error("operation {0} violates causal ordering")]
    CausalityViolation(OpId),
    /// An undo targets an operation that is not in its causal past.
    #[error("undo {0} targets an operation outside its causal past")]
    BadUndoTarget(OpId),
    /// An operation uses the reserved [`OpId::NULL`] id.
    #[error("operation uses the reserved null id")]
    NullId,
    /// An operation with this id already exists with different content.
    #[error("conflicting duplicate operation {0}")]
    ConflictingDuplicate(OpId),
    /// The serialized log could not be decoded.
    #[error("decode error: {0}")]
    Decode(String),
    /// The serialized log could not be encoded.
    #[error("encode error: {0}")]
    Encode(String),
    /// The serialized log has an unknown format or version.
    #[error("unsupported log format: {0}")]
    UnsupportedFormat(String),
}

/// Result of integrating a remote operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Integrated {
    /// The operation (and possibly buffered ones it unblocked) was added.
    Added,
    /// The operation was already present.
    Duplicate,
    /// Some dependencies are missing; the op is buffered until they arrive.
    Pending,
}

/// The append-only operation log.
#[derive(Clone, Debug)]
pub struct OpLog<P> {
    replica: ReplicaId,
    clock: u64,
    ops: BTreeMap<OpId, Op<P>>,
    heads: BTreeSet<OpId>,
    pending: Vec<Op<P>>,
    rejected: Vec<(OpId, OpLogError)>,
    undo_stack: Vec<OpId>,
    redo_stack: Vec<OpId>,
}

/// Magic bytes and version of the serialized log container.
const MAGIC: &str = "OpenDatum-oplog";
/// The magic written before the product was renamed (ADR 0017); such
/// documents are read as they are, since the container is unchanged.
const LEGACY_MAGIC: &str = "OpenShape-oplog";
const FORMAT_VERSION: u32 = 1;

#[derive(Serialize, Deserialize)]
struct Container<P> {
    magic: String,
    version: u32,
    ops: Vec<Op<P>>,
}

impl<P: Clone + PartialEq> OpLog<P> {
    /// Creates an empty log written by `replica`.
    pub fn new(replica: ReplicaId) -> Self {
        Self {
            replica,
            clock: 0,
            ops: BTreeMap::new(),
            heads: BTreeSet::new(),
            pending: Vec::new(),
            rejected: Vec::new(),
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
        }
    }

    /// The replica that writes local operations into this log.
    pub fn replica(&self) -> ReplicaId {
        self.replica
    }

    /// Number of integrated operations.
    pub fn len(&self) -> usize {
        self.ops.len()
    }

    /// True if the log holds no operations.
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// The current heads (operations no other operation depends on).
    pub fn heads(&self) -> Vec<OpId> {
        self.heads.iter().copied().collect()
    }

    /// Looks up an operation.
    pub fn get(&self, id: OpId) -> Option<&Op<P>> {
        self.ops.get(&id)
    }

    /// All operations in total (replay) order.
    pub fn ops(&self) -> impl Iterator<Item = &Op<P>> {
        self.ops.values()
    }

    /// Operations received but not yet integrated (missing dependencies).
    pub fn pending_len(&self) -> usize {
        self.pending.len()
    }

    /// Buffered operations that turned out to be invalid once their
    /// dependencies arrived (e.g. an undo whose target is not in its causal
    /// past). They are dropped and recorded here; the integration that
    /// unblocked them still succeeds.
    pub fn rejected(&self) -> &[(OpId, OpLogError)] {
        &self.rejected
    }

    fn next_id(&mut self) -> OpId {
        self.clock += 1;
        OpId {
            lamport: self.clock,
            replica: self.replica,
        }
    }

    fn insert_local(&mut self, kind: OpKind<P>) -> OpId {
        let id = self.next_id();
        let deps = self.heads();
        self.heads.clear();
        self.heads.insert(id);
        self.ops.insert(id, Op { id, deps, kind });
        id
    }

    /// Appends a local user action. Clears the redo stack.
    pub fn append(&mut self, payload: P) -> OpId {
        let id = self.insert_local(OpKind::Apply(payload));
        self.undo_stack.push(id);
        self.redo_stack.clear();
        id
    }

    /// True if this replica has an action to undo.
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// The number of this replica's actions that can be undone.
    pub fn undo_depth(&self) -> usize {
        self.undo_stack.len()
    }

    /// True if this replica has an undone action to redo.
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    /// Undoes this replica's most recent action by appending an undo op.
    /// Returns the new op id, or `None` if there is nothing to undo.
    pub fn undo(&mut self) -> Option<OpId> {
        let target = self.undo_stack.pop()?;
        let id = self.insert_local(OpKind::Undo { target });
        self.redo_stack.push(id);
        Some(id)
    }

    /// Redoes the most recently undone action by undoing its undo op.
    /// Returns the new op id, or `None` if there is nothing to redo.
    pub fn redo(&mut self) -> Option<OpId> {
        let target = self.redo_stack.pop()?;
        let id = self.insert_local(OpKind::Undo { target });
        self.undo_stack.push(id);
        Some(id)
    }

    /// Rebuilds this replica's undo and redo stacks from its own operations.
    /// Used after loading a log that this replica wrote earlier.
    pub fn rebuild_undo_stacks(&mut self) {
        self.undo_stack.clear();
        self.redo_stack.clear();
        for op in self.ops.values().filter(|o| o.id.replica == self.replica) {
            match &op.kind {
                OpKind::Apply(_) => {
                    self.undo_stack.push(op.id);
                    self.redo_stack.clear();
                }
                OpKind::Undo { target } => {
                    if self.undo_stack.last() == Some(target) {
                        self.undo_stack.pop();
                        self.redo_stack.push(op.id);
                    } else if self.redo_stack.last() == Some(target) {
                        self.redo_stack.pop();
                        self.undo_stack.push(op.id);
                    }
                }
            }
        }
    }

    /// The causal closure of `heads`: every op reachable through `deps`.
    /// Unknown ids are ignored.
    pub fn causal_closure(&self, heads: &[OpId]) -> BTreeSet<OpId> {
        let mut seen = BTreeSet::new();
        let mut stack: Vec<OpId> = heads
            .iter()
            .copied()
            .filter(|h| self.ops.contains_key(h))
            .collect();
        while let Some(id) = stack.pop() {
            if seen.insert(id) {
                if let Some(op) = self.ops.get(&id) {
                    stack.extend(op.deps.iter().copied());
                }
            }
        }
        seen
    }

    /// Computes which ops of `subset` are active (not cancelled by an active
    /// undo). Undo ops always sort after their target, so a single descending
    /// pass decides every op's state before its target is reached.
    fn active_within<'a>(
        &'a self,
        subset: impl DoubleEndedIterator<Item = &'a Op<P>>,
    ) -> BTreeSet<OpId> {
        let mut cancelled = BTreeSet::new();
        let mut active = BTreeSet::new();
        for op in subset.rev() {
            if cancelled.contains(&op.id) {
                continue;
            }
            active.insert(op.id);
            if let OpKind::Undo { target } = op.kind {
                cancelled.insert(target);
            }
        }
        active
    }

    /// The effective payloads at the current heads, in replay order.
    pub fn effective(&self) -> Vec<(OpId, &P)> {
        let active = self.active_within(self.ops.values());
        self.ops
            .values()
            .filter(|o| active.contains(&o.id))
            .filter_map(|o| match &o.kind {
                OpKind::Apply(p) => Some((o.id, p)),
                OpKind::Undo { .. } => None,
            })
            .collect()
    }

    /// The effective payloads as of an earlier (or concurrent) set of heads.
    /// This is the basis of history restore, versions and branches.
    pub fn effective_at(&self, heads: &[OpId]) -> Vec<(OpId, &P)> {
        let closure = self.causal_closure(heads);
        let subset: Vec<&Op<P>> = self
            .ops
            .values()
            .filter(|o| closure.contains(&o.id))
            .collect();
        let active = self.active_within(subset.iter().copied());
        subset
            .into_iter()
            .filter(|o| active.contains(&o.id))
            .filter_map(|o| match &o.kind {
                OpKind::Apply(p) => Some((o.id, p)),
                OpKind::Undo { .. } => None,
            })
            .collect()
    }

    /// Operations this log has that a peer whose heads are `remote_heads`
    /// lacks, in an order the peer can integrate directly.
    pub fn missing_for(&self, remote_heads: &[OpId]) -> Vec<Op<P>> {
        let have = self.causal_closure(remote_heads);
        self.ops
            .values()
            .filter(|o| !have.contains(&o.id))
            .cloned()
            .collect()
    }

    /// Checks that need nothing but the op itself.
    fn check_self(op: &Op<P>) -> Result<(), OpLogError> {
        if op.id.lamport == 0 {
            return Err(OpLogError::NullId);
        }
        for d in &op.deps {
            if d.lamport >= op.id.lamport {
                return Err(OpLogError::CausalityViolation(op.id));
            }
        }
        Ok(())
    }

    fn check(&self, op: &Op<P>) -> Result<(), OpLogError> {
        Self::check_self(op)?;
        if let OpKind::Undo { target } = op.kind {
            if !self.causal_closure(&op.deps).contains(&target) {
                return Err(OpLogError::BadUndoTarget(op.id));
            }
        }
        Ok(())
    }

    fn try_insert(&mut self, op: Op<P>) -> Result<(), OpLogError> {
        self.check(&op)?;
        for d in &op.deps {
            self.heads.remove(d);
        }
        // An op is a head unless something already integrated depends on it,
        // which cannot happen because dependents are integrated after it.
        self.heads.insert(op.id);
        self.clock = self.clock.max(op.id.lamport);
        self.ops.insert(op.id, op);
        Ok(())
    }

    /// Integrates an operation received from another replica. Operations whose
    /// dependencies have not arrived yet are buffered.
    pub fn integrate(&mut self, op: Op<P>) -> Result<Integrated, OpLogError> {
        if let Some(existing) = self.ops.get(&op.id) {
            return if *existing == op {
                Ok(Integrated::Duplicate)
            } else {
                Err(OpLogError::ConflictingDuplicate(op.id))
            };
        }
        Self::check_self(&op)?;
        if op.deps.iter().any(|d| !self.ops.contains_key(d)) {
            if !self.pending.iter().any(|p| p.id == op.id) {
                self.pending.push(op);
            }
            return Ok(Integrated::Pending);
        }
        self.try_insert(op)?;
        // Drain any buffered ops that are now ready.
        loop {
            let ready = self
                .pending
                .iter()
                .position(|p| p.deps.iter().all(|d| self.ops.contains_key(d)));
            match ready {
                Some(i) => {
                    let p = self.pending.swap_remove(i);
                    if !self.ops.contains_key(&p.id) {
                        let id = p.id;
                        if let Err(e) = self.try_insert(p) {
                            self.rejected.push((id, e));
                        }
                    }
                }
                None => break,
            }
        }
        Ok(Integrated::Added)
    }
}

impl<P: Clone + PartialEq + Serialize + DeserializeOwned> OpLog<P> {
    /// Serializes every integrated operation into the versioned CBOR container
    /// described in `docs/NATIVE_FORMAT.md`.
    pub fn encode(&self) -> Result<Vec<u8>, OpLogError> {
        let c = Container {
            magic: MAGIC.to_string(),
            version: FORMAT_VERSION,
            ops: self.ops.values().cloned().collect(),
        };
        let mut out = Vec::new();
        ciborium::into_writer(&c, &mut out).map_err(|e| OpLogError::Encode(e.to_string()))?;
        Ok(out)
    }

    /// Decodes a log written by [`OpLog::encode`]. Local writes go to
    /// `replica`; pass the replica that wrote the file to resume its undo
    /// history, or a fresh one to start a new editing session.
    pub fn decode(bytes: &[u8], replica: ReplicaId) -> Result<Self, OpLogError> {
        let c: Container<P> =
            ciborium::from_reader(bytes).map_err(|e| OpLogError::Decode(e.to_string()))?;
        if c.magic != MAGIC && c.magic != LEGACY_MAGIC {
            return Err(OpLogError::UnsupportedFormat(c.magic));
        }
        if c.version != FORMAT_VERSION {
            return Err(OpLogError::UnsupportedFormat(format!(
                "version {}",
                c.version
            )));
        }
        let mut log = OpLog::new(replica);
        for op in c.ops {
            log.integrate(op)?;
        }
        if log.pending_len() > 0 {
            return Err(OpLogError::Decode(format!(
                "{} operations have missing dependencies",
                log.pending_len()
            )));
        }
        log.rebuild_undo_stacks();
        Ok(log)
    }
}
