//! Origin tables: how each face, edge and vertex of a body arose from the
//! inputs of the operation that made it (ADR 0018, design note 0013). The
//! naming layer composes them into persistent names; the kernel knows
//! nothing of features or sketches.

use serde::{Deserialize, Serialize};

/// A face of a body: its position in the body's face order, which is also
/// the order of [`crate::BodyMesh::faces`]. Meaningful only for its body.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct FaceIx(pub u32);

/// An edge of a body, in the order of [`crate::BodyMesh::edges`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EdgeIx(pub u32);

/// A vertex of a body.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct VertexIx(pub u32);

/// A piece of a profile curve: the kernel splits arcs and circles into
/// pieces of at most π and closed splines in two (ADR 0008); other curves
/// are one piece.
#[derive(
    Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
)]
pub struct ProfilePiece {
    /// The curve's index in [`crate::Profile::curves`].
    pub curve: u32,
    /// The piece's index along the curve's own direction.
    pub piece: u32,
}

/// One end of a profile piece, in the curve's own direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct PieceEnd {
    /// The piece.
    pub piece: ProfilePiece,
    /// `false` for its start, `true` for its end.
    pub end: bool,
}

/// A vertex of a profile: the two piece ends that meet there, the lesser
/// first (so it does not depend on loop direction or curve order).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ProfileVertex {
    /// The piece ends, sorted.
    pub ends: [PieceEnd; 2],
}

/// One of the two bodies of a boolean (ADR 0020): `a` and `b` of
/// [`crate::Kernel::boolean`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Operand {
    /// The first body, `a`.
    First,
    /// The second body, `b`.
    Second,
}

/// How a face arose.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FaceOrigin {
    /// Swept by a profile piece between sweep positions `span` and
    /// `span + 1`.
    Side {
        /// The piece.
        piece: ProfilePiece,
        /// The span.
        span: u32,
    },
    /// The profile region `region` at sweep position `at` (the first or
    /// the last).
    Cap {
        /// The position.
        at: u32,
        /// The region, an index into [`Origins::regions`].
        region: u32,
    },
    /// A face of a boolean's operand, or a piece of one.
    Operand {
        /// The operand.
        operand: Operand,
        /// Its face.
        face: FaceIx,
    },
}

/// How an edge arose.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum EdgeOrigin {
    /// A copy of a profile piece at sweep position `at`. A piece lying on
    /// a revolve's axis is one edge, at position 0.
    Profile {
        /// The piece.
        piece: ProfilePiece,
        /// The position.
        at: u32,
    },
    /// The line or arc a profile vertex sweeps between positions `span`
    /// and `span + 1`.
    Swept {
        /// The vertex.
        vertex: ProfileVertex,
        /// The span.
        span: u32,
    },
    /// An edge of a boolean's operand, or a piece of one.
    Operand {
        /// The operand.
        operand: Operand,
        /// Its edge.
        edge: EdgeIx,
    },
    /// Where face `first` of a boolean's first operand meets face `second`
    /// of its second.
    Meet {
        /// The first operand's face.
        first: FaceIx,
        /// The second operand's face.
        second: FaceIx,
    },
}

/// How a vertex arose.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum VertexOrigin {
    /// A copy of a profile vertex at sweep position `at`. A vertex on a
    /// revolve's axis is one vertex, at position 0.
    Profile {
        /// The profile vertex.
        vertex: ProfileVertex,
        /// The position.
        at: u32,
    },
    /// A vertex of a boolean's operand.
    Operand {
        /// The operand.
        operand: Operand,
        /// Its vertex.
        vertex: VertexIx,
    },
    /// A new vertex of a boolean, where the operands meet.
    Meet,
}

/// A further origin of an entity made of several (ADR 0020): faces joined
/// where they continue one another, edges joined end to end, an edge or
/// vertex of both operands.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Merged {
    /// A face and one more of its origins.
    Face(FaceIx, FaceOrigin),
    /// An edge and one more of its origins.
    Edge(EdgeIx, EdgeOrigin),
    /// A vertex and one more of its origins.
    Vertex(VertexIx, VertexOrigin),
}

/// The origins of a body's entities, in handle order (ADR 0018).
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Origins {
    /// Per face.
    pub faces: Vec<FaceOrigin>,
    /// Per edge.
    pub edges: Vec<EdgeOrigin>,
    /// Per vertex.
    pub vertices: Vec<VertexOrigin>,
    /// The sweep's positions: 2 for an extrude (start and end), one more
    /// than its pieces for a revolve. A full revolve's last position is its
    /// first, so no entity is at it.
    pub positions: u32,
    /// The profile regions, in the kernel's order: each the pieces of its
    /// loops, outer loop first, each loop in its order of travel
    /// (counter-clockwise about the profile plane's normal for the outer
    /// loop, clockwise for holes).
    pub regions: Vec<Vec<Vec<ProfilePiece>>>,
    /// The further origins of entities made of several, sorted; empty but
    /// for booleans.
    #[serde(default)]
    pub merged: Vec<Merged>,
}

/// Which entities meet which, for naming edges by their faces and
/// vertices by their edges (design note 0013 §4).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Adjacency {
    /// Per edge, the faces it bounds, in face order (two on a closed
    /// solid).
    pub edge_faces: Vec<Vec<FaceIx>>,
    /// Per vertex, the edges that end at it, in edge order.
    pub vertex_edges: Vec<Vec<EdgeIx>>,
    /// Per edge, whether it is a seam: two pieces of one surface meet
    /// there with no crease (as [`crate::BodyMesh::seams`]), so the faces
    /// it parts look like one face and it like no edge.
    pub seams: Vec<bool>,
}
