//! Persistent names of faces, edges and vertices, and references to them
//! (design note 0013). Names are data stored in feature definitions, so
//! they live here beside the sketch keys they are built from; the
//! `opendatum-naming` crate makes them from the kernel's origin tables and
//! resolves references against them.
//!
//! A name says how an entity was made, never where it is or what index it
//! has: the feature, the role in that feature's operation and the sketch
//! items it was made from. Names identify *logical* entities: the kernel
//! splits periodic faces at π (ADR 0008), and the pieces of one face, edge
//! or vertex that the user sees as one share one name.

use std::collections::BTreeSet;

use glam::DVec3;
use serde::{Deserialize, Serialize};

use crate::sketch::EdgeKey;
use crate::FeatureId;

/// What a feature's operation made a face as.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Role {
    /// A side of a box: its axis in the box's frame (0 = x, 1 = y, 2 = z)
    /// and whether it faces along it.
    BoxSide {
        /// The axis.
        axis: u8,
        /// Facing along the axis.
        positive: bool,
    },
    /// The face a sketch edge sweeps (extrude, revolve).
    Side,
    /// A cap of a sweep, at its start or its end.
    Cap {
        /// The end cap (else the start cap).
        end: bool,
    },
}

/// What a face was made from.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Source {
    /// A sketch edge (the piece of a sketch curve between two vertices).
    SketchEdge(EdgeKey),
    /// A sketch region the feature selects, by the outer loop of the key
    /// the feature stores for it (design note 0006 §5): its sketch edges,
    /// starting at the smallest. A cap carries one per selected region it
    /// covers (regions that share edges become one cap).
    Region(Vec<EdgeKey>),
}

/// The name of a face: how it was made.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct FaceName {
    /// The feature that made it.
    pub feature: FeatureId,
    /// Its role in the feature's operation.
    pub role: Role,
    /// What it was made from.
    pub sources: Vec<Source>,
}

/// A point of a profile: the start or end of a piece of a sketch edge,
/// pieces counted along the sketch curve's own direction (the kernel
/// splits arcs at most every π, closed splines in two; ADR 0018).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ProfilePoint {
    /// The sketch edge.
    pub edge: EdgeKey,
    /// The piece of it.
    pub piece: u32,
    /// The piece's end (else its start).
    pub end: bool,
}

/// Where along its faces an edge or vertex lies, when its faces alone do
/// not tell it from another (design note 0013 §4): the sketch-level form
/// of the kernel's origin of it.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Along {
    /// The copy of a sketch edge at the start (`end` false) or end of the
    /// sweep; for a partial revolve, an edge on its axis.
    Profile {
        /// The sketch edge.
        edge: EdgeKey,
        /// At the sweep's end (else its start).
        end: bool,
    },
    /// The edge a profile point sweeps.
    Swept(ProfilePoint),
    /// A profile point at the start (`end` false) or end of the sweep; for
    /// a revolve, a point on its axis.
    Point {
        /// The point.
        point: ProfilePoint,
        /// At the sweep's end (else its start).
        end: bool,
    },
    /// The position among several edges or vertices that nothing else
    /// tells apart, in the order of their kernel origins. Weak: resolution
    /// treats it like a geometric tie-break (design note 0013 §4).
    Ordinal(u32),
}

/// The name of an edge: the faces it bounds, plus where along them when
/// two edges bound the same faces.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EdgeName {
    /// The two faces (one for an edge a single face bounds on both sides).
    pub faces: BTreeSet<FaceName>,
    /// The tie-break, when needed.
    pub along: Option<Along>,
}

/// The name of a vertex: the faces around it, plus where when two
/// vertices have the same faces.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct VertexName {
    /// The faces around it.
    pub faces: BTreeSet<FaceName>,
    /// The tie-break, when needed.
    pub along: Option<Along>,
}

/// The name of a face, edge or vertex.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Name {
    /// A face.
    Face(FaceName),
    /// An edge.
    Edge(EdgeName),
    /// A vertex.
    Vertex(VertexName),
}

/// What a referring feature needs of the entity (design note 0013 §6).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Use {
    /// Only the surface or curve the entity lies on (sketch on a face,
    /// extrude up to a face, a direction along an edge): any piece of a
    /// split entity will do.
    Carrier,
    /// The entity as selected, all its pieces (edges to fillet, faces to
    /// remove).
    Set,
    /// One bounded entity: of a split entity, the piece the hint lies on.
    Single,
}

/// Where the user picked the entity, in model coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Hint {
    /// The picked point.
    pub point: DVec3,
    /// The surface normal (face) or curve tangent (edge) there; zero for a
    /// vertex.
    pub direction: DVec3,
}

/// A feature's reference to an entity of an earlier feature's body.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Reference {
    /// The entity's name when it was picked.
    pub name: Name,
    /// Where it was picked.
    pub hint: Hint,
    /// What the referring feature needs of it.
    pub usage: Use,
}
