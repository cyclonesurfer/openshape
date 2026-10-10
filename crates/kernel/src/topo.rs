//! Boundary-representation topology.
//!
//! Entities live in flat arenas inside a [`Brep`] and refer to each other by
//! index. A solid owns shells; a shell owns faces; a face owns loops (the
//! first is the outer loop); a loop is a cycle of coedges; a coedge is one
//! use of an edge in a given direction; an edge is a bounded piece of a curve
//! between two vertices.
//!
//! Orientation conventions:
//! * A face's outward normal is the surface normal `S_u × S_v`, negated when
//!   `forward` is false. In a solid, outward normals point out of the material.
//! * Loops run counter-clockwise when seen from outside (against the outward
//!   normal): material is on the left of every coedge.
//! * In a closed manifold shell every edge is used by exactly two coedges, in
//!   opposite directions.

use glam::DVec3;
use serde::{Deserialize, Serialize};

use crate::geom::{Curve3, Surface};

macro_rules! id_type {
    ($name:ident) => {
        #[derive(
            Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
        )]
        pub(crate) struct $name(pub u32);
        impl $name {
            pub fn idx(self) -> usize {
                self.0 as usize
            }
        }
    };
}

id_type!(VertexId);
id_type!(EdgeId);
id_type!(FaceId);
id_type!(ShellId);

/// A vertex.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Vertex {
    pub point: DVec3,
}

/// An edge: `curve` restricted to `[t0, t1]`, running from `v0` to `v1`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Edge {
    pub curve: Curve3,
    pub t0: f64,
    pub t1: f64,
    pub v0: VertexId,
    pub v1: VertexId,
}

/// One directed use of an edge in a loop.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Coedge {
    pub edge: EdgeId,
    /// True if the loop traverses the edge from `v0` to `v1`.
    pub forward: bool,
}

/// A closed cycle of coedges.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Loop {
    pub coedges: Vec<Coedge>,
}

/// A bounded region of a surface.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Face {
    pub surface: Surface,
    /// True if the outward normal equals `S_u × S_v`.
    pub forward: bool,
    /// Outer loop first, then inner loops (holes).
    pub loops: Vec<Loop>,
}

/// A connected set of faces.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Shell {
    pub faces: Vec<FaceId>,
}

/// A solid region: an outer shell and optional void shells.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Solid {
    pub shells: Vec<ShellId>,
}

/// A complete B-rep body (possibly several disjoint solids, or lumps).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub(crate) struct Brep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub solids: Vec<Solid>,
}

impl Brep {
    pub fn add_vertex(&mut self, point: DVec3) -> VertexId {
        self.vertices.push(Vertex { point });
        VertexId(self.vertices.len() as u32 - 1)
    }

    pub fn add_edge(
        &mut self,
        curve: Curve3,
        t0: f64,
        t1: f64,
        v0: VertexId,
        v1: VertexId,
    ) -> EdgeId {
        self.edges.push(Edge {
            curve,
            t0,
            t1,
            v0,
            v1,
        });
        EdgeId(self.edges.len() as u32 - 1)
    }

    pub fn add_face(&mut self, surface: Surface, forward: bool, loops: Vec<Loop>) -> FaceId {
        self.faces.push(Face {
            surface,
            forward,
            loops,
        });
        FaceId(self.faces.len() as u32 - 1)
    }

    pub fn add_shell(&mut self, faces: Vec<FaceId>) -> ShellId {
        self.shells.push(Shell { faces });
        ShellId(self.shells.len() as u32 - 1)
    }

    pub fn add_solid(&mut self, shells: Vec<ShellId>) {
        self.solids.push(Solid { shells });
    }

    pub fn vertex(&self, id: VertexId) -> &Vertex {
        &self.vertices[id.idx()]
    }

    pub fn edge(&self, id: EdgeId) -> &Edge {
        &self.edges[id.idx()]
    }

    /// Start vertex of a coedge in its loop direction.
    pub fn coedge_start(&self, c: Coedge) -> VertexId {
        let e = self.edge(c.edge);
        if c.forward {
            e.v0
        } else {
            e.v1
        }
    }

    /// End vertex of a coedge in its loop direction.
    pub fn coedge_end(&self, c: Coedge) -> VertexId {
        let e = self.edge(c.edge);
        if c.forward {
            e.v1
        } else {
            e.v0
        }
    }
}
