//! Persistent naming (design note 0013): names for the faces, edges and
//! vertices of a body, made from the kernel's origin table (ADR 0018) and
//! the sketch items its profile came from, and the resolution of stored
//! references against them.
//!
//! One face the user sees can be several kernel faces: the kernel splits
//! periodic faces at π (ADR 0008), and a sweep makes a face of each piece
//! of a profile, so a sketch line split by a vertex gives two faces of one
//! plane. Names belong to *logical* entities, what the user sees and
//! picks: a logical face is the kernel faces joined across seams (edges
//! where pieces of one surface meet with no crease), named by all their
//! sources when those differ; a logical edge is the kernel edges between
//! the same two logical faces joined end to end where nothing else meets
//! them; a logical vertex is where logical edges end. A seam, a vertex
//! inside a logical edge and one on no logical edge at all (a sphere's
//! pole) have no names.
//!
//! Name and reference types are data stored in feature definitions and
//! live in `opendatum_core_model::naming`.

use std::collections::{BTreeMap, BTreeSet};

use opendatum_core_model::naming::{
    Along, EdgeName, FaceName, Hint, Name, ProfilePoint, Reference, Role, Source, Use, VertexName,
};
use opendatum_core_model::sketch::EdgeKey;
use opendatum_core_model::{FeatureId, RegionKey};
use opendatum_kernel::{
    Adjacency, EdgeIx, EdgeOrigin, FaceIx, FaceOrigin, Origins, PieceEnd, ProfileVertex, VertexIx,
    VertexOrigin,
};

/// What made a body, beyond its origin table.
#[derive(Clone, Copy, Debug)]
pub enum Made<'a> {
    /// `Kernel::make_box`: an extrude of the box's rectangle, whose curves
    /// run along −y, +x, +y and −x of the box's frame; `z_positive` when
    /// the box extends along +z from its frame.
    Box {
        /// The box's z size is positive.
        z_positive: bool,
    },
    /// An extrude or revolve of the sketch regions `regions`, whose
    /// profile's curve `i` is the sketch edge `curves[i]` (a closed sketch
    /// edge can be two curves, cut in two by the region profile).
    Sweep {
        /// The sketch edge of each profile curve, as the boundary of the
        /// region it bounds runs along it.
        curves: &'a [EdgeKey],
        /// The regions the feature selects.
        regions: &'a [PickedRegion<'a>],
    },
}

/// A sketch region a feature selects.
#[derive(Clone, Copy, Debug)]
pub struct PickedRegion<'a> {
    /// The key the feature's definition stores for it. Its outer loop
    /// names the caps the region becomes, so a cap keeps its name for as
    /// long as the stored key keeps finding the region (design note 0006
    /// §5), whatever else changes in the sketch.
    pub stored: &'a RegionKey,
    /// The region's key in the sketch as it is now: the sketch edges of
    /// its loops, which the profile is made of.
    pub current: &'a RegionKey,
}

/// The feature that made a body and how.
#[derive(Clone, Copy, Debug)]
pub struct BodySource<'a> {
    /// The feature.
    pub feature: FeatureId,
    /// How.
    pub made: Made<'a>,
}

/// Why a body could not be named.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum NamingError {
    /// The origin table and adjacency do not describe one body, or would
    /// give one name to entities that are not one.
    #[error("origin table does not fit the body: {0}")]
    Inconsistent(String),
    /// A profile curve without a sketch edge.
    #[error("profile curve {0} has no sketch edge")]
    MissingCurve(u32),
    /// An origin this version cannot name.
    #[error("cannot name an entity made as {0}")]
    Unsupported(String),
}

/// A logical entity: its name and its kernel pieces, in handle order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Logical<N, I> {
    /// The name.
    pub name: N,
    /// The kernel pieces.
    pub pieces: Vec<I>,
}

/// A kernel entity of a body.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Pick {
    /// A face.
    Face(FaceIx),
    /// An edge.
    Edge(EdgeIx),
    /// A vertex.
    Vertex(VertexIx),
}

/// The names of a body's logical entities.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BodyNames {
    /// Logical faces, in the order of their first pieces.
    pub faces: Vec<Logical<FaceName, FaceIx>>,
    /// Logical edges, in the order of their first pieces.
    pub edges: Vec<Logical<EdgeName, EdgeIx>>,
    /// Logical vertices, in the order of their first pieces.
    pub vertices: Vec<Logical<VertexName, VertexIx>>,
    face_of: Vec<usize>,
    edge_of: Vec<Option<usize>>,
    vertex_of: Vec<Option<usize>>,
}

/// A resolved reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Resolved {
    /// The kernel pieces of the entity, or of all the entities it now is
    /// (a carrier or set reference to a sketch edge that was split).
    pub pieces: Vec<Pick>,
    /// The stored hint lies on none of them: the entity moved with the
    /// model (informational, design note 0013 §6).
    pub moved: bool,
}

/// Why a reference did not resolve.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ResolveError {
    /// Nothing carries the name or descends from the named entity.
    #[error("no longer exists")]
    Missing,
    /// Several entities descend from the named one and the reference needs
    /// one; the hint lies on `0` of them or on more than one.
    #[error("is ambiguous: {0} candidates")]
    Ambiguous(usize),
}

/// The root of `i`'s set (union–find with path halving).
fn find(p: &mut [usize], mut i: usize) -> usize {
    while p[i] != i {
        p[i] = p[p[i]];
        i = p[i];
    }
    i
}

fn union(p: &mut [usize], a: usize, b: usize) {
    let (a, b) = (find(p, a), find(p, b));
    p[a] = b;
}

/// The name of kernel faces that look like one face: their name if they
/// share it, else, as for faces merged by an operation (design note 0013
/// §5), their feature and role with all their sources. Faces of one
/// operation made in different roles never look like one; if they did,
/// the name could not say so.
fn merged(names: BTreeSet<&FaceName>) -> Result<FaceName, NamingError> {
    let mut all = names.iter();
    let Some(&first) = all.next() else {
        return Err(NamingError::Inconsistent("a face without a name".into()));
    };
    if names.len() == 1 {
        return Ok(first.clone());
    }
    if let Some(other) = all.find(|n| n.feature != first.feature || n.role != first.role) {
        return Err(NamingError::Unsupported(format!(
            "one face made as {first:?} and {other:?}"
        )));
    }
    let mut sources: Vec<Source> = names
        .iter()
        .flat_map(|n| n.sources.iter().cloned())
        .collect();
    sources.sort();
    sources.dedup();
    Ok(FaceName {
        feature: first.feature,
        role: first.role.clone(),
        sources,
    })
}

/// The sketch edge as the region on its other side runs along it.
fn flipped(k: &EdgeKey) -> EdgeKey {
    EdgeKey {
        entity: k.entity,
        forward: !k.forward,
        from: k.to.clone(),
        to: k.from.clone(),
        part: k.part,
    }
}

/// The sources of the caps of each profile region (in the kernel's
/// order): the stored outer loops of the selected regions it covers.
/// Selected regions that share sketch edges become one profile region,
/// bounded by the edges they do not share (less any the region profile
/// leaves out as below its resolution); every profile region must be
/// bounded by edges of exactly one such group, and every group must
/// become exactly one profile region.
fn region_sources(
    origins: &Origins,
    curves: &[EdgeKey],
    regions: &[PickedRegion],
) -> Result<Vec<Vec<Source>>, NamingError> {
    let bad = |m: String| NamingError::Inconsistent(m);
    // Group the selected regions by the sketch edges they share. A key in
    // two entries means the same region selected twice.
    let mut owner: BTreeMap<&EdgeKey, usize> = BTreeMap::new();
    let mut parent: Vec<usize> = (0..regions.len()).collect();
    for (i, r) in regions.iter().enumerate() {
        for k in r.current.loops.iter().flatten() {
            if let Some(&j) = owner.get(k) {
                union(&mut parent, i, j);
            } else {
                owner.insert(k, i);
            }
        }
    }
    for (i, r) in regions.iter().enumerate() {
        for k in r.current.loops.iter().flatten() {
            if let Some(&j) = owner.get(&flipped(k)) {
                union(&mut parent, i, j);
            }
        }
    }
    let mut boundary: BTreeMap<usize, BTreeSet<&EdgeKey>> = BTreeMap::new();
    for (i, r) in regions.iter().enumerate() {
        let g = find(&mut parent, i);
        let b = boundary.entry(g).or_default();
        for k in r.current.loops.iter().flatten() {
            if !owner.contains_key(&flipped(k)) {
                b.insert(k);
            }
        }
    }
    let mut used: BTreeMap<usize, usize> = BTreeMap::new();
    let mut out = Vec::with_capacity(origins.regions.len());
    for (q, loops) in origins.regions.iter().enumerate() {
        let mut keys: BTreeSet<&EdgeKey> = BTreeSet::new();
        for p in loops.iter().flatten() {
            keys.insert(
                curves
                    .get(p.curve as usize)
                    .ok_or(NamingError::MissingCurve(p.curve))?,
            );
        }
        let groups: Vec<usize> = boundary
            .iter()
            .filter(|(_, b)| keys.iter().any(|k| b.contains(k)))
            .map(|(&g, _)| g)
            .collect();
        let [g] = groups[..] else {
            return Err(bad(format!(
                "profile region {q} is bounded by {} groups of selected regions",
                groups.len()
            )));
        };
        if !keys.iter().all(|k| boundary[&g].contains(k)) {
            return Err(bad(format!(
                "profile region {q} has an edge that bounds no selected region"
            )));
        }
        if let Some(p) = used.insert(g, q) {
            return Err(bad(format!(
                "profile regions {p} and {q} come from one group of selected regions"
            )));
        }
        let mut sources = Vec::new();
        for (i, r) in regions.iter().enumerate() {
            if find(&mut parent, i) == g {
                let outer = r
                    .stored
                    .loops
                    .first()
                    .ok_or_else(|| bad(format!("selected region {i} has no stored loop")))?;
                sources.push(Source::Region(outer.clone()));
            }
        }
        sources.sort();
        sources.dedup();
        out.push(sources);
    }
    if used.len() != boundary.len() {
        return Err(bad(format!(
            "{} groups of selected regions for {} profile regions",
            boundary.len(),
            origins.regions.len()
        )));
    }
    Ok(out)
}

/// Chooses `along` for entities named by their faces: none where the faces
/// alone are unique, the origin's form where they are not, an ordinal
/// where that does not tell them apart either.
fn tie_breaks(keys: &[(BTreeSet<usize>, Option<Along>)]) -> Vec<Option<Along>> {
    let mut by_faces: BTreeMap<&BTreeSet<usize>, Vec<usize>> = BTreeMap::new();
    for (i, (f, _)) in keys.iter().enumerate() {
        by_faces.entry(f).or_default().push(i);
    }
    let mut out = vec![None; keys.len()];
    for same in by_faces.values() {
        if same.len() == 1 {
            continue;
        }
        let alongs: Vec<&Option<Along>> = same.iter().map(|&i| &keys[i].1).collect();
        let distinct = alongs.iter().all(|a| a.is_some())
            && alongs.iter().collect::<BTreeSet<_>>().len() == alongs.len();
        for (n, &i) in same.iter().enumerate() {
            out[i] = if distinct {
                keys[i].1.clone()
            } else {
                Some(Along::Ordinal(n as u32))
            };
        }
    }
    out
}

impl BodyNames {
    /// Names the entities of a body from its origin table and adjacency.
    /// Fails, rather than give one name to two entities, when the table,
    /// the adjacency and `source` do not describe one body.
    pub fn new(
        origins: &Origins,
        adjacency: &Adjacency,
        source: &BodySource,
    ) -> Result<BodyNames, NamingError> {
        let n_faces = origins.faces.len() as u32;
        if adjacency.edge_faces.len() != origins.edges.len()
            || adjacency.seams.len() != origins.edges.len()
            || adjacency.vertex_edges.len() != origins.vertices.len()
            || adjacency
                .edge_faces
                .iter()
                .flatten()
                .any(|f| f.0 >= n_faces)
            || adjacency
                .vertex_edges
                .iter()
                .flatten()
                .any(|e| e.0 as usize >= origins.edges.len())
        {
            return Err(NamingError::Inconsistent(format!(
                "the adjacency is not of the table's {} faces, {} edges and {} vertices",
                origins.faces.len(),
                origins.edges.len(),
                origins.vertices.len(),
            )));
        }
        let key = |curve: u32| -> Result<EdgeKey, NamingError> {
            match source.made {
                Made::Sweep { curves, .. } => curves
                    .get(curve as usize)
                    .cloned()
                    .ok_or(NamingError::MissingCurve(curve)),
                Made::Box { .. } => Err(NamingError::MissingCurve(curve)),
            }
        };
        let end = |p: &PieceEnd| -> Result<ProfilePoint, NamingError> {
            Ok(ProfilePoint {
                edge: key(p.piece.curve)?,
                piece: p.piece.piece,
                end: p.end,
            })
        };
        // A profile vertex is where two piece ends meet; its name is the
        // lesser of the two in sketch terms, so it depends neither on the
        // order of the profile's curves nor on the direction of its loops.
        let point = |v: ProfileVertex| -> Result<ProfilePoint, NamingError> {
            let (a, b) = (end(&v.ends[0])?, end(&v.ends[1])?);
            Ok(a.min(b))
        };
        let last = origins.positions.saturating_sub(1);
        let cap_sources = match source.made {
            Made::Sweep { curves, regions } => region_sources(origins, curves, regions)?,
            Made::Box { .. } => Vec::new(),
        };

        // Faces.
        let mut face_names = Vec::with_capacity(origins.faces.len());
        for (i, f) in origins.faces.iter().enumerate() {
            let face = |role: Role, sources: Vec<Source>| FaceName {
                feature: source.feature,
                role,
                sources,
            };
            let name = match (*f, source.made) {
                (FaceOrigin::Side { piece, .. }, Made::Box { .. }) => {
                    let (axis, positive) = match piece.curve {
                        0 => (1, false),
                        1 => (0, true),
                        2 => (1, true),
                        3 => (0, false),
                        c => return Err(NamingError::MissingCurve(c)),
                    };
                    face(Role::BoxSide { axis, positive }, vec![])
                }
                (FaceOrigin::Cap { at, .. }, Made::Box { z_positive }) => face(
                    Role::BoxSide {
                        axis: 2,
                        positive: (at != 0) == z_positive,
                    },
                    vec![],
                ),
                (FaceOrigin::Side { piece, .. }, Made::Sweep { .. }) => {
                    face(Role::Side, vec![Source::SketchEdge(key(piece.curve)?)])
                }
                (FaceOrigin::Cap { at, region }, Made::Sweep { .. }) => {
                    if at != 0 && at != last {
                        return Err(NamingError::Inconsistent(format!(
                            "face {i}: a cap at position {at} of {}",
                            origins.positions
                        )));
                    }
                    let sources = cap_sources.get(region as usize).cloned().ok_or_else(|| {
                        NamingError::Inconsistent(format!("face {i}: no region {region}"))
                    })?;
                    face(Role::Cap { end: at != 0 }, sources)
                }
                (other, _) => return Err(NamingError::Unsupported(format!("{other:?}"))),
            };
            face_names.push(name);
        }
        // Logical faces: the kernel faces joined across seams, which look
        // like one face. Faces with one name must be one logical face.
        let nf = face_names.len();
        let mut fp: Vec<usize> = (0..nf).collect();
        for (e, fs) in adjacency.edge_faces.iter().enumerate() {
            if adjacency.seams[e] {
                for w in fs.windows(2) {
                    union(&mut fp, w[0].0 as usize, w[1].0 as usize);
                }
            }
        }
        let mut face_of = vec![0; nf];
        let mut groups: Vec<Vec<FaceIx>> = Vec::new();
        let mut group_at: BTreeMap<usize, usize> = BTreeMap::new();
        for (f, at) in face_of.iter_mut().enumerate() {
            let root = find(&mut fp, f);
            let g = *group_at.entry(root).or_insert(groups.len());
            if g == groups.len() {
                groups.push(Vec::new());
            }
            groups[g].push(FaceIx(f as u32));
            *at = g;
        }
        let mut named: BTreeMap<&FaceName, usize> = BTreeMap::new();
        for (f, n) in face_names.iter().enumerate() {
            if named.insert(n, face_of[f]).is_some_and(|g| g != face_of[f]) {
                return Err(NamingError::Inconsistent(format!(
                    "separate faces made as {n:?}"
                )));
            }
        }
        let faces = groups
            .into_iter()
            .map(|pieces| {
                let names = pieces.iter().map(|f| &face_names[f.0 as usize]).collect();
                Ok(Logical {
                    name: merged(names)?,
                    pieces,
                })
            })
            .collect::<Result<Vec<_>, NamingError>>()?;

        // Logical edges: kernel edges between two logical faces (others are
        // seams, inside one), joined end to end where no other edge meets
        // them.
        let ne = origins.edges.len();
        let face_sets: Vec<BTreeSet<usize>> = adjacency
            .edge_faces
            .iter()
            .map(|fs| fs.iter().map(|f| face_of[f.0 as usize]).collect())
            .collect();
        let is_edge = |e: usize| face_sets[e].len() >= 2;
        let mut ep: Vec<usize> = (0..ne).collect();
        for es in &adjacency.vertex_edges {
            let at: Vec<usize> = es
                .iter()
                .map(|e| e.0 as usize)
                .filter(|&e| is_edge(e))
                .collect();
            if let [a, b] = at[..] {
                if face_sets[a] == face_sets[b] {
                    union(&mut ep, a, b);
                }
            }
        }
        let mut edge_of = vec![None; ne];
        let mut chains: Vec<Vec<EdgeIx>> = Vec::new();
        let mut chain_at: BTreeMap<usize, usize> = BTreeMap::new();
        for (e, at) in edge_of.iter_mut().enumerate() {
            if !is_edge(e) {
                continue;
            }
            let root = find(&mut ep, e);
            let g = *chain_at.entry(root).or_insert(chains.len());
            if g == chains.len() {
                chains.push(Vec::new());
            }
            chains[g].push(EdgeIx(e as u32));
            *at = Some(g);
        }
        let mut keys = Vec::with_capacity(chains.len());
        for pieces in &chains {
            // The least of the pieces' origins, in sketch terms.
            let mut along: Option<Along> = None;
            for p in pieces {
                let a = match (origins.edges[p.0 as usize], source.made) {
                    (_, Made::Box { .. }) => continue,
                    (EdgeOrigin::Profile { piece, at }, _) => Along::Profile {
                        edge: key(piece.curve)?,
                        end: at != 0,
                    },
                    (EdgeOrigin::Swept { vertex, .. }, _) => Along::Swept(point(vertex)?),
                    (other, _) => return Err(NamingError::Unsupported(format!("{other:?}"))),
                };
                along = Some(match along {
                    Some(b) => a.min(b),
                    None => a,
                });
            }
            keys.push((face_sets[pieces[0].0 as usize].clone(), along));
        }
        let alongs = tie_breaks(&keys);
        let edges: Vec<Logical<EdgeName, EdgeIx>> = chains
            .into_iter()
            .zip(keys)
            .zip(alongs)
            .map(|((pieces, (fs, _)), along)| Logical {
                name: EdgeName {
                    faces: fs.iter().map(|&f| faces[f].name.clone()).collect(),
                    along,
                },
                pieces,
            })
            .collect();

        // Logical vertices: where logical edges end, not where one runs
        // through.
        let mut vertex_keys = Vec::new();
        let mut vertex_ixs = Vec::new();
        for (v, o) in origins.vertices.iter().enumerate() {
            let incident = &adjacency.vertex_edges[v];
            let at: Vec<usize> = incident
                .iter()
                .filter_map(|e| edge_of[e.0 as usize])
                .collect();
            if at.is_empty() || (at.len() == 2 && at[0] == at[1]) {
                continue;
            }
            let fs: BTreeSet<usize> = incident
                .iter()
                .flat_map(|e| face_sets[e.0 as usize].iter().copied())
                .collect();
            let along = match (*o, source.made) {
                (_, Made::Box { .. }) => None,
                (VertexOrigin::Profile { vertex, at }, _) => Some(Along::Point {
                    point: point(vertex)?,
                    end: at != 0,
                }),
                (other, _) => return Err(NamingError::Unsupported(format!("{other:?}"))),
            };
            vertex_keys.push((fs, along));
            vertex_ixs.push(VertexIx(v as u32));
        }
        let alongs = tie_breaks(&vertex_keys);
        let mut vertex_of = vec![None; origins.vertices.len()];
        let mut vertices = Vec::with_capacity(vertex_ixs.len());
        for (i, ((v, (fs, _)), along)) in vertex_ixs
            .into_iter()
            .zip(vertex_keys)
            .zip(alongs)
            .enumerate()
        {
            vertex_of[v.0 as usize] = Some(i);
            vertices.push(Logical {
                name: VertexName {
                    faces: fs.iter().map(|&f| faces[f].name.clone()).collect(),
                    along,
                },
                pieces: vec![v],
            });
        }

        let names = BodyNames {
            faces,
            edges,
            vertices,
            face_of,
            edge_of,
            vertex_of,
        };
        names.check_unique()?;
        Ok(names)
    }

    /// No two logical entities share a name.
    fn check_unique(&self) -> Result<(), NamingError> {
        let mut seen = BTreeSet::new();
        for (n, _) in self.entities() {
            if !seen.insert(n.clone()) {
                return Err(NamingError::Inconsistent(format!(
                    "two entities named {n:?}"
                )));
            }
        }
        Ok(())
    }

    /// The logical face a kernel face belongs to.
    pub fn face(&self, f: FaceIx) -> &Logical<FaceName, FaceIx> {
        &self.faces[self.face_of[f.0 as usize]]
    }

    /// The logical edge a kernel edge belongs to (none for a seam).
    pub fn edge(&self, e: EdgeIx) -> Option<&Logical<EdgeName, EdgeIx>> {
        self.edge_of[e.0 as usize].map(|g| &self.edges[g])
    }

    /// The logical vertex a kernel vertex is (none inside an edge).
    pub fn vertex(&self, v: VertexIx) -> Option<&Logical<VertexName, VertexIx>> {
        self.vertex_of[v.0 as usize].map(|g| &self.vertices[g])
    }

    /// The name of what a pick selects: the logical entity of the kernel
    /// entity (none for a seam or a point inside an edge).
    pub fn name(&self, pick: Pick) -> Option<Name> {
        match pick {
            Pick::Face(f) => Some(Name::Face(self.face(f).name.clone())),
            Pick::Edge(e) => self.edge(e).map(|l| Name::Edge(l.name.clone())),
            Pick::Vertex(v) => self.vertex(v).map(|l| Name::Vertex(l.name.clone())),
        }
    }

    /// Every logical entity with its name and pieces.
    fn entities(&self) -> impl Iterator<Item = (Name, Vec<Pick>)> + '_ {
        let faces = self.faces.iter().map(|l| {
            (
                Name::Face(l.name.clone()),
                l.pieces.iter().map(|&p| Pick::Face(p)).collect(),
            )
        });
        let edges = self.edges.iter().map(|l| {
            (
                Name::Edge(l.name.clone()),
                l.pieces.iter().map(|&p| Pick::Edge(p)).collect(),
            )
        });
        let vertices = self.vertices.iter().map(|l| {
            (
                Name::Vertex(l.name.clone()),
                l.pieces.iter().map(|&p| Pick::Vertex(p)).collect(),
            )
        });
        faces.chain(edges).chain(vertices)
    }

    /// Resolves a reference (design note 0013 §6): by its exact name; else
    /// among the entities that descend from the named one (names that
    /// differ only where sketch edges were split or extended, keeping
    /// their curves, or where selected regions merged into one cap or
    /// parted): one is the entity; several are its pieces, all of them for
    /// a carrier or set reference and the one the hint lies on for a
    /// single one; none means it no longer exists. Never another entity
    /// that merely shares the feature and role. A name with an ordinal
    /// tie-break is matched like a geometric tie-break: the candidate the
    /// hint lies on, if exactly one. `lies_on` says whether the hint lies
    /// on a kernel entity. Splits and merges by later operations (§5) come
    /// with booleans.
    pub fn resolve(
        &self,
        reference: &Reference,
        lies_on: &dyn Fn(Pick, &Hint) -> bool,
    ) -> Result<Resolved, ResolveError> {
        let hint = &reference.hint;
        let on = |pieces: &[Pick]| pieces.iter().any(|&p| lies_on(p, hint));
        let exact = self.entities().find(|(n, _)| *n == reference.name);
        if weak(&reference.name) {
            let similar: Vec<Vec<Pick>> = self
                .entities()
                .filter(|(n, _)| same_but_ordinal(n, &reference.name))
                .map(|(_, p)| p)
                .collect();
            let hit: Vec<&Vec<Pick>> = similar.iter().filter(|p| on(p)).collect();
            if hit.len() == 1 {
                return Ok(Resolved {
                    pieces: hit[0].clone(),
                    moved: false,
                });
            }
            if hit.len() > 1 {
                return Err(ResolveError::Ambiguous(hit.len()));
            }
        }
        if let Some((_, pieces)) = exact {
            let moved = !on(&pieces);
            return Ok(Resolved { pieces, moved });
        }
        let candidates: Vec<Vec<Pick>> = self
            .entities()
            .filter(|(n, _)| same_items(n, &reference.name))
            .map(|(_, p)| p)
            .collect();
        match candidates.len() {
            0 => Err(ResolveError::Missing),
            // The one entity descending from the referenced one (a sketch
            // edge extended or shortened).
            1 => {
                let pieces = candidates.into_iter().next().unwrap_or_default();
                let moved = !on(&pieces);
                Ok(Resolved { pieces, moved })
            }
            // Pieces of the referenced entity (a sketch edge split).
            n => match reference.usage {
                Use::Carrier | Use::Set => {
                    let pieces: Vec<Pick> = candidates.into_iter().flatten().collect();
                    let moved = !on(&pieces);
                    Ok(Resolved { pieces, moved })
                }
                Use::Single => {
                    let mut hit: Vec<Vec<Pick>> =
                        candidates.into_iter().filter(|p| on(p)).collect();
                    if hit.len() == 1 {
                        Ok(Resolved {
                            pieces: hit.pop().unwrap_or_default(),
                            moved: false,
                        })
                    } else {
                        Err(ResolveError::Ambiguous(n))
                    }
                }
            },
        }
    }
}

/// The name ends in an ordinal tie-break.
fn weak(name: &Name) -> bool {
    matches!(
        name,
        Name::Edge(EdgeName {
            along: Some(Along::Ordinal(_)),
            ..
        }) | Name::Vertex(VertexName {
            along: Some(Along::Ordinal(_)),
            ..
        })
    )
}

/// Equal apart from ordinal tie-breaks.
fn same_but_ordinal(a: &Name, b: &Name) -> bool {
    let ordinal = |x: &Option<Along>| matches!(x, Some(Along::Ordinal(_)));
    match (a, b) {
        (Name::Edge(a), Name::Edge(b)) => {
            a.faces == b.faces && ordinal(&a.along) && ordinal(&b.along)
        }
        (Name::Vertex(a), Name::Vertex(b)) => {
            a.faces == b.faces && ordinal(&a.along) && ordinal(&b.along)
        }
        _ => false,
    }
}

/// Two sketch edges are pieces of one sketch item: the same curve (and
/// text part) in the same direction, whatever their end vertices.
fn same_item(a: &EdgeKey, b: &EdgeKey) -> bool {
    a.entity == b.entity && a.forward == b.forward && a.part == b.part
}

/// Two sources are, or descend from, the same thing: a sketch edge by
/// its item; a selected region by its stored key, or, once the region was
/// selected anew, by the items of its outer loop.
fn same_source(a: &Source, b: &Source) -> bool {
    match (a, b) {
        (Source::SketchEdge(a), Source::SketchEdge(b)) => same_item(a, b),
        (Source::Region(a), Source::Region(b)) => {
            let items = |l: &[EdgeKey]| -> BTreeSet<_> {
                l.iter().map(|k| (k.entity, k.forward, k.part)).collect()
            };
            a == b || items(a) == items(b)
        }
        _ => false,
    }
}

/// One face descends from the other: the same feature and role, made
/// from one of the same sketch items or selected regions. Faces merge
/// when they come to look like one (the pieces of a sketch edge split by
/// a vertex, selected regions that come to share an edge) and part when
/// they no longer do, so any shared source will do.
fn same_face(a: &FaceName, b: &FaceName) -> bool {
    a.feature == b.feature
        && a.role == b.role
        && ((a.sources.is_empty() && b.sources.is_empty())
            || a.sources
                .iter()
                .any(|x| b.sources.iter().any(|y| same_source(x, y))))
}

/// Every face of one set descends from a face of the other, and the sets
/// are the same size.
fn same_faces(a: &BTreeSet<FaceName>, b: &BTreeSet<FaceName>) -> bool {
    a.len() == b.len()
        && a.iter().all(|x| b.iter().any(|y| same_face(x, y)))
        && b.iter().all(|y| a.iter().any(|x| same_face(x, y)))
}

fn same_point(a: &ProfilePoint, b: &ProfilePoint) -> bool {
    same_item(&a.edge, &b.edge) && a.piece == b.piece && a.end == b.end
}

/// Tie-breaks that do not tell the entities apart: either is absent (a
/// tie that came or went with the edit) or weak, or both say the same in
/// terms of sketch items.
fn same_along(a: &Option<Along>, b: &Option<Along>) -> bool {
    match (a, b) {
        (None, _) | (_, None) => true,
        (Some(Along::Ordinal(_)), _) | (_, Some(Along::Ordinal(_))) => true,
        (Some(Along::Profile { edge: a, end: x }), Some(Along::Profile { edge: b, end: y })) => {
            same_item(a, b) && x == y
        }
        (Some(Along::Swept(a)), Some(Along::Swept(b))) => same_point(a, b),
        (Some(Along::Point { point: a, end: x }), Some(Along::Point { point: b, end: y })) => {
            same_point(a, b) && x == y
        }
        _ => false,
    }
}

/// One entity descends from the other: names that differ at most where
/// sketch edges were split or extended, or selected regions merged or
/// parted.
fn same_items(a: &Name, b: &Name) -> bool {
    match (a, b) {
        (Name::Face(a), Name::Face(b)) => same_face(a, b),
        (Name::Edge(a), Name::Edge(b)) => {
            same_faces(&a.faces, &b.faces) && same_along(&a.along, &b.along)
        }
        (Name::Vertex(a), Name::Vertex(b)) => {
            same_faces(&a.faces, &b.faces) && same_along(&a.along, &b.along)
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use opendatum_core_model::sketch::{SketchId, VertexKey};
    use opendatum_kernel::DVec3;
    use opendatum_oplog::{OpId, ReplicaId};

    fn op(n: u64) -> OpId {
        OpId {
            lamport: n,
            replica: ReplicaId(1),
        }
    }

    fn key(i: u32) -> EdgeKey {
        let vk = VertexKey {
            entities: vec![],
            ordinal: i,
            parts: vec![],
        };
        EdgeKey {
            entity: SketchId {
                op: op(2),
                index: i,
            },
            forward: true,
            from: vk.clone(),
            to: vk,
            part: 0,
        }
    }

    fn face(i: u32) -> FaceName {
        FaceName {
            feature: FeatureId(op(1)),
            role: Role::Side,
            sources: vec![Source::SketchEdge(key(i))],
        }
    }

    #[test]
    fn tie_breaks_are_origins_where_they_differ_else_ordinals() {
        let fs = |v: &[usize]| v.iter().copied().collect::<BTreeSet<usize>>();
        let at = |i: u32, end| Along::Profile { edge: key(i), end };
        let keys = vec![
            (fs(&[0, 1]), Some(at(0, false))),
            (fs(&[0, 2]), Some(at(1, false))),
            (fs(&[0, 2]), Some(at(1, true))),
            (fs(&[1, 2]), Some(at(2, false))),
            (fs(&[1, 2]), Some(at(2, false))),
            (fs(&[1, 2]), Some(at(3, false))),
            (fs(&[2, 3]), Some(at(4, false))),
            (fs(&[2, 3]), None),
        ];
        assert_eq!(
            tie_breaks(&keys),
            vec![
                None,
                Some(at(1, false)),
                Some(at(1, true)),
                Some(Along::Ordinal(0)),
                Some(Along::Ordinal(1)),
                Some(Along::Ordinal(2)),
                Some(Along::Ordinal(0)),
                Some(Along::Ordinal(1)),
            ]
        );
    }

    /// Two edges between faces 0 and 1 that only ordinals tell apart; a
    /// hint at x = 0 lies on edge 0, at x = 1 on edge 1, at x = 2 on both.
    fn ordinal_body() -> (BodyNames, impl Fn(Pick, &Hint) -> bool) {
        let edge = |n| Logical {
            name: EdgeName {
                faces: [face(0), face(1)].into(),
                along: Some(Along::Ordinal(n)),
            },
            pieces: vec![EdgeIx(n)],
        };
        let names = BodyNames {
            faces: vec![
                Logical {
                    name: face(0),
                    pieces: vec![FaceIx(0)],
                },
                Logical {
                    name: face(1),
                    pieces: vec![FaceIx(1)],
                },
            ],
            edges: vec![edge(0), edge(1)],
            vertices: vec![],
            face_of: vec![0, 1],
            edge_of: vec![Some(0), Some(1)],
            vertex_of: vec![],
        };
        let lies_on = |p: Pick, h: &Hint| match p {
            Pick::Edge(EdgeIx(e)) => h.point.x == e as f64 || h.point.x == 2.0,
            _ => false,
        };
        (names, lies_on)
    }

    fn reference(along: u32, x: f64, usage: Use) -> Reference {
        Reference {
            name: Name::Edge(EdgeName {
                faces: [face(0), face(1)].into(),
                along: Some(Along::Ordinal(along)),
            }),
            hint: Hint {
                point: DVec3::new(x, 0.0, 0.0),
                direction: DVec3::ZERO,
            },
            usage,
        }
    }

    #[test]
    fn ordinals_are_weak() {
        let (names, lies_on) = ordinal_body();
        let edge = |e| Pick::Edge(EdgeIx(e));
        let got = |r: &Reference| names.resolve(r, &lies_on);
        // The hint wins over the ordinal.
        let r = got(&reference(0, 1.0, Use::Single)).unwrap();
        assert_eq!((r.pieces, r.moved), (vec![edge(1)], false));
        // No hint on either: the ordinal, reported as moved.
        let r = got(&reference(0, 5.0, Use::Single)).unwrap();
        assert_eq!((r.pieces, r.moved), (vec![edge(0)], true));
        // A hint on both decides nothing.
        assert_eq!(
            got(&reference(1, 2.0, Use::Single)),
            Err(ResolveError::Ambiguous(2))
        );
        // An ordinal that no longer exists: the hint, else both for a
        // carrier or set and nothing for a single entity.
        let r = got(&reference(7, 1.0, Use::Single)).unwrap();
        assert_eq!(r.pieces, vec![edge(1)]);
        let r = got(&reference(7, 5.0, Use::Set)).unwrap();
        assert_eq!(r.pieces, vec![edge(0), edge(1)]);
        assert_eq!(
            got(&reference(7, 5.0, Use::Single)),
            Err(ResolveError::Ambiguous(2))
        );
    }
}
