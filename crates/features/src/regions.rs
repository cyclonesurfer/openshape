//! Naming sketch regions so that profile selections survive sketch edits
//! (design note 0006 §5; the sketch-level case of ADR 0003).

use std::collections::{BTreeMap, BTreeSet};

use opendatum_core_model::sketch::{EdgeKey, VertexKey};
use opendatum_core_model::{RegionKey, SketchId};
use opendatum_kernel::{Kernel, ProfileEdge, RegionLoop, SketchArrangement};

use crate::SketchRegions;

/// The key of every region of `arrangement`, whose curve `i` belongs to the
/// sketch curve `owners[i]`. An entity with several curves (text) numbers
/// them in order as parts (ADR 0015).
pub fn region_keys(arrangement: &SketchArrangement, owners: &[SketchId]) -> Vec<RegionKey> {
    let parts = parts_of(owners);
    let vertices = vertex_keys(arrangement, owners, &parts);
    arrangement
        .regions
        .iter()
        .map(|r| {
            let mut holes: Vec<Vec<EdgeKey>> = r
                .holes
                .iter()
                .map(|h| loop_key(h, arrangement, owners, &parts, &vertices))
                .collect();
            holes.sort();
            let mut loops = vec![loop_key(&r.outer, arrangement, owners, &parts, &vertices)];
            loops.extend(holes);
            RegionKey {
                loops,
                interior: r.interior_point,
            }
        })
        .collect()
}

/// Each curve's index among the curves of its owner.
fn parts_of(owners: &[SketchId]) -> Vec<u32> {
    let mut count: BTreeMap<SketchId, u32> = BTreeMap::new();
    owners
        .iter()
        .map(|o| {
            let n = count.entry(*o).or_insert(0);
            *n += 1;
            *n - 1
        })
        .collect()
}

/// Names every vertex by the sketch curves (and parts) through it, with an
/// ordinal for vertices sharing the same curves, ordered along the first of
/// them.
fn vertex_keys(
    arrangement: &SketchArrangement,
    owners: &[SketchId],
    parts: &[u32],
) -> Vec<VertexKey> {
    let sets: Vec<Vec<(SketchId, u32)>> = arrangement
        .vertices
        .iter()
        .map(|v| {
            let mut s: Vec<(SketchId, u32)> =
                v.curves.iter().map(|&c| (owners[c], parts[c])).collect();
            s.sort();
            s.dedup();
            s
        })
        .collect();
    let mut groups: BTreeMap<&[(SketchId, u32)], Vec<usize>> = BTreeMap::new();
    for (i, s) in sets.iter().enumerate() {
        groups.entry(s.as_slice()).or_default().push(i);
    }
    let mut ordinal = vec![0u32; sets.len()];
    for (set, mut members) in groups {
        if members.len() > 1 {
            let first = set.first().copied();
            let along = |v: usize| -> f64 {
                let vx = &arrangement.vertices[v];
                vx.curves
                    .iter()
                    .zip(&vx.params)
                    .find(|(c, _)| Some((owners[**c], parts[**c])) == first)
                    .map(|(_, t)| *t)
                    .unwrap_or(0.0)
            };
            members.sort_by(|&a, &b| {
                let (pa, pb) = (arrangement.vertices[a].point, arrangement.vertices[b].point);
                along(a)
                    .total_cmp(&along(b))
                    .then(pa.x.total_cmp(&pb.x))
                    .then(pa.y.total_cmp(&pb.y))
            });
        }
        for (k, v) in members.into_iter().enumerate() {
            ordinal[v] = k as u32;
        }
    }
    sets.into_iter()
        .zip(ordinal)
        .map(|(set, ordinal)| {
            let entities = set.iter().map(|(e, _)| *e).collect();
            let parts = if set.iter().any(|(_, p)| *p != 0) {
                set.iter().map(|(_, p)| *p).collect()
            } else {
                Vec::new()
            };
            VertexKey {
                entities,
                ordinal,
                parts,
            }
        })
        .collect()
}

/// The key of each arrangement edge of a region profile (in the direction
/// the profile runs along it), as [`region_keys`] names the same edge in
/// the loops of the regions it bounds: the sketch edge behind each curve
/// of [`Kernel::region_profile_edges`].
pub fn edge_keys(
    arrangement: &SketchArrangement,
    owners: &[SketchId],
    edges: &[ProfileEdge],
) -> Vec<EdgeKey> {
    let parts = parts_of(owners);
    let vertices = vertex_keys(arrangement, owners, &parts);
    edges
        .iter()
        .map(|e| edge_key(arrangement, owners, &parts, &vertices, e.edge, e.forward))
        .collect()
}

fn edge_key(
    arrangement: &SketchArrangement,
    owners: &[SketchId],
    parts: &[u32],
    vertices: &[VertexKey],
    e: usize,
    forward: bool,
) -> EdgeKey {
    let edge = &arrangement.edges[e];
    let (from, to) = if forward {
        (edge.start, edge.end)
    } else {
        (edge.end, edge.start)
    };
    EdgeKey {
        entity: owners[edge.curve],
        forward,
        from: vertices[from].clone(),
        to: vertices[to].clone(),
        part: parts[edge.curve],
    }
}

fn loop_key(
    l: &RegionLoop,
    arrangement: &SketchArrangement,
    owners: &[SketchId],
    parts: &[u32],
    vertices: &[VertexKey],
) -> Vec<EdgeKey> {
    let mut keys: Vec<EdgeKey> = l
        .edges
        .iter()
        .map(|&(e, forward)| edge_key(arrangement, owners, parts, vertices, e, forward))
        .collect();
    // Canonical rotation: start at the smallest key.
    if let Some(start) = (0..keys.len()).min_by(|&a, &b| keys[a].cmp(&keys[b])) {
        keys.rotate_left(start);
    }
    keys
}

fn entities(key: &RegionKey) -> BTreeSet<SketchId> {
    key.loops.iter().flatten().map(|e| e.entity).collect()
}

/// Finds the region a stored key refers to: an exact key match, otherwise
/// the region containing the key's interior point if it is still bounded by
/// at least one of the same sketch curves. `None` means the region no longer
/// exists (never a guess).
pub fn resolve_region<K: Kernel>(
    kernel: &K,
    regions: &SketchRegions,
    key: &RegionKey,
) -> Option<usize> {
    if let Some(i) = regions.keys.iter().position(|k| k.loops == key.loops) {
        return Some(i);
    }
    let i = kernel.region_containing(&regions.arrangement, &regions.curves, key.interior)?;
    let wanted = entities(key);
    entities(&regions.keys[i])
        .intersection(&wanted)
        .next()
        .map(|_| i)
}
