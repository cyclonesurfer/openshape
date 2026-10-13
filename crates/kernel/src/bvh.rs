//! A bounding-volume hierarchy: a binary tree of axis-aligned boxes over a
//! set of items, so the kernel's spatial queries (which items' boxes meet a
//! box, a segment, or the neighbourhood of a point) cost a logarithmic walk
//! rather than a pass over every item.
//!
//! The tree is built top down, splitting each node's items at the median of
//! their box centres along the widest extent of those centres (Ericson,
//! *Real-Time Collision Detection*, 2005, §6.1–6.2), so its depth is
//! logarithmic whatever the input. The build is deterministic: ties are
//! broken by item index, so equal inputs give equal trees and every query
//! visits items in the same order on every machine.
//!
//! The tree only prunes. A query reports every item whose box meets the
//! query region (and possibly more, where a caller pads a box); callers then
//! decide each reported item with their own exact test, so using the tree
//! never changes an answer, only how many items are looked at.

use glam::DVec3;

/// Items per leaf, at most.
const LEAF: usize = 4;

#[derive(Clone, Copy, Debug)]
struct Node {
    lo: DVec3,
    hi: DVec3,
    /// A leaf: its first item in `order`. An inner node: the index of its
    /// right child (the left child is the next node).
    first: u32,
    /// The items in a leaf; zero for an inner node.
    count: u32,
}

/// A bounding-volume hierarchy over item boxes (see the module docs).
pub(crate) struct Bvh {
    nodes: Vec<Node>,
    /// Item indices, grouped by leaf.
    order: Vec<u32>,
    /// Each item's box, by item index.
    boxes: Vec<(DVec3, DVec3)>,
}

impl Bvh {
    /// A tree over the items with the given boxes (`lo`, `hi`), indexed by
    /// their position in `boxes`.
    pub fn new(boxes: Vec<(DVec3, DVec3)>) -> Bvh {
        let n = boxes.len();
        assert!(n < u32::MAX as usize, "too many items for a tree");
        let mut order: Vec<u32> = (0..n as u32).collect();
        let centres: Vec<DVec3> = boxes.iter().map(|(lo, hi)| (*lo + *hi) * 0.5).collect();
        let mut nodes = Vec::with_capacity(2 * n.div_ceil(LEAF));
        if n > 0 {
            build(&mut nodes, &mut order, 0, n, &boxes, &centres);
        }
        Bvh {
            nodes,
            order,
            boxes,
        }
    }

    /// A tree over points, each its own (empty) box.
    pub fn of_points(points: &[DVec3]) -> Bvh {
        Bvh::new(points.iter().map(|p| (*p, *p)).collect())
    }

    /// The box around every item, if there are any.
    pub fn bounds(&self) -> Option<(DVec3, DVec3)> {
        self.nodes.first().map(|n| (n.lo, n.hi))
    }

    /// Walks the tree depth first: a node or an item is entered only if
    /// `test` accepts its box, and `item` is called with each accepted item.
    /// The walk stops as soon as `item` returns false; the result is false
    /// if it stopped early.
    pub fn visit(
        &self,
        mut test: impl FnMut(DVec3, DVec3) -> bool,
        mut item: impl FnMut(usize) -> bool,
    ) -> bool {
        if self.nodes.is_empty() {
            return true;
        }
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);
        while let Some(i) = stack.pop() {
            let n = &self.nodes[i as usize];
            if !test(n.lo, n.hi) {
                continue;
            }
            if n.count == 0 {
                // Right child pushed first, so the left is walked first.
                stack.push(n.first);
                stack.push(i + 1);
                continue;
            }
            for &k in &self.order[n.first as usize..(n.first + n.count) as usize] {
                let (lo, hi) = self.boxes[k as usize];
                if test(lo, hi) && !item(k as usize) {
                    return false;
                }
            }
        }
        true
    }

    /// The items whose boxes meet the closed box `lo`–`hi`, in increasing
    /// order, into `out` (cleared first).
    pub fn overlapping(&self, lo: DVec3, hi: DVec3, out: &mut Vec<usize>) {
        out.clear();
        self.visit(
            |a, b| a.cmple(hi).all() && lo.cmple(b).all(),
            |k| {
                out.push(k);
                true
            },
        );
        out.sort_unstable();
    }
}

/// Builds the subtree over `order[start..end]` and returns its root.
fn build(
    nodes: &mut Vec<Node>,
    order: &mut [u32],
    start: usize,
    end: usize,
    boxes: &[(DVec3, DVec3)],
    centres: &[DVec3],
) -> u32 {
    let (mut lo, mut hi) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
    let (mut clo, mut chi) = (lo, hi);
    for &k in &order[start..end] {
        let (a, b) = boxes[k as usize];
        lo = lo.min(a);
        hi = hi.max(b);
        clo = clo.min(centres[k as usize]);
        chi = chi.max(centres[k as usize]);
    }
    let i = nodes.len();
    nodes.push(Node {
        lo,
        hi,
        first: start as u32,
        count: (end - start) as u32,
    });
    if end - start <= LEAF {
        return i as u32;
    }
    let ext = chi - clo;
    let axis = if ext.x >= ext.y && ext.x >= ext.z {
        0
    } else if ext.y >= ext.z {
        1
    } else {
        2
    };
    let mid = start + (end - start) / 2;
    order[start..end].select_nth_unstable_by(mid - start, |&a, &b| {
        centres[a as usize][axis]
            .total_cmp(&centres[b as usize][axis])
            .then(a.cmp(&b))
    });
    build(nodes, order, start, mid, boxes, centres);
    let right = build(nodes, order, mid, end, boxes, centres);
    nodes[i].first = right;
    nodes[i].count = 0;
    i as u32
}

/// True if the segment `p`–`q` meets the box `lo`–`hi` grown by `pad` on
/// every side (the slab test, Ericson §5.3.3). A positive `pad` makes the
/// test conservative against the rounding of the slab arithmetic: a
/// segment that meets the box is never rejected.
pub(crate) fn segment_meets_box(p: DVec3, q: DVec3, lo: DVec3, hi: DVec3, pad: f64) -> bool {
    let (lo, hi) = (lo - DVec3::splat(pad), hi + DVec3::splat(pad));
    let d = q - p;
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for k in 0..3 {
        if d[k] == 0.0 {
            if p[k] < lo[k] || p[k] > hi[k] {
                return false;
            }
            continue;
        }
        let (a, b) = ((lo[k] - p[k]) / d[k], (hi[k] - p[k]) / d[k]);
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        t0 = t0.max(a);
        t1 = t1.min(b);
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// True if `p` lies in the box `lo`–`hi` grown by `r` on every side.
pub(crate) fn near_box(p: DVec3, lo: DVec3, hi: DVec3, r: f64) -> bool {
    (lo - DVec3::splat(r)).cmple(p).all() && p.cmple(hi + DVec3::splat(r)).all()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic pseudo-random sequence in [0, 1).
    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*state >> 11) as f64 / (1u64 << 53) as f64
    }

    fn random_boxes(n: usize, seed: u64) -> Vec<(DVec3, DVec3)> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                let c = DVec3::new(lcg(&mut s), lcg(&mut s), lcg(&mut s)) * 10.0;
                let h = DVec3::new(lcg(&mut s), lcg(&mut s), lcg(&mut s)) * 0.5;
                (c - h, c + h)
            })
            .collect()
    }

    #[test]
    fn box_queries_match_a_full_scan() {
        let boxes = random_boxes(500, 1);
        let t = Bvh::new(boxes.clone());
        let mut s = 7;
        let mut out = Vec::new();
        for _ in 0..200 {
            let c = DVec3::new(lcg(&mut s), lcg(&mut s), lcg(&mut s)) * 10.0;
            let h = DVec3::splat(lcg(&mut s));
            let (lo, hi) = (c - h, c + h);
            t.overlapping(lo, hi, &mut out);
            let want: Vec<usize> = (0..boxes.len())
                .filter(|&k| boxes[k].0.cmple(hi).all() && lo.cmple(boxes[k].1).all())
                .collect();
            assert_eq!(out, want);
        }
    }

    #[test]
    fn segment_queries_match_a_full_scan() {
        let boxes = random_boxes(300, 2);
        let t = Bvh::new(boxes.clone());
        let mut s = 9;
        for _ in 0..200 {
            let p = DVec3::new(lcg(&mut s), lcg(&mut s), lcg(&mut s)) * 12.0 - 1.0;
            let q = DVec3::new(lcg(&mut s), lcg(&mut s), lcg(&mut s)) * 12.0 - 1.0;
            let mut got = Vec::new();
            t.visit(
                |lo, hi| segment_meets_box(p, q, lo, hi, 0.0),
                |k| {
                    got.push(k);
                    true
                },
            );
            got.sort_unstable();
            let want: Vec<usize> = (0..boxes.len())
                .filter(|&k| segment_meets_box(p, q, boxes[k].0, boxes[k].1, 0.0))
                .collect();
            assert_eq!(got, want);
        }
    }

    #[test]
    fn slab_test_cases() {
        let (lo, hi) = (DVec3::ZERO, DVec3::ONE);
        let x = DVec3::X;
        // Through, ending short, parallel outside and inside, touching.
        assert!(segment_meets_box(
            -x,
            x * 2.0 + DVec3::splat(0.5),
            lo,
            hi,
            0.0
        ));
        assert!(!segment_meets_box(
            DVec3::new(-2.0, 0.5, 0.5),
            DVec3::new(-1.0, 0.5, 0.5),
            lo,
            hi,
            0.0
        ));
        assert!(!segment_meets_box(
            DVec3::new(-1.0, 2.0, 0.5),
            DVec3::new(2.0, 2.0, 0.5),
            lo,
            hi,
            0.0
        ));
        assert!(segment_meets_box(
            DVec3::new(-1.0, 1.0, 0.5),
            DVec3::new(2.0, 1.0, 0.5),
            lo,
            hi,
            0.0
        ));
        assert!(!segment_meets_box(
            DVec3::new(-1.0, 1.0 + 1e-9, 0.5),
            DVec3::new(2.0, 1.0 + 1e-9, 0.5),
            lo,
            hi,
            0.0
        ));
        assert!(segment_meets_box(
            DVec3::new(-1.0, 1.0 + 1e-9, 0.5),
            DVec3::new(2.0, 1.0 + 1e-9, 0.5),
            lo,
            hi,
            1e-8
        ));
        assert!(near_box(DVec3::new(1.5, 0.5, 0.5), lo, hi, 0.5));
        assert!(!near_box(DVec3::new(1.6, 0.5, 0.5), lo, hi, 0.5));
    }

    #[test]
    fn equal_boxes_and_empty_trees() {
        let t = Bvh::new(vec![(DVec3::ZERO, DVec3::ONE); 100]);
        let mut out = Vec::new();
        t.overlapping(DVec3::splat(0.5), DVec3::splat(0.5), &mut out);
        assert_eq!(out, (0..100).collect::<Vec<_>>());
        let e = Bvh::new(Vec::new());
        assert_eq!(e.bounds(), None);
        e.overlapping(DVec3::ZERO, DVec3::ONE, &mut out);
        assert!(out.is_empty());
        let p = Bvh::of_points(&[DVec3::X, DVec3::Y]);
        p.overlapping(DVec3::ZERO, DVec3::new(1.0, 0.5, 0.0), &mut out);
        assert_eq!(out, vec![0]);
    }

    #[test]
    fn early_stop() {
        let t = Bvh::new(random_boxes(100, 3));
        let mut seen = 0;
        let done = t.visit(
            |_, _| true,
            |_| {
                seen += 1;
                seen < 10
            },
        );
        assert!(!done);
        assert_eq!(seen, 10);
    }
}
