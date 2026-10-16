//! A flat face's loops in its plane, with trees over their sides and over
//! the loops' boxes built when a point first needs them: where a point lies
//! against the face is decided by the sides and loops near it, not by a
//! walk around every loop (`assemble::check_apart`, the strip search of
//! `imprint::held_strip`, design note 0012 §20, and the placing of points
//! in long cycles of `planar::regions`).

use std::cell::OnceCell;

use glam::{DVec2, DVec3};

use super::shape::{Face, Shape};
use super::{slot, Pair};
use crate::bvh::Bvh;
use crate::origin::Operand;

/// The loops of a flat face in its plane: their sides, loop after loop,
/// with trees over them built when first needed.
pub(super) struct FlatOutline {
    sides: Vec<(DVec2, DVec2)>,
    /// Where each loop's sides start in `sides`, then their end.
    starts: Vec<usize>,
    /// Trees of every side's box and of every loop's box.
    trees: OnceCell<(Bvh, Bvh)>,
    /// A tree of each long loop's sides.
    long: Vec<OnceCell<Bvh>>,
}

/// Loops with more sides than this are searched through a tree of them.
pub(super) const LONG_LOOP: usize = 32;

fn lift(v: DVec2) -> DVec3 {
    DVec3::new(v.x, v.y, 0.0)
}

fn side_boxes(sides: &[(DVec2, DVec2)]) -> Vec<(DVec3, DVec3)> {
    (sides.iter())
        .map(|&(a, b)| (lift(a.min(b)), lift(a.max(b))))
        .collect()
}

impl FlatOutline {
    /// The outline of the loops whose sides are `sides`, those of loop `k`
    /// from `starts[k]` to `starts[k + 1]`.
    pub(super) fn new(sides: Vec<(DVec2, DVec2)>, starts: Vec<usize>) -> FlatOutline {
        FlatOutline {
            long: (1..starts.len()).map(|_| OnceCell::new()).collect(),
            sides,
            starts,
            trees: OnceCell::new(),
        }
    }

    /// The outline of face `f` of `shape`, in the face's frame
    /// ([`Face::to2`]): its loops' straight sides from corner to corner (a
    /// curved edge by its chord).
    pub(super) fn of_face(shape: &Shape, f: &Face) -> FlatOutline {
        let mut sides = Vec::with_capacity(f.coedges());
        let mut starts = vec![0];
        for vs in &f.verts {
            let q: Vec<DVec2> = vs.iter().map(|&v| f.to2(shape.points[v])).collect();
            sides.extend((0..q.len()).map(|k| (q[k], q[(k + 1) % q.len()])));
            starts.push(sides.len());
        }
        FlatOutline::new(sides, starts)
    }

    /// How many sides its loops have.
    pub(super) fn len(&self) -> usize {
        self.sides.len()
    }

    fn trees(&self) -> &(Bvh, Bvh) {
        self.trees.get_or_init(|| {
            let loops = (self.starts.windows(2))
                .map(|w| {
                    (side_boxes(&self.sides[w[0]..w[1]]).into_iter()).fold(
                        (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
                        |(lo, hi), (a, b)| (lo.min(a), hi.max(b)),
                    )
                })
                .collect();
            (Bvh::new(side_boxes(&self.sides)), Bvh::new(loops))
        })
    }

    /// Whether `q` is within `eps` of a side, adding to `seen` the sides
    /// measured.
    pub(super) fn near(&self, q: DVec2, eps: f64, seen: &mut usize) -> bool {
        let r = DVec3::new(eps, eps, 0.0);
        let (lo, hi) = (lift(q) - r, lift(q) + r);
        !self.trees().0.visit(
            |a, b| a.cmple(hi).all() && lo.cmple(b).all(),
            |k| {
                *seen += 1;
                let (a, b) = self.sides[k];
                let d = b - a;
                let t =
                    ((q - a).dot(d) / d.length_squared().max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
                (q - (a + d * t)).length() > eps
            },
        )
    }

    /// Whether `q` is inside the loops: whether the ray from it toward +x
    /// crosses an odd number of sides, adding to `seen` the loops and the
    /// sides looked at. A loop whose box does not hold `q` is crossed an
    /// even number of times and is skipped, so a point among many islands
    /// or holes meets only the loops around it (found by review: a plate
    /// with 1,600 flush pins, whose rays crossed every hole in their row).
    pub(super) fn inside(&self, q: DVec2, seen: &mut usize) -> bool {
        let crosses = |&(a, b): &(DVec2, DVec2)| {
            (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x)
        };
        let mut inside = false;
        let mut held = Vec::new();
        self.trees().1.overlapping(lift(q), lift(q), &mut held);
        *seen += held.len();
        for k in held {
            let sides = &self.sides[self.starts[k]..self.starts[k + 1]];
            if sides.len() <= LONG_LOOP {
                *seen += sides.len();
                inside ^= sides.iter().filter(|s| crosses(s)).count() % 2 == 1;
                continue;
            }
            let tree = self.long[k].get_or_init(|| Bvh::new(side_boxes(sides)));
            tree.visit(
                |lo, hi| lo.y <= q.y && q.y <= hi.y && q.x <= hi.x,
                |i| {
                    *seen += 1;
                    inside ^= crosses(&sides[i]);
                    true
                },
            );
        }
        inside
    }
}

/// The outlines of both operands' faces, each built when first needed (and
/// the table of them when the first is), so a boolean that needs none
/// spends nothing on them.
pub(super) struct FlatOutlines(OnceCell<[Vec<OnceCell<FlatOutline>>; 2]>);

impl FlatOutlines {
    /// No outline built yet.
    pub(super) fn new() -> FlatOutlines {
        FlatOutlines(OnceCell::new())
    }

    /// The outline of face `f` of `op` (flat), one of `pair`'s faces, and
    /// whether it was built by this call.
    pub(super) fn get(&self, pair: &Pair, op: Operand, f: usize) -> (&FlatOutline, bool) {
        let cells = |s: &Shape| (0..s.faces.len()).map(|_| OnceCell::new()).collect();
        let table = self
            .0
            .get_or_init(|| [cells(&pair.first), cells(&pair.second)]);
        let cell = &table[slot(op)][f];
        let fresh = cell.get().is_none();
        let shape = pair.get(op);
        (
            cell.get_or_init(|| FlatOutline::of_face(shape, &shape.faces[f])),
            fresh,
        )
    }
}
