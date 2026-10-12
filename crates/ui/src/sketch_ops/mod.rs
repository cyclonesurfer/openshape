//! Sketch editing operations (design note 0007): pure functions from a
//! sketch and the user's input to one [`SketchEdit`] batch. The editor in
//! [`crate::sketch_mode`] decides *what* the user asked for; these functions
//! compute the exact geometry and the constraints that keep it.

pub mod fillet;
pub mod mirror;
pub mod offset;
pub mod shapes;
pub mod trim;

use std::collections::{BTreeMap, BTreeSet};

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, ConstraintEntry, Entity, Geometry, TangentKind};
use openshape_core_model::{SketchData, SketchEdit, SketchId};
use openshape_sketch_solver::ResolvedCurve;

use crate::sketch_mode::{Snap, SnapTarget};

/// Positions closer than this are the same point (metres).
pub(crate) const SAME_POINT: f64 = 1e-7;

/// Builds one sketch edit batch with batch-local ids.
#[derive(Default)]
pub struct Batch {
    /// The edits so far.
    pub edits: Vec<SketchEdit>,
    next: u32,
}

impl Batch {
    /// Adds an entity.
    pub fn entity(&mut self, geometry: Geometry<SketchId>, construction: bool) -> SketchId {
        let id = SketchId::new_in_batch(self.next);
        self.next += 1;
        self.edits.push(SketchEdit::AddEntity(Entity {
            id,
            geometry,
            construction,
        }));
        id
    }

    /// Adds a point.
    pub fn point(&mut self, at: DVec2) -> SketchId {
        self.entity(Geometry::Point { at }, false)
    }

    /// Adds a line with two new end points.
    pub fn line(&mut self, a: DVec2, b: DVec2, construction: bool) -> Ends {
        let start = self.point(a);
        let end = self.point(b);
        let curve = self.entity(Geometry::Line { start, end }, construction);
        Ends { curve, start, end }
    }

    /// Adds a counter-clockwise arc about `center` from `a` to `b` with new
    /// points (`a` and `b` should be equally far from `center`).
    pub fn arc(&mut self, center: DVec2, a: DVec2, b: DVec2, construction: bool) -> Arc {
        let c = self.point(center);
        let start = self.point(a);
        let end = self.point(b);
        let curve = self.entity(
            Geometry::Arc {
                center: c,
                start,
                end,
            },
            construction,
        );
        Arc {
            curve,
            center: c,
            start,
            end,
        }
    }

    /// Adds a constraint.
    pub fn constraint(&mut self, constraint: Constraint<SketchId>) -> SketchId {
        let id = SketchId::new_in_batch(self.next);
        self.next += 1;
        self.edits.push(SketchEdit::AddConstraint(ConstraintEntry {
            id,
            constraint,
        }));
        id
    }

    /// Makes two points coincident.
    pub fn coincident(&mut self, a: SketchId, b: SketchId) -> SketchId {
        self.constraint(Constraint::Coincident { a, b })
    }

    /// Makes two curves tangent where they meet at `point` (an end point of
    /// both, or coincident with one of each).
    pub fn tangent_at(&mut self, a: SketchId, b: SketchId, point: SketchId) -> SketchId {
        self.constraint(Constraint::Tangent {
            a,
            b,
            kind: TangentKind::AtPoint { point },
        })
    }

    /// A new point at a snap, constrained to what it snapped to.
    pub fn snapped_point(&mut self, s: &Snap) -> SketchId {
        let p = self.point(s.at);
        self.attach(p, s);
        p
    }

    /// Constrains an existing new point to what `s` snapped to.
    pub fn attach(&mut self, p: SketchId, s: &Snap) {
        match s.target {
            SnapTarget::Free => {}
            SnapTarget::Point(q) => {
                self.coincident(q, p);
            }
            SnapTarget::Curve(c) => {
                self.constraint(Constraint::PointOnCurve { point: p, curve: c });
            }
            SnapTarget::Origin => {
                self.constraint(Constraint::Fix {
                    point: p,
                    at: DVec2::ZERO,
                });
            }
        }
    }
}

/// A line (or another curve with two ends) added to a [`Batch`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ends {
    /// The curve.
    pub curve: SketchId,
    /// Its start point.
    pub start: SketchId,
    /// Its end point.
    pub end: SketchId,
}

/// An arc added to a [`Batch`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Arc {
    /// The arc.
    pub curve: SketchId,
    /// Its centre point.
    pub center: SketchId,
    /// Its start point.
    pub start: SketchId,
    /// Its end point.
    pub end: SketchId,
}

/// The items to delete for deleting `selected`: the selection plus the
/// points that only selected curves use (the delete edit removes
/// constraints on any of them).
pub fn deletion(sketch: &SketchData, selected: &BTreeSet<SketchId>) -> Vec<SketchId> {
    let kept_points: BTreeSet<SketchId> = sketch
        .entities
        .iter()
        .filter(|e| !selected.contains(&e.id))
        .flat_map(|e| e.geometry.referenced_points())
        .collect();
    let mut items = selected.clone();
    for e in &sketch.entities {
        if selected.contains(&e.id) {
            items.extend(
                e.geometry
                    .referenced_points()
                    .into_iter()
                    .filter(|p| !kept_points.contains(p)),
            );
        }
    }
    items.into_iter().collect()
}

/// The left normal of `v` (turned a quarter counter-clockwise).
pub fn perp(v: DVec2) -> DVec2 {
    DVec2::new(-v.y, v.x)
}

/// The curve that has `point` as its start or end point, with `true` when
/// it is the end: the curve a new tangent curve continues.
pub fn curve_ending_at(sketch: &SketchData, point: SketchId) -> Option<(SketchId, bool)> {
    sketch.entities.iter().find_map(|e| {
        let (first, last) = match &e.geometry {
            Geometry::Line { start, end } | Geometry::Arc { start, end, .. } => (*start, *end),
            Geometry::ControlSpline { control: p, .. }
            | Geometry::InterpolatedSpline { through: p, .. } => (*p.first()?, *p.last()?),
            _ => return None,
        };
        if last == point {
            Some((e.id, true))
        } else if first == point {
            Some((e.id, false))
        } else {
            None
        }
    })
}

/// The unit direction in which a curve leaves its end (`at_end`) or start
/// (`!at_end`): the direction a tangent continuation starts in.
pub fn leaving_direction(curve: &ResolvedCurve, at_end: bool) -> Option<DVec2> {
    let d = match curve {
        ResolvedCurve::Line { start, end } => {
            if at_end {
                *end - *start
            } else {
                *start - *end
            }
        }
        ResolvedCurve::Arc {
            start_angle,
            sweep_angle,
            ..
        } => {
            // Counter-clockwise tangent at an angle a is (−sin a, cos a).
            if at_end {
                perp(DVec2::from_angle(start_angle + sweep_angle))
            } else {
                -perp(DVec2::from_angle(*start_angle))
            }
        }
        ResolvedCurve::Spline { points, .. } => {
            // A clamped B-spline is tangent to its first and last legs.
            let n = points.len();
            if n < 2 {
                return None;
            }
            if at_end {
                points[n - 1] - points[n - 2]
            } else {
                points[0] - points[1]
            }
        }
        ResolvedCurve::Circle { .. } => return None,
    };
    let len = d.length();
    (len.is_finite() && len > 0.0).then(|| d / len)
}

/// The line or circle a curve lies on: a line through `p` in the unit
/// direction `d`, or a circle.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Carrier {
    /// A line through `p` in the unit direction `d`.
    Line {
        /// A point of the line.
        p: DVec2,
        /// Its unit direction.
        d: DVec2,
    },
    /// A circle.
    Circle {
        /// Centre.
        c: DVec2,
        /// Radius.
        r: f64,
    },
}

/// Where two carriers meet (none for parallel lines, concentric circles).
pub(crate) fn meet(a: Carrier, b: Carrier) -> Vec<DVec2> {
    match (a, b) {
        (Carrier::Line { p, d }, Carrier::Line { p: q, d: e }) => {
            let den = d.perp_dot(e);
            if den.abs() < 1e-12 {
                return Vec::new();
            }
            vec![p + d * ((q - p).perp_dot(e) / den)]
        }
        (Carrier::Line { p, d }, Carrier::Circle { c, r })
        | (Carrier::Circle { c, r }, Carrier::Line { p, d }) => {
            // |p + d s − c|² = r².
            let w = p - c;
            let (b, k) = (w.dot(d), w.length_squared() - r * r);
            let disc = b * b - k;
            if disc < 0.0 {
                return Vec::new();
            }
            let h = disc.sqrt();
            vec![p + d * (-b - h), p + d * (-b + h)]
        }
        (Carrier::Circle { c, r }, Carrier::Circle { c: k, r: s }) => {
            let v = k - c;
            let dist = v.length();
            if dist < 1e-12 || dist > r + s || dist < (r - s).abs() {
                return Vec::new();
            }
            let a = (r * r - s * s + dist * dist) / (2.0 * dist);
            let h = (r * r - a * a).max(0.0).sqrt();
            let m = c + v * (a / dist);
            let n = perp(v / dist);
            vec![m + n * h, m - n * h]
        }
    }
}

/// A small union-find over sketch ids.
#[derive(Default)]
pub(crate) struct UnionFind(BTreeMap<SketchId, SketchId>);

impl UnionFind {
    pub(crate) fn find(&mut self, x: SketchId) -> SketchId {
        let mut root = x;
        while let Some(&p) = self.0.get(&root) {
            if p == root {
                break;
            }
            root = p;
        }
        self.0.insert(x, root);
        root
    }

    pub(crate) fn union(&mut self, a: SketchId, b: SketchId) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0.insert(ra, rb);
        }
    }
}

/// Each point's coincidence class (points in no coincidence are absent).
pub(crate) fn coincidence_classes(sketch: &SketchData) -> BTreeMap<SketchId, usize> {
    let mut uf = UnionFind::default();
    for c in &sketch.constraints {
        if let Constraint::Coincident { a, b } = c.constraint {
            uf.union(a, b);
        }
    }
    let points: Vec<SketchId> = uf.0.keys().copied().collect();
    let mut roots: BTreeMap<SketchId, usize> = BTreeMap::new();
    let mut out = BTreeMap::new();
    for p in points {
        let r = uf.find(p);
        let n = roots.len();
        let class = *roots.entry(r).or_insert(n);
        out.insert(p, class);
    }
    out
}
