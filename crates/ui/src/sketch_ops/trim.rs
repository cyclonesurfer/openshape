//! Trim, extend and split (design note 0007 §7): cut a curve where the
//! other curves of the sketch meet it, and keep every constraint that still
//! means the same thing.
//!
//! The planners are pure: they read the sketch and its solved curves and
//! return one edit batch ([`Plan`]) for the editor to commit. Where a cut
//! end is attached to the curve that cut it, a point-on-curve constraint is
//! used rather than a coincidence: at a T junction the cutter's end often
//! already lies on the cut curve, and a coincidence would repeat that.

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, Geometry, TangentKind};
use openshape_core_model::{SketchData, SketchEdit, SketchId};
use openshape_features::kernel_curve;
use openshape_kernel::{Kernel, ProfileCurve};
use openshape_sketch_solver::ResolvedCurve;

use super::{deletion, Batch};

/// Cuts closer than this to one another, or to the cut curve's ends, are
/// one cut (metres): a piece a micron long is not a piece anyone meant.
pub const CUT_GAP: f64 = 1e-6;

/// Curves meeting at an angle whose sine is below this (about 0.6°)
/// touch rather than cross, see [`Cut::cutters`].
pub const TANGENT_SIN: f64 = 1e-2;

/// Where other curves cut a curve.
#[derive(Clone, Debug, PartialEq)]
pub struct Cut {
    /// Parameter along the cut curve (the arrangement's conventions: a
    /// line's fraction, the angle travelled along an arc from its start, a
    /// circle's angle from +x).
    pub t: f64,
    /// The point.
    pub point: DVec2,
    /// The curves crossing the cut curve there, which its new end can be
    /// attached to. Curves touching it there (meeting it at less than
    /// [`TANGENT_SIN`], running along it on the same line or circle) still
    /// cut it but are left out: there a point-on-curve constraint is
    /// ill-conditioned, or repeats what a tangency already says.
    pub cutters: Vec<SketchId>,
}

/// The cuts of one curve, sorted along it.
#[derive(Clone, Debug, PartialEq)]
pub struct Cuts {
    /// The cut curve.
    pub target: SketchId,
    /// Its solved geometry.
    pub curve: ResolvedCurve,
    /// The cuts, away from the curve's ends.
    pub cuts: Vec<Cut>,
}

/// One trim, extend or split, ready to commit.
#[derive(Clone, Debug, PartialEq)]
pub struct Plan {
    /// The edit batch.
    pub edits: Vec<SketchEdit>,
    /// Constraints deleted because what they constrained no longer exists
    /// (an end that moved, a segment's length).
    pub removed: usize,
    /// What to show while hovering.
    pub preview: Preview,
}

/// What a plan shows before it is committed.
#[derive(Clone, Debug, PartialEq)]
pub enum Preview {
    /// The piece that goes (trim) or comes (extend).
    Piece(ResolvedCurve),
    /// Where the curve is cut (split).
    Point(DVec2),
}

/// Finds where the other curves cut `target` (`curves` are the sketch's
/// solved curves, construction geometry included).
pub fn cuts(
    kernel: &(impl Kernel + ?Sized),
    curves: &[(SketchId, ResolvedCurve)],
    target: SketchId,
) -> Result<Cuts, String> {
    let curve = curves
        .iter()
        .find(|(id, _)| *id == target)
        .map(|(_, c)| c.clone())
        .ok_or("the curve no longer exists")?;
    let a = kernel_curve(&curve);
    let ends = end_points(&curve);
    let mut raw = Vec::new();
    for (id, c) in curves {
        if *id == target {
            continue;
        }
        for x in contacts(kernel, &a, c)? {
            let at_end = ends.is_some_and(|(s, e)| {
                (x.point - s).length() <= CUT_GAP || (x.point - e).length() <= CUT_GAP
            });
            if !at_end {
                raw.push(x.cut(*id));
            }
        }
    }
    raw.sort_by(|x, y| x.t.total_cmp(&y.t));
    let mut out: Vec<Cut> = Vec::with_capacity(raw.len());
    for c in raw {
        match out.last_mut() {
            Some(l) if (l.point - c.point).length() <= CUT_GAP => merge(l, c),
            _ => out.push(c),
        }
    }
    if matches!(curve, ResolvedCurve::Circle { .. }) && out.len() >= 2 {
        // The first and last cut of a circle can meet across its seam.
        if (out[0].point - out[out.len() - 1].point).length() <= CUT_GAP {
            let last = out.pop().expect("two cuts");
            merge(&mut out[0], last);
        }
    }
    Ok(Cuts {
        target,
        curve,
        cuts: out,
    })
}

fn merge(into: &mut Cut, c: Cut) {
    for id in c.cutters {
        if !into.cutters.contains(&id) {
            into.cutters.push(id);
        }
    }
}

/// Where curve `a` meets another curve.
struct Contact {
    /// Parameter along `a`.
    t: f64,
    point: DVec2,
    /// The curves cross there rather than touch (see [`Cut::cutters`]).
    crossing: bool,
}

impl Contact {
    /// The cut it makes, by curve `id`.
    fn cut(&self, id: SketchId) -> Cut {
        Cut {
            t: self.t,
            point: self.point,
            cutters: if self.crossing { vec![id] } else { Vec::new() },
        }
    }
}

/// The kernel's contacts of `a` with a sketch curve; a degenerate sketch
/// curve (a zero-length line, say) meets nothing.
fn contacts(
    kernel: &(impl Kernel + ?Sized),
    a: &ProfileCurve,
    b: &ResolvedCurve,
) -> Result<Vec<Contact>, String> {
    let b = kernel_curve(b);
    if kernel.curve_point(&b, 0.0).is_err() {
        return Ok(Vec::new());
    }
    let found = kernel
        .curve_contacts(a, &b)
        .map_err(|e| format!("cannot intersect the curves: {e}"))?;
    Ok(found
        .into_iter()
        .map(|x| {
            let crossing = match (
                direction(kernel, a, x.ta, x.point),
                direction(kernel, &b, x.tb, x.point),
            ) {
                (Some(u), Some(v)) => u.perp_dot(v).abs() > TANGENT_SIN,
                _ => false,
            };
            Contact {
                t: x.ta,
                point: x.point,
                crossing,
            }
        })
        .collect())
}

/// The unit direction of a curve at parameter `t` (point `p`), when it
/// has one.
fn direction(kernel: &(impl Kernel + ?Sized), c: &ProfileCurve, t: f64, p: DVec2) -> Option<DVec2> {
    let d = match c {
        ProfileCurve::Line { start, end } => *end - *start,
        ProfileCurve::Arc { center, .. } | ProfileCurve::Circle { center, .. } => {
            super::perp(p - *center)
        }
        ProfileCurve::Spline { knots, .. } => {
            // A central difference over a millionth of the domain.
            let h = 1e-6 * (knots.last()? - knots.first()?);
            kernel.curve_point(c, t + h).ok()? - kernel.curve_point(c, t - h).ok()?
        }
    };
    let n = d.length();
    (n.is_finite() && n > 0.0).then(|| d / n)
}

/// Plans trimming the piece of `cuts.target` nearest `near` (design note
/// 0007 §7.2).
pub fn trim(sketch: &SketchData, cuts: &Cuts, near: DVec2) -> Result<Plan, String> {
    let t = cuts.target;
    let geometry = &sketch
        .entity(t)
        .ok_or("the curve no longer exists")?
        .geometry;
    let k = cuts.cuts.len();
    if k == 0 {
        // Nothing cuts it: the whole curve goes.
        let items = deletion(sketch, &BTreeSet::from([t]));
        let removed = sketch
            .constraints
            .iter()
            .filter(|c| c.constraint.referenced().iter().any(|e| items.contains(e)))
            .count();
        return Ok(Plan {
            edits: vec![SketchEdit::Delete { items }],
            removed,
            preview: Preview::Piece(cuts.curve.clone()),
        });
    }
    match (geometry, &cuts.curve) {
        (Geometry::Circle { center, .. }, ResolvedCurve::Circle { .. }) => {
            trim_circle(sketch, cuts, *center, near)
        }
        (Geometry::Line { start, end }, _) | (Geometry::Arc { start, end, .. }, _) => {
            let center = match geometry {
                Geometry::Arc { center, .. } => Some(*center),
                _ => None,
            };
            let tn = param(&cuts.curve, near);
            let i = cuts.cuts.iter().filter(|c| c.t < tn).count();
            let mut r = Rewrite::new(sketch, t);
            r.drop_segment_measures();
            if i == 0 || i == k {
                // An end piece: that end moves to the cut next to it.
                let (p, cut) = if i == 0 {
                    (*start, &cuts.cuts[0])
                } else {
                    (*end, &cuts.cuts[k - 1])
                };
                r.detach(p)?;
                r.moves.push((p, cut.point));
                r.attach(p, cut);
                let from = sketch.point(p).ok_or("missing end point")?;
                let preview = piece(&cuts.curve, param(&cuts.curve, from), cut.t);
                Ok(r.finish(Preview::Piece(preview)))
            } else {
                // A middle piece: the curve keeps the part before it; a new
                // curve takes the part after it and the old end.
                let (lo, hi) = (&cuts.cuts[i - 1], &cuts.cuts[i]);
                let d = r.continuation(*end, center, hi.point)?;
                r.moves.push((*end, lo.point));
                r.attach(*end, lo);
                r.attach(d.start, hi);
                match center {
                    None => {
                        r.pending.push(Constraint::PointOnCurve {
                            point: d.start,
                            curve: t,
                        });
                        r.pending.push(Constraint::PointOnCurve {
                            point: d.end,
                            curve: t,
                        });
                    }
                    Some(_) => r.pending.push(Constraint::Equal { a: t, b: d.curve }),
                }
                Ok(r.finish(Preview::Piece(piece(&cuts.curve, lo.t, hi.t))))
            }
        }
        (Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. }, _) => Err(
            "A spline can only be trimmed whole: its pieces are not exact splines of the \
             sketch's kinds (design note 0007 §7.5). Delete it, or redraw it in pieces"
                .into(),
        ),
        _ => Err("this curve cannot be trimmed".into()),
    }
}

/// A circle cut at two or more places becomes an arc on the same centre
/// point with every constraint of the circle (they all use the full
/// circle).
fn trim_circle(
    sketch: &SketchData,
    cuts: &Cuts,
    center: SketchId,
    near: DVec2,
) -> Result<Plan, String> {
    let k = cuts.cuts.len();
    if k < 2 {
        return Err(
            "A circle needs two crossings to trim: one crossing leaves nothing to cut away".into(),
        );
    }
    let t = cuts.target;
    let tn = param(&cuts.curve, near);
    // The removed piece runs counter-clockwise from cut `lo` to cut `hi`.
    let i = cuts.cuts.iter().filter(|c| c.t <= tn).count();
    let (lo, hi) = if i == 0 || i == k {
        (&cuts.cuts[k - 1], &cuts.cuts[0])
    } else {
        (&cuts.cuts[i - 1], &cuts.cuts[i])
    };
    let mut r = Rewrite::new(sketch, t);
    let construction = sketch.entity(t).is_some_and(|e| e.construction);
    let start = r.batch.point(hi.point);
    let end = r.batch.point(lo.point);
    let arc = r
        .batch
        .entity(Geometry::Arc { center, start, end }, construction);
    for c in &sketch.constraints {
        if c.constraint.referenced().contains(&t) {
            let moved = c.constraint.map_keys(|k| if k == t { arc } else { k });
            r.replace(c.id, moved);
        }
    }
    r.delete.insert(t);
    r.attach_for(start, hi, t);
    r.attach_for(end, lo, t);
    Ok(r.finish(Preview::Piece(piece(&cuts.curve, lo.t, hi.t))))
}

/// Plans extending the end of line or arc `target` nearest `near` to the
/// first curve beyond it (design note 0007 §7.3).
pub fn extend(
    kernel: &(impl Kernel + ?Sized),
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    target: SketchId,
    near: DVec2,
) -> Result<Plan, String> {
    let curve = curves
        .iter()
        .find(|(id, _)| *id == target)
        .map(|(_, c)| c.clone())
        .ok_or("the curve no longer exists")?;
    let (start, end) = match &sketch
        .entity(target)
        .ok_or("the curve no longer exists")?
        .geometry
    {
        Geometry::Line { start, end } | Geometry::Arc { start, end, .. } => (*start, *end),
        _ => return Err("Only lines and arcs can be extended".into()),
    };
    let (s, e) = end_points(&curve).ok_or("Only lines and arcs can be extended")?;
    let at_end = (near - e).length() <= (near - s).length();
    let (p, from, other) = if at_end { (end, e, s) } else { (start, s, e) };
    // The way the end can go: a long ray, or the rest of the arc's circle.
    let (proxy, length) = match &curve {
        ResolvedCurve::Line { .. } => {
            let d = (from - other).normalize();
            let reach = 4.0 * extent(curves).max(from.length()).max(1e-3);
            (
                ProfileCurve::Line {
                    start: from,
                    end: from + d * reach,
                },
                reach,
            )
        }
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let rest = TAU - sweep_angle;
            if rest * radius <= 2.0 * CUT_GAP {
                return Err("The arc is already a full circle".into());
            }
            let proxy = if at_end {
                ProfileCurve::Arc {
                    center: *center,
                    radius: *radius,
                    start_angle: start_angle + sweep_angle,
                    sweep_angle: rest,
                }
            } else {
                ProfileCurve::Arc {
                    center: *center,
                    radius: *radius,
                    start_angle: *start_angle,
                    sweep_angle: -rest,
                }
            };
            (proxy, *radius)
        }
        _ => return Err("Only lines and arcs can be extended".into()),
    };
    // The nearest contact beyond the end (not the arc's own other end).
    let mut best: Option<(f64, Cut)> = None;
    for (id, c) in curves {
        if *id == target {
            continue;
        }
        for x in contacts(kernel, &proxy, c)? {
            let dist = x.t * length;
            if dist <= CUT_GAP || (x.point - other).length() <= CUT_GAP {
                continue;
            }
            match &mut best {
                Some((_, b)) if (b.point - x.point).length() <= CUT_GAP => merge(b, x.cut(*id)),
                Some((d, _)) if *d <= dist => {}
                _ => best = Some((dist, x.cut(*id))),
            }
        }
    }
    let (dist, cut) = best.ok_or("Nothing to extend to: no curve lies beyond this end")?;
    let mut r = Rewrite::new(sketch, target);
    r.drop_segment_measures();
    r.detach(p)?;
    r.moves.push((p, cut.point));
    r.attach(p, &cut);
    let preview = match &curve {
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let turn = dist / radius;
            ResolvedCurve::Arc {
                center: *center,
                radius: *radius,
                start_angle: if at_end {
                    start_angle + sweep_angle
                } else {
                    start_angle - turn
                },
                sweep_angle: turn,
            }
        }
        _ => ResolvedCurve::Line {
            start: from,
            end: cut.point,
        },
    };
    Ok(r.finish(Preview::Piece(preview)))
}

/// Plans splitting line or arc `target` at its point nearest `near` into
/// two curves joined there (design note 0007 §7.4).
pub fn split(
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    target: SketchId,
    near: DVec2,
) -> Result<Plan, String> {
    let curve = curves
        .iter()
        .find(|(id, _)| *id == target)
        .map(|(_, c)| c.clone())
        .ok_or("the curve no longer exists")?;
    let (end, center) = match &sketch
        .entity(target)
        .ok_or("the curve no longer exists")?
        .geometry
    {
        Geometry::Line { end, .. } => (*end, None),
        Geometry::Arc { end, center, .. } => (*end, Some(*center)),
        Geometry::Circle { .. } => {
            return Err("A circle has no ends to split between; trim it instead".into())
        }
        _ => return Err("Only lines and arcs can be split".into()),
    };
    let (s, e) = end_points(&curve).ok_or("Only lines and arcs can be split")?;
    let q = point_at(&curve, param(&curve, near));
    if (q - s).length() <= CUT_GAP || (q - e).length() <= CUT_GAP {
        return Err("Too close to an end to split".into());
    }
    let mut r = Rewrite::new(sketch, target);
    r.drop_segment_measures();
    let d = r.continuation(end, center, q)?;
    r.moves.push((end, q));
    r.pending
        .push(Constraint::Coincident { a: end, b: d.start });
    if center.is_none() {
        r.pending.push(Constraint::PointOnCurve {
            point: d.end,
            curve: target,
        });
    }
    Ok(r.finish(Preview::Point(q)))
}

/// A curve added after the cut, see [`Rewrite::continuation`].
struct Continuation {
    curve: SketchId,
    start: SketchId,
    end: SketchId,
}

/// Builds one plan: deletions, point moves, new entities, then constraints.
struct Rewrite<'a> {
    sketch: &'a SketchData,
    /// The curve being cut.
    target: SketchId,
    /// Coincidence classes of the sketch before the edit.
    classes: BTreeMap<SketchId, usize>,
    delete: BTreeSet<SketchId>,
    moves: Vec<(SketchId, DVec2)>,
    batch: Batch,
    /// Constraints to add once every new entity exists.
    pending: Vec<Constraint<SketchId>>,
    removed: usize,
}

impl<'a> Rewrite<'a> {
    fn new(sketch: &'a SketchData, target: SketchId) -> Rewrite<'a> {
        Rewrite {
            sketch,
            target,
            classes: coincidence_classes(sketch),
            delete: BTreeSet::new(),
            moves: Vec::new(),
            batch: Batch::default(),
            pending: Vec::new(),
            removed: 0,
        }
    }

    fn live(
        &self,
    ) -> impl Iterator<Item = &'a openshape_core_model::sketch::ConstraintEntry<SketchId>> + '_
    {
        self.sketch
            .constraints
            .iter()
            .filter(|c| !self.delete.contains(&c.id))
    }

    fn same_point(&self, a: SketchId, b: SketchId) -> bool {
        a == b
            || matches!(
                (self.classes.get(&a), self.classes.get(&b)),
                (Some(x), Some(y)) if x == y
            )
    }

    fn drop(&mut self, id: SketchId) {
        if self.delete.insert(id) {
            self.removed += 1;
        }
    }

    fn replace(&mut self, id: SketchId, with: Constraint<SketchId>) {
        self.delete.insert(id);
        self.pending.push(with);
    }

    /// Drops what measured a line as a segment: its length, equal lengths
    /// and midpoints (design note 0007 §7.2). Arcs keep everything: their
    /// constraints use the full circle.
    fn drop_segment_measures(&mut self) {
        let t = self.target;
        if !matches!(
            self.sketch.entity(t).map(|e| &e.geometry),
            Some(Geometry::Line { .. })
        ) {
            return;
        }
        let gone: Vec<SketchId> = self
            .live()
            .filter(|c| match &c.constraint {
                Constraint::Length { line, .. } | Constraint::Midpoint { line, .. } => *line == t,
                Constraint::Equal { a, b } => *a == t || *b == t,
                _ => false,
            })
            .map(|c| c.id)
            .collect();
        for id in gone {
            self.drop(id);
        }
    }

    /// Refuses to move end point `p` when another curve uses it too.
    fn check_own(&self, p: SketchId) -> Result<(), String> {
        let shared = self
            .sketch
            .entities
            .iter()
            .any(|e| e.id != self.target && e.geometry.referenced_points().contains(&p));
        if shared {
            Err("This end is shared with another curve and cannot move on its own".into())
        } else {
            Ok(())
        }
    }

    /// Frees end point `p` of the target, which is about to move: every
    /// constraint on it goes, curves that met there stay joined to one
    /// another, and a tangency at that end goes with it.
    fn detach(&mut self, p: SketchId) -> Result<(), String> {
        self.check_own(p)?;
        let t = self.target;
        let neighbours: Vec<SketchId> = self
            .live()
            .filter_map(|c| match c.constraint {
                Constraint::Coincident { a, b } if a == p => Some(b),
                Constraint::Coincident { a, b } if b == p => Some(a),
                _ => None,
            })
            .collect();
        let affected: Vec<(SketchId, Constraint<SketchId>)> = self
            .live()
            .filter(|c| {
                c.constraint.referenced().contains(&p)
                    || matches!(c.constraint, Constraint::Tangent {
                        a, b, kind: TangentKind::AtPoint { point } }
                        if (a == t || b == t) && self.same_point(point, p))
            })
            .map(|c| (c.id, c.constraint.clone()))
            .collect();
        for (id, c) in affected {
            match c {
                // Two other curves tangent where they met the target stay
                // tangent at the joint that remains.
                Constraint::Tangent {
                    a,
                    b,
                    kind: TangentKind::AtPoint { point },
                } if point == p && a != t && b != t && !neighbours.is_empty() => {
                    self.replace(
                        id,
                        Constraint::Tangent {
                            a,
                            b,
                            kind: TangentKind::AtPoint {
                                point: neighbours[0],
                            },
                        },
                    );
                }
                _ => self.drop(id),
            }
        }
        self.rejoin(&neighbours);
        Ok(())
    }

    /// Keeps `points` (which were all coincident through a point that is
    /// leaving) coincident with one another.
    fn rejoin(&mut self, points: &[SketchId]) {
        let Some((first, rest)) = points.split_first() else {
            return;
        };
        let mut uf = UnionFind::default();
        let links: Vec<(SketchId, SketchId)> = self
            .live()
            .filter_map(|c| match c.constraint {
                Constraint::Coincident { a, b } => Some((a, b)),
                _ => None,
            })
            .collect();
        for (a, b) in links {
            uf.union(a, b);
        }
        for q in rest {
            if uf.find(*first) != uf.find(*q) {
                uf.union(*first, *q);
                self.pending
                    .push(Constraint::Coincident { a: *first, b: *q });
                self.removed -= 1;
            }
        }
    }

    /// Adds the curve that continues the target from `at` to its old end
    /// `end`, which hands every constraint over to the new curve's end (and
    /// a tangency at that end over to the new curve). Arcs share the
    /// target's centre point.
    fn continuation(
        &mut self,
        end: SketchId,
        center: Option<SketchId>,
        at: DVec2,
    ) -> Result<Continuation, String> {
        self.check_own(end)?;
        let t = self.target;
        let old = self.sketch.point(end).ok_or("missing end point")?;
        let construction = self.sketch.entity(t).is_some_and(|e| e.construction);
        let start = self.batch.point(at);
        let new_end = self.batch.point(old);
        let geometry = match center {
            Some(center) => Geometry::Arc {
                center,
                start,
                end: new_end,
            },
            None => Geometry::Line {
                start,
                end: new_end,
            },
        };
        let curve = self.batch.entity(geometry, construction);
        let moved: Vec<(SketchId, Constraint<SketchId>)> = self
            .live()
            .filter_map(|c| {
                let joint = matches!(c.constraint, Constraint::Tangent {
                    a, b, kind: TangentKind::AtPoint { point } }
                    if (a == t || b == t) && self.same_point(point, end));
                (joint || c.constraint.referenced().contains(&end)).then(|| {
                    let m = c.constraint.map_keys(|k| {
                        if k == end {
                            new_end
                        } else if joint && k == t {
                            curve
                        } else {
                            k
                        }
                    });
                    (c.id, m)
                })
            })
            .collect();
        for (id, m) in moved {
            self.replace(id, m);
        }
        Ok(Continuation {
            curve,
            start,
            end: new_end,
        })
    }

    /// Attaches the target's new end `p` to a curve through `cut`.
    fn attach(&mut self, p: SketchId, cut: &Cut) {
        self.attach_for(p, cut, self.target);
    }

    /// Puts `p` on a curve through `cut` — one that is not constrained
    /// tangent to `curve` (there the two equations would have a double
    /// root, design note 0007 §7.2). With no such curve the end stays free
    /// along the cut curve.
    fn attach_for(&mut self, p: SketchId, cut: &Cut, curve: SketchId) {
        let tangent = |x: SketchId| {
            self.sketch.constraints.iter().any(|c| {
                matches!(c.constraint, Constraint::Tangent { a, b, .. }
                    if (a == curve && b == x) || (a == x && b == curve))
            })
        };
        if let Some(x) = cut.cutters.iter().copied().find(|x| !tangent(*x)) {
            self.pending
                .push(Constraint::PointOnCurve { point: p, curve: x });
        }
    }

    fn finish(mut self, preview: Preview) -> Plan {
        let mut edits = Vec::new();
        if !self.delete.is_empty() {
            edits.push(SketchEdit::Delete {
                items: self.delete.iter().copied().collect(),
            });
        }
        if !self.moves.is_empty() {
            edits.push(SketchEdit::MovePoints {
                points: self.moves.clone(),
            });
        }
        for c in std::mem::take(&mut self.pending) {
            self.batch.constraint(c);
        }
        edits.extend(self.batch.edits);
        Plan {
            edits,
            removed: self.removed,
            preview,
        }
    }
}

/// A small union-find over sketch ids.
#[derive(Default)]
struct UnionFind(BTreeMap<SketchId, SketchId>);

impl UnionFind {
    fn find(&mut self, x: SketchId) -> SketchId {
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

    fn union(&mut self, a: SketchId, b: SketchId) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0.insert(ra, rb);
        }
    }
}

/// Each point's coincidence class (points in no coincidence are absent).
fn coincidence_classes(sketch: &SketchData) -> BTreeMap<SketchId, usize> {
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

/// A line's or arc's start and end points.
fn end_points(c: &ResolvedCurve) -> Option<(DVec2, DVec2)> {
    match c {
        ResolvedCurve::Line { start, end } => Some((*start, *end)),
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => Some((
            *center + DVec2::from_angle(*start_angle) * *radius,
            *center + DVec2::from_angle(start_angle + sweep_angle) * *radius,
        )),
        _ => None,
    }
}

/// The parameter of the point of a line, arc or circle nearest `p`.
fn param(c: &ResolvedCurve, p: DVec2) -> f64 {
    match c {
        ResolvedCurve::Line { start, end } => {
            let v = *end - *start;
            ((p - *start).dot(v) / v.length_squared()).clamp(0.0, 1.0)
        }
        ResolvedCurve::Arc {
            center,
            start_angle,
            sweep_angle,
            ..
        } => {
            let phi = ((p - *center).to_angle() - start_angle).rem_euclid(TAU);
            if phi <= *sweep_angle {
                phi
            } else if phi - sweep_angle < TAU - phi {
                *sweep_angle
            } else {
                0.0
            }
        }
        ResolvedCurve::Circle { center, .. } => (p - *center).to_angle().rem_euclid(TAU),
        ResolvedCurve::Spline { .. } => 0.0,
    }
}

/// The point of a line, arc or circle at parameter `t`.
fn point_at(c: &ResolvedCurve, t: f64) -> DVec2 {
    match c {
        ResolvedCurve::Line { start, end } => start.lerp(*end, t),
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            ..
        } => *center + DVec2::from_angle(start_angle + t) * *radius,
        ResolvedCurve::Circle { center, radius } => *center + DVec2::from_angle(t) * *radius,
        ResolvedCurve::Spline { points, .. } => points[0],
    }
}

/// The piece of a line, arc or circle between parameters `a` and `b`
/// (counter-clockwise from `a` to `b` on a circle).
fn piece(c: &ResolvedCurve, a: f64, b: f64) -> ResolvedCurve {
    match c {
        ResolvedCurve::Line { .. } => ResolvedCurve::Line {
            start: point_at(c, a.min(b)),
            end: point_at(c, a.max(b)),
        },
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            ..
        } => ResolvedCurve::Arc {
            center: *center,
            radius: *radius,
            start_angle: start_angle + a.min(b),
            sweep_angle: (b - a).abs(),
        },
        ResolvedCurve::Circle { center, radius } => ResolvedCurve::Arc {
            center: *center,
            radius: *radius,
            start_angle: a,
            sweep_angle: (b - a).rem_euclid(TAU),
        },
        ResolvedCurve::Spline { .. } => c.clone(),
    }
}

/// The diagonal of the box around all curves.
fn extent(curves: &[(SketchId, ResolvedCurve)]) -> f64 {
    let mut lo = DVec2::splat(f64::INFINITY);
    let mut hi = DVec2::splat(f64::NEG_INFINITY);
    let mut add = |p: DVec2| {
        lo = lo.min(p);
        hi = hi.max(p);
    };
    for (_, c) in curves {
        match c {
            ResolvedCurve::Line { start, end } => {
                add(*start);
                add(*end);
            }
            ResolvedCurve::Arc { center, radius, .. }
            | ResolvedCurve::Circle { center, radius } => {
                add(*center - DVec2::splat(*radius));
                add(*center + DVec2::splat(*radius));
            }
            ResolvedCurve::Spline { points, .. } => points.iter().for_each(|p| add(*p)),
        }
    }
    if lo.x <= hi.x {
        (hi - lo).length()
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(a: (f64, f64), b: (f64, f64)) -> ResolvedCurve {
        ResolvedCurve::Line {
            start: DVec2::new(a.0, a.1),
            end: DVec2::new(b.0, b.1),
        }
    }

    #[test]
    fn parameters_round_trip() {
        let arc = ResolvedCurve::Arc {
            center: DVec2::ZERO,
            radius: 2.0,
            start_angle: 3.0,
            sweep_angle: 1.5,
        };
        for t in [0.0, 0.4, 1.5] {
            assert!((param(&arc, point_at(&arc, t)) - t).abs() < 1e-12);
        }
        // Off the arc: the nearer end by angle.
        assert_eq!(param(&arc, DVec2::from_angle(3.0 + 1.6)), 1.5);
        assert_eq!(param(&arc, DVec2::from_angle(2.9)), 0.0);
        let l = line((0.0, 0.0), (4.0, 0.0));
        assert_eq!(param(&l, DVec2::new(1.0, 3.0)), 0.25);
        assert_eq!(param(&l, DVec2::new(-1.0, 3.0)), 0.0);
        let c = ResolvedCurve::Circle {
            center: DVec2::ZERO,
            radius: 1.0,
        };
        assert!((param(&c, DVec2::new(0.0, -1.0)) - 1.5 * std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn circle_pieces_run_counter_clockwise() {
        let c = ResolvedCurve::Circle {
            center: DVec2::ZERO,
            radius: 1.0,
        };
        match piece(&c, 5.0, 1.0) {
            ResolvedCurve::Arc {
                start_angle,
                sweep_angle,
                ..
            } => {
                assert_eq!(start_angle, 5.0);
                assert!((sweep_angle - (TAU - 4.0)).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn union_find_joins_classes() {
        let id = SketchId::new_in_batch;
        let mut uf = UnionFind::default();
        uf.union(id(1), id(2));
        uf.union(id(3), id(4));
        assert_ne!(uf.find(id(1)), uf.find(id(3)));
        uf.union(id(2), id(4));
        assert_eq!(uf.find(id(1)), uf.find(id(3)));
        assert_eq!(uf.find(id(9)), id(9));
    }

    #[test]
    fn extent_covers_every_curve() {
        let curves = vec![
            (SketchId::new_in_batch(0), line((0.0, 0.0), (3.0, 0.0))),
            (
                SketchId::new_in_batch(1),
                ResolvedCurve::Circle {
                    center: DVec2::new(3.0, 4.0),
                    radius: 1.0,
                },
            ),
        ];
        assert!((extent(&curves) - DVec2::new(4.0, 5.0).length()).abs() < 1e-12);
        assert_eq!(extent(&[]), 0.0);
    }
}
