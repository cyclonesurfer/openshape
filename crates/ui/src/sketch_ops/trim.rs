//! Trim, extend and split (design note 0007 §7): cut a curve where the
//! other curves of the sketch meet it, and keep every constraint that still
//! means the same thing. Every piece of a line, circle or ellipse is a line,
//! arc or elliptical arc on the same points, and every piece of a conic is
//! a conic, so cuts are exact (design notes 0008 §8.1, 0009 §8.1).
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
use openshape_kernel::{Kernel, ProfileCurve, Tolerance};
use openshape_sketch_solver::{conic, ResolvedCurve, RHO_MIN};

use super::shapes::{closest_anomaly, eccentric_anomaly, ellipse_point};
use super::{coincidence_classes, deletion, own_axes, text, Batch, UnionFind};

/// Cuts closer than this to one another, or to the cut curve's ends, are
/// one cut (metres): a piece a micron long is not a piece anyone meant.
pub const CUT_GAP: f64 = 1e-6;

/// Curves meeting at an angle whose sine is below this (about 0.6°)
/// touch rather than cross, see [`Cut::cutters`].
pub const TANGENT_SIN: f64 = 1e-2;

/// Where other curves cut a curve.
#[derive(Clone, Debug, PartialEq)]
pub struct Cut {
    /// Parameter along the cut curve: a line's fraction, the angle
    /// travelled along an arc from its start, a circle's angle from +x (the
    /// arrangement's conventions), the eccentric anomaly turned
    /// counter-clockwise in the plane from an elliptical arc's start or an
    /// ellipse's first axis end, and a conic's own parameter in [0, 1].
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
    /// New curves and points (mirror images).
    Shapes {
        /// Curves.
        curves: Vec<ResolvedCurve>,
        /// Points.
        points: Vec<DVec2>,
    },
}

/// Finds where the other curves of `sketch` cut `target` (`curves` are its
/// solved curves, construction geometry included). A text's letters cut
/// but cannot hold a cut end (ADR 0015); an ellipse's own semi-axes do not
/// cut it ([`own_axes`]).
pub fn cuts(
    kernel: &(impl Kernel + ?Sized),
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    target: SketchId,
) -> Result<Cuts, String> {
    let curve = curves
        .iter()
        .find(|(id, _)| *id == target)
        .map(|(_, c)| c.clone())
        .ok_or("the curve no longer exists")?;
    let a = kernel_curve(&curve).map_err(|e| e.to_string())?;
    let ends = end_points(&curve);
    let loose = text::texts(sketch);
    let own = own_axes(sketch, target);
    let mut raw = Vec::new();
    for (id, c) in curves {
        if *id == target || own.contains(id) {
            continue;
        }
        for x in contacts(kernel, &a, c)? {
            let at_end = ends.is_some_and(|(s, e)| {
                (x.point - s).length() <= CUT_GAP || (x.point - e).length() <= CUT_GAP
            });
            if !at_end {
                let mut cut = x.cut(*id, &loose);
                // The kernel's parameter on an ellipse's rational form runs
                // the same way as the anomaly, but not evenly.
                if let ResolvedCurve::Ellipse { .. } | ResolvedCurve::EllipticalArc { .. } = curve {
                    cut.t = turned_to(&curve, cut.point).expect("an ellipse");
                }
                raw.push(cut);
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
    if matches!(
        curve,
        ResolvedCurve::Circle { .. } | ResolvedCurve::Ellipse { .. }
    ) && out.len() >= 2
    {
        // The first and last cut of a closed curve can meet across its
        // seam.
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
    /// The cut it makes, by curve `id` (which can hold the cut end unless
    /// it is `loose`).
    fn cut(&self, id: SketchId, loose: &BTreeSet<SketchId>) -> Cut {
        Cut {
            t: self.t,
            point: self.point,
            cutters: if self.crossing && !loose.contains(&id) {
                vec![id]
            } else {
                Vec::new()
            },
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
    let Ok(b) = kernel_curve(b) else {
        return Ok(Vec::new());
    };
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
            let center = *center;
            trim_closed(sketch, cuts, near, "A circle", |start, end| Geometry::Arc {
                center,
                start,
                end,
            })
        }
        (
            &Geometry::Ellipse {
                center,
                major,
                minor,
            },
            ResolvedCurve::Ellipse { .. },
        ) => trim_closed(sketch, cuts, near, "An ellipse", |start, end| {
            Geometry::EllipticalArc {
                center,
                major,
                minor,
                start,
                end,
            }
        }),
        (
            Geometry::Line { start, end }
            | Geometry::Arc { start, end, .. }
            | Geometry::EllipticalArc { start, end, .. },
            _,
        ) => {
            let carrier = Carrier::of(geometry).ok_or("this curve cannot be trimmed")?;
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
                let d = r.continuation(*end, carrier, hi.point)?;
                r.moves.push((*end, lo.point));
                r.attach(*end, lo);
                r.attach(d.start, hi);
                match carrier {
                    Carrier::Line => {
                        r.pending.push(Constraint::PointOnCurve {
                            point: d.start,
                            curve: t,
                        });
                        r.pending.push(Constraint::PointOnCurve {
                            point: d.end,
                            curve: t,
                        });
                    }
                    Carrier::Circle(_) => r.pending.push(Constraint::Equal { a: t, b: d.curve }),
                    // The new piece's own points define its ellipse
                    // (conics are trimmed by `trim_conic`).
                    Carrier::Ellipse { .. } | Carrier::Conic { .. } => {}
                }
                Ok(r.finish(Preview::Piece(piece(&cuts.curve, lo.t, hi.t))))
            }
        }
        (Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. }, _) => Err(
            "A spline can only be trimmed whole: its pieces are not exact splines of the \
             sketch's kinds (design note 0007 §7.5). Delete it, or redraw it in pieces"
                .into(),
        ),
        (&Geometry::Conic { control, end, .. }, &ResolvedCurve::Conic { .. }) => {
            trim_conic(sketch, cuts, near, control, end)
        }
        _ => Err("this curve cannot be trimmed".into()),
    }
}

/// Trims a conic (design note 0009 §8.1): what is kept is one or two exact
/// pieces of it ([`conic::piece`]), each a conic with its own control
/// point (where its end tangents meet) and ρ. The old control point moves
/// to the kept piece's; constraints on it described the old corner and
/// go, as does a ρ dimension.
fn trim_conic(
    sketch: &SketchData,
    cuts: &Cuts,
    near: DVec2,
    control: SketchId,
    end: SketchId,
) -> Result<Plan, String> {
    let t = cuts.target;
    let (p, w) = conic_of(&cuts.curve).ok_or("the conic no longer exists")?;
    let k = cuts.cuts.len();
    let tn = param(&cuts.curve, near);
    let i = cuts.cuts.iter().filter(|c| c.t < tn).count();
    let mut r = Rewrite::new(sketch, t);
    r.drop_rho();
    r.detach(control)?;
    if i == 0 || i == k {
        // An end piece: that end moves to the cut next to it.
        let (moved, cut, kept) = if i == 0 {
            let start = match sketch.entity(t).map(|e| &e.geometry) {
                Some(&Geometry::Conic { start, .. }) => start,
                _ => return Err("the conic no longer exists".into()),
            };
            (start, &cuts.cuts[0], (cuts.cuts[0].t, 1.0))
        } else {
            (end, &cuts.cuts[k - 1], (0.0, cuts.cuts[k - 1].t))
        };
        let (q, wq) = conic_piece(p, w, kept.0, kept.1)?;
        r.detach(moved)?;
        r.moves.push((moved, cut.point));
        r.moves.push((control, q[1]));
        r.set_rho(t, wq)?;
        r.attach(moved, cut);
        r.points_on_conic(p, w, &[(kept.0, kept.1, t)]);
        let gone = if i == 0 { (0.0, cut.t) } else { (cut.t, 1.0) };
        Ok(r.finish(Preview::Piece(piece(&cuts.curve, gone.0, gone.1))))
    } else {
        // A middle piece: the conic keeps the part before it; a new conic
        // takes the part after it and the old end.
        let (lo, hi) = (&cuts.cuts[i - 1], &cuts.cuts[i]);
        let (qa, wa) = conic_piece(p, w, 0.0, lo.t)?;
        let (qb, wb) = conic_piece(p, w, hi.t, 1.0)?;
        let carrier = Carrier::Conic {
            control: qb[1],
            rho: rho_of_weight(wb)?,
        };
        let d = r.continuation(end, carrier, hi.point)?;
        r.moves.push((end, lo.point));
        r.moves.push((control, qa[1]));
        r.set_rho(t, wa)?;
        r.attach(end, lo);
        r.attach(d.start, hi);
        r.points_on_conic(p, w, &[(0.0, lo.t, t), (hi.t, 1.0, d.curve)]);
        Ok(r.finish(Preview::Piece(piece(&cuts.curve, lo.t, hi.t))))
    }
}

const NO_PIECE: &str = "This piece of the conic is not a conic arc of its own";

/// The piece of conic arc (`p`, `w`) from parameter `a` to `b`
/// ([`conic::piece`]) as a conic of the sketch: refused when it is no conic
/// arc of its own, or when its control point lies within the model
/// tolerance of its chord. Such a piece (a few microns of a gently curved
/// conic) is a straight segment to the tolerance: its ρ means nothing and
/// the solver can flatten it altogether.
fn conic_piece(p: [DVec2; 3], w: f64, a: f64, b: f64) -> Result<([DVec2; 3], f64), String> {
    let (q, wq) = conic::piece(p, w, a, b).ok_or(NO_PIECE)?;
    let chord = q[2] - q[0];
    let height = chord.perp_dot(q[1] - q[0]).abs() / chord.length();
    if height > Tolerance::default().linear {
        Ok((q, wq))
    } else {
        Err(
            "The conic's piece would be too short to keep its shape: its corner would lie \
             within the model tolerance of its chord"
                .into(),
        )
    }
}

/// A conic's control points and middle weight.
fn conic_of(c: &ResolvedCurve) -> Option<([DVec2; 3], f64)> {
    match *c {
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => Some(([start, control, end], weight)),
        _ => None,
    }
}

/// How far past its end (t = 1) conic arc (`p`, `w`) can grow as one conic
/// arc whose control point and end stay within `reach` of its start: the
/// last t of a geometric scan of t − 1 (from 1e-4, by factors of 1.25)
/// before the first where the arc from 0 to t is no conic arc of the
/// sketch (it would turn half a turn, pass through infinity, or have a ρ
/// out of range, [`rho_of_weight`]) or reaches too far. `None` when it
/// cannot grow by 1e-4 of its parameter.
fn conic_reach(p: [DVec2; 3], w: f64, reach: f64) -> Option<f64> {
    let near = |x: DVec2| (x - p[0]).length() <= reach;
    let ok = |t: f64| {
        matches!(conic::piece(p, w, 0.0, t),
            Some((q, wq)) if near(q[1]) && near(q[2]) && rho_of_weight(wq).is_ok())
            && conic::piece(p, w, 1.0, t).is_some()
    };
    let mut best = None;
    let mut h = 1e-4;
    while h < 1e12 && ok(1.0 + h) {
        best = Some(1.0 + h);
        h *= 1.25;
    }
    best
}

/// The ρ of middle weight `w`, when a conic of the sketch can have it
/// (design note 0009 §4.6): a piece cut very short is all but straight.
fn rho_of_weight(w: f64) -> Result<f64, String> {
    let rho = w / (1.0 + w);
    if (RHO_MIN..=1.0 - RHO_MIN).contains(&rho) {
        Ok(rho)
    } else {
        Err(format!(
            "The conic's piece would be too {} to keep (ρ {rho:.1e})",
            if rho < 0.5 { "flat" } else { "sharp" }
        ))
    }
}

/// A closed curve (circle or ellipse) cut at two or more places becomes an
/// arc of it (`arc(start, end)`, on the same defining points) with every
/// constraint of the closed curve: they all use the full curve. `noun`
/// names the curve in the refusal.
fn trim_closed(
    sketch: &SketchData,
    cuts: &Cuts,
    near: DVec2,
    noun: &str,
    arc: impl FnOnce(SketchId, SketchId) -> Geometry<SketchId>,
) -> Result<Plan, String> {
    let k = cuts.cuts.len();
    if k < 2 {
        return Err(format!(
            "{noun} needs two crossings to trim: one crossing leaves nothing to cut away"
        ));
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
    let arc = r.batch.entity(arc(start, end), construction);
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

/// Plans extending the end of line, arc, elliptical arc or conic `target`
/// nearest `near` to the first curve beyond it (design note 0007 §7.3); the
/// new end is not held on a text's letters (see [`cuts`]).
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
    let geometry = &sketch
        .entity(target)
        .ok_or("the curve no longer exists")?
        .geometry;
    let (start, end) = match geometry {
        Geometry::Line { start, end }
        | Geometry::Arc { start, end, .. }
        | Geometry::EllipticalArc { start, end, .. }
        | Geometry::Conic { start, end, .. } => (*start, *end),
        g => return Err(cannot("extended", g)),
    };
    let (s, e) = end_points(&curve).ok_or("this curve has no ends to extend")?;
    let at_end = (near - e).length() <= (near - s).length();
    let (p, from, other) = if at_end { (end, e, s) } else { (start, s, e) };
    let reach = 4.0 * extent(curves).max(from.length()).max(1e-3);
    // A conic arc grows at its end as the arc of the same conic from its
    // other end; at its start, as the reversed arc does at its end.
    let grown = conic_of(&curve).map(|(c, w)| {
        let c = if at_end { c } else { [c[2], c[1], c[0]] };
        (c, w)
    });
    // The way the end can go: a long ray, or the rest of the arc's circle
    // or ellipse, starting at the end.
    let proxy = match &curve {
        ResolvedCurve::Line { .. } => {
            let d = (from - other).normalize();
            ProfileCurve::Line {
                start: from,
                end: from + d * reach,
            }
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
            if at_end {
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
            }
        }
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            start_param,
            sweep_param,
        } => {
            let rest = TAU - sweep_param.abs();
            if rest * major.length().min(minor.length()) <= 2.0 * CUT_GAP {
                return Err("The elliptical arc is already a full ellipse".into());
            }
            let h = sweep_param.signum();
            let (t0, sweep) = if at_end {
                (start_param + sweep_param, h * rest)
            } else {
                (*start_param, -h * rest)
            };
            ProfileCurve::elliptical_arc(*center, *major, *minor, t0, sweep)
                .map_err(|e| e.to_string())?
        }
        ResolvedCurve::Conic { .. } => {
            let (c, w) = grown.expect("a conic");
            let far = conic_reach(c, w, reach).ok_or(
                "This conic cannot be extended: a longer arc of it would not be one conic arc \
                 (it would turn half a turn or more, or run off to infinity)",
            )?;
            let (q, wq) = conic::piece(c, w, 1.0, far).ok_or(NO_PIECE)?;
            ProfileCurve::conic(q[0], q[1], q[2], wq).map_err(|e| e.to_string())?
        }
        _ => return Err("this curve has no ends to extend".into()),
    };
    // The nearest contact beyond the end (not the curve's own other end),
    // by the parameter along the proxy, which runs away from the end.
    let own = own_axes(sketch, target);
    let loose = text::texts(sketch);
    let mut best: Option<(f64, Cut)> = None;
    for (id, c) in curves {
        if *id == target || own.contains(id) {
            continue;
        }
        for x in contacts(kernel, &proxy, c)? {
            if (x.point - from).length() <= CUT_GAP || (x.point - other).length() <= CUT_GAP {
                continue;
            }
            match &mut best {
                Some((_, b)) if (b.point - x.point).length() <= CUT_GAP => {
                    merge(b, x.cut(*id, &loose))
                }
                Some((t, _)) if *t <= x.t => {}
                _ => best = Some((x.t, x.cut(*id, &loose))),
            }
        }
    }
    let (_, cut) = best.ok_or("Nothing to extend to: no curve lies beyond this end")?;
    let mut r = Rewrite::new(sketch, target);
    r.drop_segment_measures();
    r.detach(p)?;
    r.moves.push((p, cut.point));
    r.attach(p, &cut);
    if let (Some((c, w)), &Geometry::Conic { control, .. }) = (grown, geometry) {
        // The longer arc of the same conic: its control point moves along
        // the kept end's tangent to where the new end's tangent meets it.
        let tx = conic::parameter_of(c, w, cut.point);
        let (q, wq) = conic_piece(c, w, 0.0, tx)?;
        let (added, wa) = conic::piece(c, w, 1.0, tx).ok_or(NO_PIECE)?;
        r.drop_rho();
        r.detach(control)?;
        r.moves.push((control, q[1]));
        r.set_rho(target, wq)?;
        let [start, control, end] = added;
        return Ok(r.finish(Preview::Piece(ResolvedCurve::Conic {
            start,
            control,
            end,
            weight: wa,
        })));
    }
    // The piece added, from the old end to the cut.
    let preview = match turned_to(&curve, cut.point) {
        Some(turned) if at_end => piece(&curve, span(&curve), turned),
        Some(turned) => piece(&curve, turned - TAU, 0.0),
        None => ResolvedCurve::Line {
            start: from,
            end: cut.point,
        },
    };
    Ok(r.finish(Preview::Piece(preview)))
}

/// Plans splitting line, arc, elliptical arc or conic `target` at its
/// point nearest `near` into two curves joined there (design note 0007
/// §7.4).
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
    let geometry = &sketch
        .entity(target)
        .ok_or("the curve no longer exists")?
        .geometry;
    let end = match geometry {
        Geometry::Line { end, .. }
        | Geometry::Arc { end, .. }
        | Geometry::EllipticalArc { end, .. }
        | Geometry::Conic { end, .. } => *end,
        g => return Err(cannot("split", g)),
    };
    let (s, e) = end_points(&curve).ok_or("this curve has no ends to split between")?;
    let tq = param(&curve, near);
    let q = point_at(&curve, tq);
    if (q - s).length() <= CUT_GAP || (q - e).length() <= CUT_GAP {
        return Err("Too close to an end to split".into());
    }
    let mut r = Rewrite::new(sketch, target);
    r.drop_segment_measures();
    // A conic's two pieces are conics of their own (design note 0009
    // §8.1), the first with the old control point moved to its corner.
    let conic_pieces = match (geometry, conic_of(&curve)) {
        (&Geometry::Conic { control, .. }, Some((p, w))) => {
            let (qa, wa) = conic_piece(p, w, 0.0, tq)?;
            let (qb, wb) = conic_piece(p, w, tq, 1.0)?;
            r.drop_rho();
            r.detach(control)?;
            r.moves.push((control, qa[1]));
            r.set_rho(target, wa)?;
            Some((p, w, qb[1], rho_of_weight(wb)?))
        }
        _ => None,
    };
    let carrier = match conic_pieces {
        Some((_, _, control, rho)) => Carrier::Conic { control, rho },
        None => Carrier::of(geometry).ok_or("this curve cannot be split")?,
    };
    let d = r.continuation(end, carrier, q)?;
    r.moves.push((end, q));
    // An arc's new piece shares its centre and has its own radius; an
    // elliptical arc's shares its ellipse, which holds the joint once
    // (design note 0008 §8.1); a conic's pieces are held tangent at the
    // joint, as they were.
    r.pending
        .push(Constraint::Coincident { a: end, b: d.start });
    match carrier {
        Carrier::Line => r.pending.push(Constraint::PointOnCurve {
            point: d.end,
            curve: target,
        }),
        Carrier::Conic { .. } => {
            r.pending.push(Constraint::Tangent {
                a: target,
                b: d.curve,
                kind: TangentKind::AtPoint { point: end },
            });
            let (p, w, ..) = conic_pieces.expect("a conic");
            r.points_on_conic(p, w, &[(0.0, tq, target), (tq, 1.0, d.curve)]);
        }
        Carrier::Circle(_) | Carrier::Ellipse { .. } => {}
    }
    Ok(r.finish(Preview::Point(q)))
}

/// Why `g` cannot be `verb` ("extended", "split").
fn cannot(verb: &str, g: &Geometry<SketchId>) -> String {
    let closed = match g {
        Geometry::Circle { .. } => Some("A circle"),
        Geometry::Ellipse { .. } => Some("An ellipse"),
        _ => None,
    };
    if let Some(noun) = closed {
        return if verb == "split" {
            format!("{noun} has no ends to split between; trim it instead")
        } else {
            format!("{noun} has no ends to extend")
        };
    }
    match g {
        Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. } => format!(
            "A spline cannot be {verb}: its pieces and extensions are not exact splines of the \
             sketch's kinds (design note 0007 §7.5)"
        ),
        _ => format!("Only lines, arcs, elliptical arcs and conics can be {verb}"),
    }
}

/// What a cut curve lies on, which a new piece of it shares.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Carrier {
    /// A line: a new piece is held on it by constraints.
    Line,
    /// An arc's circle: a new piece shares the centre point.
    Circle(SketchId),
    /// An elliptical arc's ellipse: a new piece shares its three points,
    /// which define the ellipse.
    Ellipse {
        center: SketchId,
        major: SketchId,
        minor: SketchId,
    },
    /// A conic: a new piece is a conic of its own with this control point
    /// and ρ, the piece of the old one (design note 0009 §8.1).
    Conic { control: DVec2, rho: f64 },
}

impl Carrier {
    fn of(g: &Geometry<SketchId>) -> Option<Carrier> {
        match *g {
            Geometry::Line { .. } => Some(Carrier::Line),
            Geometry::Arc { center, .. } => Some(Carrier::Circle(center)),
            Geometry::EllipticalArc {
                center,
                major,
                minor,
                ..
            } => Some(Carrier::Ellipse {
                center,
                major,
                minor,
            }),
            _ => None,
        }
    }

    /// A piece from `start` to `end` on it (adding any other point it
    /// needs to `batch`).
    fn piece(self, batch: &mut Batch, start: SketchId, end: SketchId) -> Geometry<SketchId> {
        match self {
            Carrier::Line => Geometry::Line { start, end },
            Carrier::Circle(center) => Geometry::Arc { center, start, end },
            Carrier::Ellipse {
                center,
                major,
                minor,
            } => Geometry::EllipticalArc {
                center,
                major,
                minor,
                start,
                end,
            },
            Carrier::Conic { control, rho } => Geometry::Conic {
                start,
                control: batch.point(control),
                end,
                rho,
            },
        }
    }
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
    /// Edits of existing entities after the moves (a conic's new ρ).
    extra: Vec<SketchEdit>,
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
            extra: Vec::new(),
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
                    || c.constraint.joint().is_some_and(|(a, b, point)| {
                        (a == t || b == t) && self.same_point(point, p)
                    })
            })
            .map(|c| (c.id, c.constraint.clone()))
            .collect();
        for (id, c) in affected {
            match c.joint() {
                // Two other curves tangent (or curvature continuous) where
                // they met the target stay so at the joint that remains.
                Some((a, b, point)) if point == p && a != t && b != t && !neighbours.is_empty() => {
                    let moved = c.map_keys(|k| if k == p { neighbours[0] } else { k });
                    self.replace(id, moved);
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
    /// target's centre point, elliptical arcs its three ellipse points.
    fn continuation(
        &mut self,
        end: SketchId,
        carrier: Carrier,
        at: DVec2,
    ) -> Result<Continuation, String> {
        self.check_own(end)?;
        let t = self.target;
        let old = self.sketch.point(end).ok_or("missing end point")?;
        let construction = self.sketch.entity(t).is_some_and(|e| e.construction);
        let start = self.batch.point(at);
        let new_end = self.batch.point(old);
        let geometry = carrier.piece(&mut self.batch, start, new_end);
        let curve = self.batch.entity(geometry, construction);
        let moved: Vec<(SketchId, Constraint<SketchId>)> =
            self.live()
                .filter_map(|c| {
                    let joint = c.constraint.joint().is_some_and(|(a, b, point)| {
                        (a == t || b == t) && self.same_point(point, end)
                    });
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

    /// Drops the target's ρ dimension: a piece of a conic has its own ρ.
    fn drop_rho(&mut self) {
        let t = self.target;
        let gone: Vec<SketchId> = self
            .live()
            .filter(|c| matches!(c.constraint, Constraint::Rho { curve, .. } if curve == t))
            .map(|c| c.id)
            .collect();
        for id in gone {
            self.drop(id);
        }
    }

    /// Sets the ρ of conic `curve` from a piece's middle weight `w`,
    /// refusing a piece too flat or too sharp to be a conic of the sketch
    /// (design note 0009 §4.6).
    fn set_rho(&mut self, curve: SketchId, w: f64) -> Result<f64, String> {
        let rho = rho_of_weight(w)?;
        self.extra.push(SketchEdit::SetRho { entity: curve, rho });
        Ok(rho)
    }

    /// Hands the target conic's point-on-curve constraints to the pieces
    /// that keep their points: `pieces` are parameter ranges of the old arc
    /// (`p`, `w`) with the curve that keeps each. A constraint whose point
    /// was on a removed part goes: on a conic it means on the arc.
    fn points_on_conic(&mut self, p: [DVec2; 3], w: f64, pieces: &[(f64, f64, SketchId)]) {
        let t = self.target;
        let on: Vec<(SketchId, SketchId)> = self
            .live()
            .filter_map(|c| match c.constraint {
                Constraint::PointOnCurve { point, curve } if curve == t => Some((c.id, point)),
                _ => None,
            })
            .collect();
        for (id, point) in on {
            let Some(at) = self.sketch.point(point) else {
                continue;
            };
            let tq = conic::closest_param(p, w, at);
            // A point at a piece's end (within the cut gap) is on it.
            let holder = pieces.iter().find(|(a, b, _)| {
                (a.min(*b)..=a.max(*b)).contains(&tq)
                    || [a, b]
                        .iter()
                        .any(|x| (conic::point(p, w, **x) - at).length() <= CUT_GAP)
            });
            match holder {
                Some(&(_, _, c)) if c == t => {}
                Some(&(_, _, c)) => self.replace(id, Constraint::PointOnCurve { point, curve: c }),
                None => self.drop(id),
            }
        }
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
                    | Constraint::Curvature { a, b, .. }
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
        edits.append(&mut self.extra);
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
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            start_param,
            sweep_param,
        } => {
            let at = |t: f64| *center + *major * t.cos() + *minor * t.sin();
            Some((at(*start_param), at(start_param + sweep_param)))
        }
        ResolvedCurve::Conic { start, end, .. } => Some((*start, *end)),
        _ => None,
    }
}

/// The parameter ([`Cut::t`]) of the point of a line, arc, circle,
/// ellipse, elliptical arc or conic nearest `p`, for a `p` near it.
fn param(c: &ResolvedCurve, p: DVec2) -> f64 {
    match c {
        ResolvedCurve::Line { start, end } => {
            let v = *end - *start;
            ((p - *start).dot(v) / v.length_squared()).clamp(0.0, 1.0)
        }
        ResolvedCurve::Arc { sweep_angle, .. } => {
            on_span(turned_to(c, p).expect("an arc"), *sweep_angle)
        }
        ResolvedCurve::Circle { center, .. } => (p - *center).to_angle().rem_euclid(TAU),
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        } => turned(c, closest_anomaly(*center, *major, *minor, p)),
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            sweep_param,
            ..
        } => on_span(
            turned(c, closest_anomaly(*center, *major, *minor, p)),
            sweep_param.abs(),
        ),
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => conic::closest_param([*start, *control, *end], *weight, p),
        ResolvedCurve::Spline { .. } => 0.0,
    }
}

/// `phi` (turned from an arc's start, in [0, 2π)) where it is on the arc
/// of extent `span`, else the nearer end by angle.
fn on_span(phi: f64, span: f64) -> f64 {
    if phi <= span {
        phi
    } else if phi - span < TAU - phi {
        span
    } else {
        0.0
    }
}

/// The angle (arcs) or eccentric anomaly (ellipses, elliptical arcs)
/// turned counter-clockwise in the plane from the curve's start to point
/// `p` on it (or on its full circle or ellipse), in [0, 2π). The start of a
/// full ellipse is its first axis end.
fn turned_to(c: &ResolvedCurve, p: DVec2) -> Option<f64> {
    match c {
        ResolvedCurve::Arc {
            center,
            start_angle,
            ..
        } => Some(((p - *center).to_angle() - start_angle).rem_euclid(TAU)),
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        }
        | ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            ..
        } => Some(turned(c, eccentric_anomaly(*center, *major, *minor, p))),
        _ => None,
    }
}

/// The anomaly turned counter-clockwise in the plane from an ellipse's
/// (or elliptical arc's) start to its point at anomaly `theta`, in
/// [0, 2π); see [`anomaly`].
fn turned(c: &ResolvedCurve, theta: f64) -> f64 {
    let (t0, h) = anomaly_frame(c);
    (h * (theta - t0)).rem_euclid(TAU)
}

/// The eccentric anomaly of an ellipse's (or elliptical arc's) point
/// turned `phi` counter-clockwise in the plane from its start.
fn anomaly(c: &ResolvedCurve, phi: f64) -> f64 {
    let (t0, h) = anomaly_frame(c);
    t0 + h * phi
}

/// An ellipse's start anomaly (0, the first axis end; an arc's start) and
/// the way anomaly runs in the plane: 1 counter-clockwise, −1 clockwise.
fn anomaly_frame(c: &ResolvedCurve) -> (f64, f64) {
    match c {
        ResolvedCurve::Ellipse { major, minor, .. } => (
            0.0,
            if major.perp_dot(*minor) < 0.0 {
                -1.0
            } else {
                1.0
            },
        ),
        ResolvedCurve::EllipticalArc {
            start_param,
            sweep_param,
            ..
        } => (*start_param, sweep_param.signum()),
        _ => (0.0, 1.0),
    }
}

/// An arc's or elliptical arc's extent in its parameter ([`Cut::t`]).
fn span(c: &ResolvedCurve) -> f64 {
    match c {
        ResolvedCurve::Arc { sweep_angle, .. } => *sweep_angle,
        ResolvedCurve::EllipticalArc { sweep_param, .. } => sweep_param.abs(),
        _ => 1.0,
    }
}

/// The point of a line, arc, circle, ellipse, elliptical arc or conic at
/// parameter `t` ([`Cut::t`]).
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
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        }
        | ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            ..
        } => ellipse_point(*center, *major, *minor, anomaly(c, t)),
        ResolvedCurve::Spline { points, .. } => points[0],
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => conic::point([*start, *control, *end], *weight, t),
    }
}

/// The piece of a line, arc, circle, ellipse, elliptical arc or conic
/// between parameters `a` and `b` (counter-clockwise from `a` to `b` on a
/// closed curve; a conic's piece is the exact conic arc, [`conic::piece`]).
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
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        } => {
            let (_, h) = anomaly_frame(c);
            ResolvedCurve::EllipticalArc {
                center: *center,
                major: *major,
                minor: *minor,
                start_param: anomaly(c, a),
                sweep_param: h * (b - a).rem_euclid(TAU),
            }
        }
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            sweep_param,
            ..
        } => ResolvedCurve::EllipticalArc {
            center: *center,
            major: *major,
            minor: *minor,
            start_param: anomaly(c, a.min(b)),
            sweep_param: sweep_param.signum() * (b - a).abs(),
        },
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => match conic::piece([*start, *control, *end], *weight, a.min(b), a.max(b)) {
            Some(([start, control, end], weight)) => ResolvedCurve::Conic {
                start,
                control,
                end,
                weight,
            },
            None => c.clone(),
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
            // Inside its control triangle.
            ResolvedCurve::Conic {
                start,
                control,
                end,
                ..
            } => [start, control, end].into_iter().for_each(|p| add(*p)),
            ResolvedCurve::Ellipse {
                center,
                major,
                minor,
            }
            | ResolvedCurve::EllipticalArc {
                center,
                major,
                minor,
                ..
            } => {
                // The box of the ellipse: its half-widths are the norms of
                // the axes' x and y components.
                let half = DVec2::new(major.x.hypot(minor.x), major.y.hypot(minor.y));
                add(*center - half);
                add(*center + half);
            }
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
