//! Mirror (design note 0007 §9): images of curves and points about a line,
//! held by one derived mirror constraint per line.
//!
//! The images get new points at the reflected positions and no constraints
//! of their own: the mirror constraint fixes every image point from its
//! source, so an image adds no degrees of freedom and a dimension on it
//! would be redundant. A source point on the mirror line is shared with the
//! image (the profile stays joined there) and held on the line; a
//! point-on-line constraint that said the same is replaced. Mirroring again
//! about the same line extends its constraint, reusing the images of points
//! already mirrored, so curves mirrored one at a time join up as they did
//! when they were sources.

use std::collections::{BTreeMap, BTreeSet};

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{SketchData, SketchEdit, SketchId};
use openshape_sketch_solver::ResolvedCurve;

use super::trim::{Plan, Preview};
use super::{coincidence_classes, Batch, SAME_POINT};

/// A mirror line: a point on it and its unit direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axis {
    /// A point of the line (its start).
    pub at: DVec2,
    /// Unit direction.
    pub dir: DVec2,
}

impl Axis {
    /// The line entity `axis` of `sketch`, or why it cannot be a mirror line.
    pub fn of(sketch: &SketchData, axis: SketchId) -> Result<Axis, String> {
        let Some(Geometry::Line { start, end }) = sketch.entity(axis).map(|e| &e.geometry) else {
            return Err("The mirror line must be a line".into());
        };
        let (a, b) = (point(sketch, *start)?, point(sketch, *end)?);
        let d = b - a;
        let len = d.length();
        if !(len.is_finite() && len > SAME_POINT) {
            return Err("The mirror line has no length".into());
        }
        Ok(Axis {
            at: a,
            dir: d / len,
        })
    }

    /// The reflection of `p`.
    pub fn reflect(&self, p: DVec2) -> DVec2 {
        let v = p - self.at;
        self.at + self.dir * (2.0 * v.dot(self.dir)) - v
    }

    /// The reflection of the vector `v` (a direction, not a position).
    pub fn reflect_vector(&self, v: DVec2) -> DVec2 {
        self.dir * (2.0 * v.dot(self.dir)) - v
    }

    /// The distance of `p` from the (infinite) line.
    pub fn distance(&self, p: DVec2) -> f64 {
        self.dir.perp_dot(p - self.at).abs()
    }

    /// The reflection of a solved curve; arcs keep running counter-clockwise
    /// (from the image of the end to the image of the start). An ellipse's
    /// image has the images of its semi-axes, so the same eccentric anomaly
    /// gives the image point, with the opposite handedness: an elliptical
    /// arc's image starts where its source ends, with the opposite sweep.
    /// A conic's image is the conic on the images of its points, with the
    /// same weight (design note 0009 §4.5).
    pub fn reflect_curve(&self, c: &ResolvedCurve) -> ResolvedCurve {
        match c {
            ResolvedCurve::Line { start, end } => ResolvedCurve::Line {
                start: self.reflect(*start),
                end: self.reflect(*end),
            },
            ResolvedCurve::Circle { center, radius } => ResolvedCurve::Circle {
                center: self.reflect(*center),
                radius: *radius,
            },
            ResolvedCurve::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => {
                let c = self.reflect(*center);
                let end = *center + DVec2::from_angle(start_angle + sweep_angle) * *radius;
                ResolvedCurve::Arc {
                    center: c,
                    radius: *radius,
                    start_angle: (self.reflect(end) - c).to_angle(),
                    sweep_angle: *sweep_angle,
                }
            }
            ResolvedCurve::Ellipse {
                center,
                major,
                minor,
            } => ResolvedCurve::Ellipse {
                center: self.reflect(*center),
                major: self.reflect_vector(*major),
                minor: self.reflect_vector(*minor),
            },
            ResolvedCurve::EllipticalArc {
                center,
                major,
                minor,
                start_param,
                sweep_param,
            } => ResolvedCurve::EllipticalArc {
                center: self.reflect(*center),
                major: self.reflect_vector(*major),
                minor: self.reflect_vector(*minor),
                start_param: start_param + sweep_param,
                sweep_param: -sweep_param,
            },
            ResolvedCurve::Spline {
                degree,
                knots,
                points,
            } => ResolvedCurve::Spline {
                degree: *degree,
                knots: knots.clone(),
                points: points.iter().map(|p| self.reflect(*p)).collect(),
            },
            ResolvedCurve::Conic {
                start,
                control,
                end,
                weight,
            } => ResolvedCurve::Conic {
                start: self.reflect(*start),
                control: self.reflect(*control),
                end: self.reflect(*end),
                weight: *weight,
            },
        }
    }
}

fn point(sketch: &SketchData, id: SketchId) -> Result<DVec2, String> {
    sketch
        .point(id)
        .ok_or_else(|| format!("point {id:?} does not exist"))
}

/// The mirror constraints about `axis` in `sketch`, with their pairs.
fn mirrors_about(
    sketch: &SketchData,
    axis: SketchId,
) -> Vec<(SketchId, Vec<(SketchId, SketchId)>)> {
    sketch
        .constraints
        .iter()
        .filter_map(|c| match &c.constraint {
            Constraint::Mirror { axis: a, pairs } if *a == axis => Some((c.id, pairs.clone())),
            _ => None,
        })
        .collect()
}

/// The points a mirror pairs up for the pair `(s, i)`, in order (an arc's
/// ends cross over, as in the solver; a point entity is its own point).
fn point_pairs(sketch: &SketchData, s: SketchId, i: SketchId) -> Vec<(SketchId, SketchId)> {
    let (Some(gs), Some(gi)) = (sketch.entity(s), sketch.entity(i)) else {
        return Vec::new();
    };
    let points = |id: SketchId, g: &Geometry<SketchId>| match g {
        Geometry::Point { .. } => vec![id],
        Geometry::Arc { center, start, end } if id == s => vec![*center, *end, *start],
        g => g.referenced_points(),
    };
    points(s, &gs.geometry)
        .into_iter()
        .zip(points(i, &gi.geometry))
        .collect()
}

/// The entities to mirror for a selection: its curves, and its points that
/// are not points of a selected curve (those come with their curve).
pub fn sources(sketch: &SketchData, selection: &[SketchId]) -> Vec<SketchId> {
    let picked: Vec<&openshape_core_model::sketch::Entity<SketchId>> = selection
        .iter()
        .filter_map(|id| sketch.entity(*id))
        .collect();
    let of_curves: BTreeSet<SketchId> = picked
        .iter()
        .filter(|e| e.geometry.is_curve())
        .flat_map(|e| e.geometry.referenced_points())
        .collect();
    let mut seen = BTreeSet::new();
    picked
        .into_iter()
        .filter(|e| e.geometry.is_curve() || !of_curves.contains(&e.id))
        .map(|e| e.id)
        .filter(|id| seen.insert(*id))
        .collect()
}

/// Mirrors `sources` (curves and points, [`sources`]) about the line `axis`:
/// one batch adding their images and the mirror constraint about `axis`
/// (replacing the existing ones about it, whose pairs it keeps).
///
/// `curves` are the sketch's solved curves, for the preview. Sources that
/// are already mirrored about this line (as a source or an image), and
/// those that are their own image (a circle centred on the line, a line
/// along it or across it symmetrically, a point on it), are skipped;
/// `share` lets images share source points on the line (see the module
/// documentation; without it they get points of their own there, which the
/// editor falls back to when sharing would make a constraint redundant).
pub fn mirror(
    sketch: &SketchData,
    curves: &[(SketchId, ResolvedCurve)],
    sources: &[SketchId],
    axis: SketchId,
    share: bool,
) -> Result<Mirrored, String> {
    let line = Axis::of(sketch, axis)?;
    let existing = mirrors_about(sketch, axis);
    let old_pairs: Vec<(SketchId, SketchId)> = existing
        .iter()
        .flat_map(|(_, p)| p.iter().copied())
        .collect();
    let mirrored: BTreeSet<SketchId> = old_pairs.iter().flat_map(|&(s, i)| [s, i]).collect();
    // Point images so far, both ways (a reflection is its own inverse).
    let mut image_of: BTreeMap<SketchId, SketchId> = BTreeMap::new();
    for &(s, i) in &old_pairs {
        for (p, q) in point_pairs(sketch, s, i) {
            image_of.insert(p, q);
            image_of.insert(q, p);
        }
    }
    // Points that may be shared: on the line, and not held there by a
    // midpoint constraint (which also places them along it: sharing would
    // say the same twice). Point-on-line constraints the mirror replaces.
    let classes = coincidence_classes(sketch);
    let class = |p: SketchId| classes.get(&p).copied();
    let same_class =
        |p: SketchId, q: SketchId| p == q || (class(p).is_some() && class(p) == class(q));
    let midpoint_on_axis: Vec<SketchId> = sketch
        .constraints
        .iter()
        .filter_map(|c| match c.constraint {
            Constraint::Midpoint { point, line } if line == axis => Some(point),
            _ => None,
        })
        .collect();
    let shareable = |p: SketchId, at: DVec2| {
        share
            && line.distance(at) <= SAME_POINT
            && !midpoint_on_axis.iter().any(|&m| same_class(m, p))
    };

    let mut batch = Batch::default();
    let mut pairs = old_pairs;
    let mut shared: BTreeSet<SketchId> = BTreeSet::new();
    let mut preview_curves = Vec::new();
    let mut preview_points = Vec::new();
    let (mut made, mut skipped) = (0, 0);
    for &src in sources {
        let Some(e) = sketch.entity(src) else {
            continue;
        };
        let point_done =
            matches!(e.geometry, Geometry::Point { .. }) && image_of.contains_key(&src);
        if src == axis
            || mirrored.contains(&src)
            || point_done
            || own_image(sketch, &e.geometry, &line)?
        {
            skipped += 1;
            continue;
        }
        // Image points: earlier images (of a source, or the source of an
        // image: the reflection is its own inverse), shared points on the
        // line, or new points.
        let mut map = |p: SketchId| -> Result<SketchId, String> {
            if let Some(&q) = image_of.get(&p) {
                return Ok(q);
            }
            let at = point(sketch, p)?;
            let q = if shareable(p, at) {
                shared.insert(p);
                p
            } else {
                batch.point(line.reflect(at))
            };
            image_of.insert(p, q);
            image_of.insert(q, p);
            Ok(q)
        };
        let image = match &e.geometry {
            Geometry::Point { at } => {
                preview_points.push(line.reflect(*at));
                Geometry::Point {
                    at: line.reflect(*at),
                }
            }
            Geometry::Line { start, end } => Geometry::Line {
                start: map(*start)?,
                end: map(*end)?,
            },
            Geometry::Circle { center, radius } => Geometry::Circle {
                center: map(*center)?,
                radius: *radius,
            },
            Geometry::Arc { center, start, end } => {
                let c = map(*center)?;
                let (s, t) = (map(*end)?, map(*start)?);
                Geometry::Arc {
                    center: c,
                    start: s,
                    end: t,
                }
            }
            Geometry::Ellipse {
                center,
                major,
                minor,
            } => Geometry::Ellipse {
                center: map(*center)?,
                major: map(*major)?,
                minor: map(*minor)?,
            },
            Geometry::EllipticalArc {
                center,
                major,
                minor,
                start,
                end,
            } => {
                // As for arcs: the image runs from the image of the end.
                let (c, m, n) = (map(*center)?, map(*major)?, map(*minor)?);
                let (s, t) = (map(*end)?, map(*start)?);
                Geometry::EllipticalArc {
                    center: c,
                    major: m,
                    minor: n,
                    start: s,
                    end: t,
                }
            }
            Geometry::ControlSpline { control, degree } => Geometry::ControlSpline {
                control: control.iter().map(|p| map(*p)).collect::<Result<_, _>>()?,
                degree: *degree,
            },
            Geometry::InterpolatedSpline { through, params } => Geometry::InterpolatedSpline {
                through: through.iter().map(|p| map(*p)).collect::<Result<_, _>>()?,
                params: params.clone(),
            },
            Geometry::Conic {
                start,
                control,
                end,
                rho,
            } => Geometry::Conic {
                start: map(*start)?,
                control: map(*control)?,
                end: map(*end)?,
                rho: *rho,
            },
        };
        if let Some((_, c)) = curves.iter().find(|(id, _)| *id == src) {
            preview_curves.push(line.reflect_curve(c));
        }
        let id = batch.entity(image, e.construction);
        if matches!(e.geometry, Geometry::Point { .. }) {
            image_of.insert(src, id);
            image_of.insert(id, src);
        }
        pairs.push((src, id));
        made += 1;
    }
    if made == 0 {
        return Err(if skipped == 0 {
            "Nothing to mirror: select curves or points".into()
        } else {
            "Already mirrored about this line, or on it (their own mirror image)".into()
        });
    }
    // The mirror constraints about this line become one, and point-on-line
    // constraints of shared points go: the mirror holds them on the line.
    let mut gone: Vec<SketchId> = existing.iter().map(|(id, _)| *id).collect();
    gone.extend(
        sketch
            .constraints
            .iter()
            .filter_map(|c| match c.constraint {
                Constraint::PointOnCurve { point, curve }
                    if curve == axis && shared.iter().any(|&p| same_class(p, point)) =>
                {
                    Some(c.id)
                }
                _ => None,
            }),
    );
    if !gone.is_empty() {
        batch.edits.push(SketchEdit::Delete { items: gone });
    }
    batch.constraint(Constraint::Mirror { axis, pairs });
    Ok(Mirrored {
        plan: Plan {
            edits: batch.edits,
            removed: 0,
            preview: Preview::Shapes {
                curves: preview_curves,
                points: preview_points,
            },
        },
        made,
        skipped,
    })
}

/// A mirror ready to commit.
#[derive(Clone, Debug, PartialEq)]
pub struct Mirrored {
    /// The batch and its preview.
    pub plan: Plan,
    /// Images made.
    pub made: usize,
    /// Sources skipped: already mirrored about the line, or their own image.
    pub skipped: usize,
}

impl Mirrored {
    /// The status-bar message after committing.
    pub fn message(&self) -> String {
        let made = match self.made {
            1 => "Mirrored 1 item".to_string(),
            n => format!("Mirrored {n} items"),
        };
        match self.skipped {
            0 => made,
            n => format!("{made}; skipped {n} already mirrored about this line or lying on it"),
        }
    }
}

/// True when the reflection of `g` is `g` itself: every reflected point
/// lands on the corresponding point, or (for a line, an arc's ends, a spline
/// and a conic's ends) on the corresponding point in reverse order.
fn own_image(sketch: &SketchData, g: &Geometry<SketchId>, line: &Axis) -> Result<bool, String> {
    let at = |p: SketchId| point(sketch, p);
    let same = |a: DVec2, b: DVec2| (a - b).length() <= SAME_POINT;
    Ok(match g {
        Geometry::Point { at: p } => line.distance(*p) <= SAME_POINT,
        Geometry::Circle { center, .. } => line.distance(at(*center)?) <= SAME_POINT,
        Geometry::Line { start, end } => {
            let (a, b) = (at(*start)?, at(*end)?);
            (same(line.reflect(a), a) && same(line.reflect(b), b))
                || (same(line.reflect(a), b) && same(line.reflect(b), a))
        }
        Geometry::Arc { center, start, end } => {
            let (c, s, e) = (at(*center)?, at(*start)?, at(*end)?);
            same(line.reflect(c), c) && same(line.reflect(s), e)
        }
        Geometry::Ellipse {
            center,
            major,
            minor,
        } => {
            // Symmetric about the line: the centre on it and each axis end
            // reflected onto itself or onto the opposite vertex (the line is
            // an axis of the ellipse).
            let c = at(*center)?;
            let fixed = |q: DVec2| same(line.reflect(q), q) || same(line.reflect(q), 2.0 * c - q);
            same(line.reflect(c), c) && fixed(at(*major)?) && fixed(at(*minor)?)
        }
        Geometry::EllipticalArc {
            center,
            major,
            minor,
            start,
            end,
        } => {
            let c = at(*center)?;
            let fixed = |q: DVec2| same(line.reflect(q), q) || same(line.reflect(q), 2.0 * c - q);
            same(line.reflect(c), c)
                && fixed(at(*major)?)
                && fixed(at(*minor)?)
                && same(line.reflect(at(*start)?), at(*end)?)
        }
        Geometry::ControlSpline { control: pts, .. }
        | Geometry::InterpolatedSpline { through: pts, .. } => {
            let p: Vec<DVec2> = pts.iter().map(|k| at(*k)).collect::<Result<_, _>>()?;
            let r: Vec<DVec2> = p.iter().map(|q| line.reflect(*q)).collect();
            r.iter().zip(&p).all(|(a, b)| same(*a, *b))
                || r.iter().zip(p.iter().rev()).all(|(a, b)| same(*a, *b))
        }
        Geometry::Conic {
            start,
            control,
            end,
            ..
        } => {
            // Its control point on the line and its ends swapped (its three
            // points are never all on a line).
            let k = at(*control)?;
            same(line.reflect(k), k) && same(line.reflect(at(*start)?), at(*end)?)
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflections_are_involutions_and_keep_arcs_counter_clockwise() {
        let line = Axis {
            at: DVec2::new(1.0, 2.0),
            dir: DVec2::new(3.0, 4.0) / 5.0,
        };
        let p = DVec2::new(-2.0, 7.5);
        assert!((line.reflect(line.reflect(p)) - p).length() < 1e-12);
        assert!(line.distance(line.at + line.dir * 3.0) < 1e-12);
        let arc = ResolvedCurve::Arc {
            center: DVec2::new(4.0, 0.0),
            radius: 2.0,
            start_angle: 0.3,
            sweep_angle: 1.2,
        };
        let ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } = line.reflect_curve(&arc)
        else {
            unreachable!()
        };
        assert_eq!(radius, 2.0);
        assert_eq!(sweep_angle, 1.2);
        // The image starts at the image of the source's end.
        let src_end = DVec2::new(4.0, 0.0) + DVec2::from_angle(1.5) * 2.0;
        let img_start = center + DVec2::from_angle(start_angle) * 2.0;
        assert!((img_start - line.reflect(src_end)).length() < 1e-12);
    }
}
