//! Shape tools (design note 0007 §4): rectangles, slots, regular polygons
//! and tangent arcs, built from lines, arcs and constraints.
//!
//! Each builder computes the exact geometry, adds it to a [`Batch`] with
//! constraints the geometry already satisfies, and refuses degenerate input
//! with a message for the status bar. The `*_outline` functions give the
//! same geometry as polylines for the rubber band.

use std::f64::consts::TAU;

use glam::DVec2;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::SketchId;

use super::{perp, Batch, Ends, SAME_POINT};
use crate::sketch_mode::Snap;

/// Fewest sides of a regular polygon.
pub const MIN_SIDES: u32 = 3;
/// Most sides of a regular polygon.
pub const MAX_SIDES: u32 = 64;

/// Adds a closed chain of lines through `corners` (end to end, the last
/// back to the first), joined by coincident end points.
fn closed_chain(b: &mut Batch, corners: &[DVec2], construction: bool) -> Vec<Ends> {
    let n = corners.len();
    let lines: Vec<Ends> = (0..n)
        .map(|i| b.line(corners[i], corners[(i + 1) % n], construction))
        .collect();
    for i in 0..n {
        b.coincident(lines[i].end, lines[(i + 1) % n].start);
    }
    lines
}

/// Corners of the axis-aligned rectangle centred on `center` with a corner
/// at `corner`, in the order of [`centre_rectangle`]'s lines.
pub fn centre_rectangle_corners(center: DVec2, corner: DVec2) -> [DVec2; 4] {
    let q = 2.0 * center - corner;
    [
        corner,
        DVec2::new(q.x, corner.y),
        q,
        DVec2::new(corner.x, q.y),
    ]
}

/// A rectangle from its centre and a corner: four lines (horizontal,
/// vertical, horizontal, vertical), two construction diagonals, and the
/// centre point as the first diagonal's midpoint. DOF 4.
pub fn centre_rectangle(
    b: &mut Batch,
    center: &Snap,
    corner: &Snap,
    construction: bool,
) -> Result<(), String> {
    let c = centre_rectangle_corners(center.at, corner.at);
    if (c[0].x - c[2].x).abs() <= SAME_POINT || (c[0].y - c[2].y).abs() <= SAME_POINT {
        return Err("The rectangle has no width or no height".into());
    }
    let sides = closed_chain(b, &c, construction);
    for (i, s) in sides.iter().enumerate() {
        b.constraint(if i % 2 == 0 {
            Constraint::Horizontal { line: s.curve }
        } else {
            Constraint::Vertical { line: s.curve }
        });
    }
    b.attach(sides[0].start, corner);
    let d1 = b.line(c[0], c[2], true);
    let d2 = b.line(c[1], c[3], true);
    b.coincident(sides[0].start, d1.start);
    b.coincident(sides[2].start, d1.end);
    b.coincident(sides[1].start, d2.start);
    b.coincident(sides[3].start, d2.end);
    let m = b.snapped_point(center);
    b.constraint(Constraint::Midpoint {
        point: m,
        line: d1.curve,
    });
    Ok(())
}

/// Corners of the rectangle with one side from `a` to `b` and the opposite
/// side through `width` (projected onto the side's normal), or `None` when
/// it would be degenerate.
pub fn aligned_rectangle_corners(a: DVec2, b: DVec2, width: DVec2) -> Option<[DVec2; 4]> {
    let u = b - a;
    let len = u.length();
    if len <= SAME_POINT {
        return None;
    }
    let n = perp(u / len);
    let h = (width - a).dot(n);
    if h.abs() <= SAME_POINT {
        return None;
    }
    Some([a, b, b + n * h, a + n * h])
}

/// A rectangle at any angle from two corners along one side and a point on
/// the opposite side: perpendicular adjacent sides, parallel opposite ones.
/// DOF 5.
pub fn aligned_rectangle(
    b: &mut Batch,
    first: &Snap,
    second: &Snap,
    width: DVec2,
    construction: bool,
) -> Result<(), String> {
    let c = aligned_rectangle_corners(first.at, second.at, width)
        .ok_or("The rectangle has no width")?;
    let s = closed_chain(b, &c, construction);
    b.constraint(Constraint::Perpendicular {
        a: s[0].curve,
        b: s[1].curve,
    });
    b.constraint(Constraint::Parallel {
        a: s[0].curve,
        b: s[2].curve,
    });
    b.constraint(Constraint::Parallel {
        a: s[1].curve,
        b: s[3].curve,
    });
    b.attach(s[0].start, first);
    b.attach(s[0].end, second);
    Ok(())
}

/// A slot's geometry: the radius, and for the two half-circles their
/// centre, start and end (counter-clockwise).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlotGeometry {
    /// Half the slot's width.
    pub radius: f64,
    /// The first centre and its arc's start and end.
    pub first: [DVec2; 3],
    /// The second centre and its arc's start and end.
    pub second: [DVec2; 3],
}

/// The slot with arc centres `c0` and `c1` whose side passes through
/// `side` (its distance from the centre line is the radius).
pub fn slot_geometry(c0: DVec2, c1: DVec2, side: DVec2) -> Option<SlotGeometry> {
    let u = c1 - c0;
    let len = u.length();
    if len <= SAME_POINT {
        return None;
    }
    let n = perp(u / len);
    let r = (side - c0).dot(n).abs();
    if r <= SAME_POINT {
        return None;
    }
    Some(SlotGeometry {
        radius: r,
        // Around c0 from +n to −n through −u; around c1 from −n to +n
        // through +u.
        first: [c0, c0 + n * r, c0 - n * r],
        second: [c1, c1 - n * r, c1 + n * r],
    })
}

/// A slot: two half-circles about the ends of a construction centre line,
/// joined by two lines, tangent at all four joints, with equal radii.
/// DOF 5 (the two centres and the radius).
pub fn slot(
    b: &mut Batch,
    c0: &Snap,
    c1: &Snap,
    side: DVec2,
    construction: bool,
) -> Result<(), String> {
    let g = slot_geometry(c0.at, c1.at, side).ok_or("The slot has no length or no width")?;
    let axis = b.line(c0.at, c1.at, true);
    b.attach(axis.start, c0);
    b.attach(axis.end, c1);
    let a0 = b.arc(g.first[0], g.first[1], g.first[2], construction);
    let a1 = b.arc(g.second[0], g.second[1], g.second[2], construction);
    b.coincident(axis.start, a0.center);
    b.coincident(axis.end, a1.center);
    // Bottom: from the end of the first arc to the start of the second;
    // top: from the end of the second arc to the start of the first.
    let bottom = b.line(g.first[2], g.second[1], construction);
    let top = b.line(g.second[2], g.first[1], construction);
    b.coincident(a0.end, bottom.start);
    b.coincident(bottom.end, a1.start);
    b.coincident(a1.end, top.start);
    b.coincident(top.end, a0.start);
    b.tangent_at(a0.curve, bottom.curve, bottom.start);
    b.tangent_at(a1.curve, bottom.curve, bottom.end);
    b.tangent_at(a1.curve, top.curve, top.start);
    b.tangent_at(a0.curve, top.curve, top.end);
    b.constraint(Constraint::Equal {
        a: a0.curve,
        b: a1.curve,
    });
    Ok(())
}

/// The outline of a slot, as one closed polyline.
pub fn slot_outline(g: &SlotGeometry) -> Vec<DVec2> {
    let mut pts = arc_polyline(g.first[0], g.first[1], g.first[2]);
    pts.extend(arc_polyline(g.second[0], g.second[1], g.second[2]));
    pts.push(g.first[1]);
    pts
}

/// The vertices of the regular polygon about `center` with a vertex at
/// `vertex`, counter-clockwise.
pub fn polygon_vertices(center: DVec2, vertex: DVec2, sides: u32) -> Vec<DVec2> {
    let v = vertex - center;
    let (r, a0) = (v.length(), v.to_angle());
    (0..sides)
        .map(|k| {
            if k == 0 {
                vertex
            } else {
                center + DVec2::from_angle(a0 + TAU * k as f64 / sides as f64) * r
            }
        })
        .collect()
}

/// A regular polygon: lines joined end to end with every vertex on a
/// construction circle about the centre, and equal sides. DOF 4 (centre,
/// size, rotation).
pub fn polygon(
    b: &mut Batch,
    center: &Snap,
    vertex: &Snap,
    sides: u32,
    construction: bool,
) -> Result<(), String> {
    if !(MIN_SIDES..=MAX_SIDES).contains(&sides) {
        return Err(format!("A polygon has {MIN_SIDES} to {MAX_SIDES} sides"));
    }
    let radius = (vertex.at - center.at).length();
    if radius <= SAME_POINT {
        return Err("The polygon has no size".into());
    }
    let c = b.snapped_point(center);
    let circle = b.entity(Geometry::Circle { center: c, radius }, true);
    let lines = closed_chain(
        b,
        &polygon_vertices(center.at, vertex.at, sides),
        construction,
    );
    for l in &lines {
        b.constraint(Constraint::PointOnCurve {
            point: l.start,
            curve: circle,
        });
    }
    for l in &lines[1..] {
        b.constraint(Constraint::Equal {
            a: lines[0].curve,
            b: l.curve,
        });
    }
    b.attach(lines[0].start, vertex);
    Ok(())
}

/// Where a tangent arc starts: the curve it continues, the end point there,
/// its position and the unit direction leaving the curve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentStart {
    /// The curve continued.
    pub curve: SketchId,
    /// Its end point the arc starts at.
    pub point: SketchId,
    /// That point's position.
    pub at: DVec2,
    /// Unit direction leaving the curve there.
    pub direction: DVec2,
}

/// The arc leaving `p` in unit direction `d` and ending at `e`: its centre
/// and whether it turns left (counter-clockwise from `p`). `None` when `e`
/// is (nearly) on the tangent line, where the circle degenerates.
pub fn tangent_arc_geometry(p: DVec2, d: DVec2, e: DVec2) -> Option<(DVec2, bool)> {
    let v = e - p;
    let n = perp(d);
    let along = n.dot(v);
    // Refuse sagittas below a thousandth of the chord: the radius would
    // exceed 500 chords.
    if v.length() <= SAME_POINT || along.abs() <= 1e-3 * v.length() {
        return None;
    }
    let k = v.length_squared() / (2.0 * along);
    Some((p + n * k, k > 0.0))
}

/// A tangent arc: continues `start`'s curve tangentially to the clicked
/// end, coincident and tangent where it starts.
pub fn tangent_arc(
    b: &mut Batch,
    start: &TangentStart,
    start_snap: &Snap,
    end: &Snap,
    construction: bool,
) -> Result<(), String> {
    let (center, left) = tangent_arc_geometry(start.at, start.direction, end.at)
        .ok_or("The end is in line with the tangent: that would be a straight line")?;
    let c = b.point(center);
    let p = b.snapped_point(start_snap);
    let e = b.snapped_point(end);
    let (s0, s1) = if left { (p, e) } else { (e, p) };
    let arc = b.entity(
        Geometry::Arc {
            center: c,
            start: s0,
            end: s1,
        },
        construction,
    );
    b.tangent_at(start.curve, arc, p);
    Ok(())
}

/// The tangent arc's outline for the rubber band (a straight line when it
/// is degenerate).
pub fn tangent_arc_outline(start: &TangentStart, end: DVec2) -> Vec<DVec2> {
    match tangent_arc_geometry(start.at, start.direction, end) {
        Some((c, true)) => arc_polyline(c, start.at, end),
        Some((c, false)) => {
            let mut pts = arc_polyline(c, end, start.at);
            pts.reverse();
            pts
        }
        None => vec![start.at, end],
    }
}

/// A counter-clockwise arc about `c` from `a` to `b` as a polyline.
pub fn arc_polyline(c: DVec2, a: DVec2, b: DVec2) -> Vec<DVec2> {
    let r = (a - c).length();
    let a0 = (a - c).to_angle();
    let sweep = ((b - c).to_angle() - a0).rem_euclid(TAU);
    let sweep = if sweep == 0.0 { TAU } else { sweep };
    let n = ((sweep / 3.75f64.to_radians()).ceil() as usize).clamp(4, 96);
    (0..=n)
        .map(|i| c + DVec2::from_angle(a0 + sweep * i as f64 / n as f64) * r)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tangent_arc_through_a_quarter_turn() {
        // Leaving the origin along +x and ending at (1, 1): the centre is
        // (0, 1), turning left.
        let (c, left) = tangent_arc_geometry(DVec2::ZERO, DVec2::X, DVec2::new(1.0, 1.0)).unwrap();
        assert!((c - DVec2::new(0.0, 1.0)).length() < 1e-15);
        assert!(left);
        // Ending at (1, −1) turns right about (0, −1).
        let (c, left) = tangent_arc_geometry(DVec2::ZERO, DVec2::X, DVec2::new(1.0, -1.0)).unwrap();
        assert!((c - DVec2::new(0.0, -1.0)).length() < 1e-15);
        assert!(!left);
        // Straight ahead is refused.
        assert!(tangent_arc_geometry(DVec2::ZERO, DVec2::X, DVec2::new(2.0, 0.0)).is_none());
    }

    #[test]
    fn polygon_vertices_are_regular() {
        let v = polygon_vertices(DVec2::new(1.0, 2.0), DVec2::new(3.0, 2.0), 6);
        assert_eq!(v.len(), 6);
        for i in 0..6 {
            let side = (v[(i + 1) % 6] - v[i]).length();
            assert!((side - 2.0).abs() < 1e-12, "side {i} is {side}");
        }
    }

    #[test]
    fn slot_geometry_has_semicircles() {
        let g = slot_geometry(DVec2::ZERO, DVec2::new(4.0, 0.0), DVec2::new(1.0, 1.5)).unwrap();
        assert_eq!(g.radius, 1.5);
        assert_eq!(
            g.first,
            [DVec2::ZERO, DVec2::new(0.0, 1.5), DVec2::new(0.0, -1.5)]
        );
        assert_eq!(
            g.second,
            [
                DVec2::new(4.0, 0.0),
                DVec2::new(4.0, -1.5),
                DVec2::new(4.0, 1.5)
            ]
        );
    }
}
