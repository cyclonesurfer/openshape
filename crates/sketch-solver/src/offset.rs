//! The validity of an offset after a solve (design note 0007 §10.5).
//!
//! The offset's equations keep every image point on its offset curve and at
//! its joint. What equations cannot say is checked here on the solved
//! geometry: every image piece lies on its offset curve at both ends (at a
//! smooth joint only one piece's curve is an equation, so this fails when
//! the sources stop meeting tangentially), a line piece keeps its source's
//! direction and an arc piece its sense (a reversed piece means the distance
//! swallowed it), and no two pieces cross or touch except consecutive pieces
//! at their joint. An offset that fails is reported, never trimmed silently.

use std::f64::consts::{PI, TAU};

use glam::DVec2;

/// A solved curve of an offset, a source or an image.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Piece {
    /// A line from `a` to `b`.
    Line {
        /// Start.
        a: DVec2,
        /// End.
        b: DVec2,
    },
    /// The counter-clockwise arc about `c` from `a` to `b` (radius `|a − c|`).
    Arc {
        /// Centre.
        c: DVec2,
        /// Start.
        a: DVec2,
        /// End.
        b: DVec2,
    },
    /// A full circle.
    Circle {
        /// Centre.
        c: DVec2,
        /// Radius.
        r: f64,
    },
}

/// One piece of an offset: its source, its image, and which side.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct OffsetPiece {
    /// The source curve.
    pub src: Piece,
    /// The image curve.
    pub img: Piece,
    /// The image lies to the source's left (towards an arc's centre,
    /// inside a circle).
    pub left: bool,
    /// The image's end point indices (start, end); none for a circle.
    pub ends: Option<(usize, usize)>,
    /// The image's start and end are held on the normal at the source's
    /// (a free end or a smooth joint), so lie beside it.
    pub held: (bool, bool),
}

/// Why an offset is invalid, naming the piece (its position in the
/// constraint) or the two pieces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Invalid {
    /// The image no longer lies on its offset curve at an end (`true` for
    /// the end, `false` for the start).
    OffCurve(usize, bool),
    /// A line piece points the other way, an arc piece turned back or a
    /// circle shrank to nothing: the distance is larger than the piece
    /// allows.
    TurnedBack(usize),
    /// Two pieces cross, touch or overlap away from a joint they share.
    Cross(usize, usize),
}

/// Angle of `v`.
fn angle(v: DVec2) -> f64 {
    v.y.atan2(v.x)
}

/// `a` wrapped into (−π, π].
fn wrap(a: f64) -> f64 {
    let w = a.rem_euclid(TAU);
    if w > PI {
        w - TAU
    } else {
        w
    }
}

/// The counter-clockwise sweep from direction `a` to direction `b`, in
/// (0, 2π].
fn sweep(a: DVec2, b: DVec2) -> f64 {
    let s = (angle(b) - angle(a)).rem_euclid(TAU);
    if s == 0.0 {
        TAU
    } else {
        s
    }
}

/// The signed distance from `p` to the offset curve of `src` (positive
/// beyond it, away from the source on the offset side).
fn off_curve(src: &Piece, left: bool, d: f64, p: DVec2) -> f64 {
    let side = if left { 1.0 } else { -1.0 };
    match *src {
        Piece::Line { a, b } => {
            let u = (b - a).normalize_or_zero();
            u.perp_dot(p - a) * side - d
        }
        Piece::Arc { c, a, .. } => {
            // Left of a counter-clockwise curve is towards its centre.
            let target = (a - c).length() - side * d;
            (p - c).length() - target
        }
        Piece::Circle { c, r } => (p - c).length() - (r - side * d),
    }
}

/// The ends of a line or arc.
fn ends(p: &Piece) -> Option<(DVec2, DVec2)> {
    match *p {
        Piece::Line { a, b } | Piece::Arc { a, b, .. } => Some((a, b)),
        Piece::Circle { .. } => None,
    }
}

/// True if `x`, a point of the piece's carrier (to within `tol`), lies on
/// the piece, ends included, to within `tol`.
fn within(p: &Piece, x: DVec2, tol: f64) -> bool {
    match *p {
        Piece::Line { a, b } => {
            let len = (b - a).length();
            if len == 0.0 {
                return (x - a).length() <= tol;
            }
            let t = (x - a).dot(b - a) / len;
            t >= -tol && t <= len + tol
        }
        Piece::Arc { c, a, b } => {
            if (x - a).length() <= tol || (x - b).length() <= tol {
                return true;
            }
            sweep(a - c, x - c) <= sweep(a - c, b - c)
        }
        Piece::Circle { .. } => true,
    }
}

/// The distance from `x` to the piece.
fn distance(p: &Piece, x: DVec2) -> f64 {
    match *p {
        Piece::Line { a, b } => {
            let v = b - a;
            let t = if v.length_squared() == 0.0 {
                0.0
            } else {
                ((x - a).dot(v) / v.length_squared()).clamp(0.0, 1.0)
            };
            (x - (a + v * t)).length()
        }
        Piece::Arc { c, a, b } => {
            let r = (a - c).length();
            let to = x - c;
            let on_circle = if to.length() == 0.0 {
                a
            } else {
                c + to * (r / to.length())
            };
            if within(p, on_circle, 0.0) {
                (x - on_circle).length()
            } else {
                (x - a).length().min((x - b).length())
            }
        }
        Piece::Circle { c, r } => ((x - c).length() - r).abs(),
    }
}

/// A line or circle carrying a piece.
#[derive(Clone, Copy, Debug)]
enum Carrier {
    Line { p: DVec2, u: DVec2 },
    Circle { c: DVec2, r: f64 },
}

fn carrier(p: &Piece) -> Option<Carrier> {
    match *p {
        Piece::Line { a, b } => {
            let u = (b - a).normalize_or_zero();
            (u != DVec2::ZERO).then_some(Carrier::Line { p: a, u })
        }
        Piece::Arc { c, a, .. } => Some(Carrier::Circle {
            c,
            r: (a - c).length(),
        }),
        Piece::Circle { c, r } => Some(Carrier::Circle { c, r }),
    }
}

/// Where two carriers meet, to within `tol` (one point where they touch);
/// `None` when they coincide (to within `tol`).
fn meet(a: Carrier, b: Carrier, tol: f64) -> Option<Vec<DVec2>> {
    Some(match (a, b) {
        (Carrier::Line { p, u }, Carrier::Line { p: q, u: v }) => {
            let den = u.perp_dot(v);
            if den.abs() <= 1e-12 {
                if u.perp_dot(q - p).abs() <= tol {
                    return None;
                }
                Vec::new()
            } else {
                vec![p + u * ((q - p).perp_dot(v) / den)]
            }
        }
        (Carrier::Line { p, u }, Carrier::Circle { c, r })
        | (Carrier::Circle { c, r }, Carrier::Line { p, u }) => {
            let foot = p + u * (c - p).dot(u);
            let h = (foot - c).length();
            if h > r + tol {
                Vec::new()
            } else if h >= r - tol {
                vec![foot]
            } else {
                let s = (r * r - h * h).sqrt();
                vec![foot - u * s, foot + u * s]
            }
        }
        (Carrier::Circle { c, r }, Carrier::Circle { c: k, r: s }) => {
            let v = k - c;
            let d = v.length();
            if d <= tol {
                if (r - s).abs() <= tol {
                    return None;
                }
                return Some(Vec::new());
            }
            let n = v / d;
            if d > r + s + tol || d < (r - s).abs() - tol {
                Vec::new()
            } else if d >= r + s - tol {
                vec![c + n * r]
            } else if d <= (r - s).abs() + tol {
                // Internal touching: on the side of the larger circle's
                // far point.
                vec![if r >= s { c + n * r } else { c - n * r }]
            } else {
                let x = (r * r - s * s + d * d) / (2.0 * d);
                let y = (r * r - x * x).max(0.0).sqrt();
                let m = c + n * x;
                let t = DVec2::new(-n.y, n.x);
                vec![m + t * y, m - t * y]
            }
        }
    })
}

/// The points where pieces `p` and `q` cross or touch, to within `tol`,
/// and the ends of their overlap if they share a carrier.
fn contacts(p: &Piece, q: &Piece, tol: f64) -> Vec<DVec2> {
    let mut out = Vec::new();
    // Ends of one piece on the other: touching, and overlap ends.
    for (x, y) in [(p, q), (q, p)] {
        if let Some((a, b)) = ends(x) {
            for e in [a, b] {
                if distance(y, e) <= tol {
                    out.push(e);
                }
            }
        }
    }
    // A collapsed line has no carrier: its point was handled above.
    if let (Some(cp), Some(cq)) = (carrier(p), carrier(q)) {
        match meet(cp, cq, tol) {
            Some(points) => out.extend(
                points
                    .into_iter()
                    .filter(|x| within(p, *x, tol) && within(q, *x, tol)),
            ),
            None => {
                // The same carrier: two full circles, or a full circle and
                // an arc of it, overlap entirely.
                if let (Piece::Circle { c, r }, other) | (other, Piece::Circle { c, r }) = (p, q) {
                    out.push(ends(other).map_or(*c + DVec2::X * *r, |(a, _)| a));
                }
            }
        }
    }
    out
}

/// The validity check of design note 0007 §10.5 for the pieces of one
/// offset at distance `d`, to within `tol` (metres). Pieces that share an
/// image end point are consecutive and may meet there.
pub(crate) fn check(pieces: &[OffsetPiece], d: f64, tol: f64) -> Result<(), Invalid> {
    for (i, p) in pieces.iter().enumerate() {
        match (p.src, p.img) {
            (Piece::Line { a, b }, Piece::Line { a: ia, b: ib }) => {
                for (end, x) in [(false, ia), (true, ib)] {
                    if off_curve(&p.src, p.left, d, x).abs() > tol {
                        return Err(Invalid::OffCurve(i, end));
                    }
                }
                if (ib - ia).dot(b - a) <= 0.0 {
                    return Err(Invalid::TurnedBack(i));
                }
            }
            (
                Piece::Arc { c, a, b },
                Piece::Arc {
                    c: ic,
                    a: ia,
                    b: ib,
                },
            ) => {
                for (end, x) in [(false, ia), (true, ib)] {
                    if off_curve(&p.src, p.left, d, x).abs() > tol {
                        return Err(Invalid::OffCurve(i, end));
                    }
                }
                // An end held on the normal through the centre must be on
                // the source end's side of it: the normal also crosses the
                // offset circle opposite.
                for (h, x, y) in [(p.held.0, ia, a), (p.held.1, ib, b)] {
                    if h && (x - ic).dot(y - c) <= 0.0 {
                        return Err(Invalid::TurnedBack(i));
                    }
                }
                // The image's ends, measured from the source's: trimmed or
                // extended by less than half a turn at each end, its sweep
                // is the source's less the start's advance plus the end's.
                let along = sweep(a - c, b - c) - wrap(angle(ia - ic) - angle(a - c))
                    + wrap(angle(ib - ic) - angle(b - c));
                let r = (ia - ic).length().max(tol);
                if along <= 0.0 || (sweep(ia - ic, ib - ic) - along).abs() * r > tol.max(1e-9 * r) {
                    return Err(Invalid::TurnedBack(i));
                }
            }
            (Piece::Circle { .. }, Piece::Circle { c, r }) => {
                // The distance swallowed the circle.
                if r <= tol {
                    return Err(Invalid::TurnedBack(i));
                }
                if off_curve(&p.src, p.left, d, c + DVec2::X * r).abs() > tol {
                    return Err(Invalid::OffCurve(i, false));
                }
            }
            _ => return Err(Invalid::TurnedBack(i)),
        }
    }
    for i in 0..pieces.len() {
        for j in i + 1..pieces.len() {
            let (p, q) = (&pieces[i], &pieces[j]);
            // Joints the two share: contacts there are the joint.
            let shared: Vec<DVec2> = match (p.ends, q.ends, ends(&p.img)) {
                (Some((ps, pe)), Some((qs, qe)), Some((a, b))) => [(ps, a), (pe, b)]
                    .into_iter()
                    .filter(|(k, _)| *k == qs || *k == qe)
                    .map(|(_, x)| x)
                    .collect(),
                _ => Vec::new(),
            };
            let near = 10.0 * tol;
            if contacts(&p.img, &q.img, tol)
                .into_iter()
                .any(|x| shared.iter().all(|s| (x - *s).length() > near))
            {
                return Err(Invalid::Cross(i, j));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64) -> DVec2 {
        DVec2::new(x, y)
    }

    fn line(a: DVec2, b: DVec2) -> Piece {
        Piece::Line { a, b }
    }

    #[test]
    fn crossing_lines_and_touching_ends_are_contacts() {
        let tol = 1e-9;
        let a = line(v(0.0, 0.0), v(2.0, 2.0));
        let b = line(v(0.0, 2.0), v(2.0, 0.0));
        assert!(contacts(&a, &b, tol)
            .iter()
            .any(|x| (*x - v(1.0, 1.0)).length() < 1e-12));
        // An end touching the middle of another line.
        let t = line(v(1.0, 0.0), v(1.0, 1.0));
        let u = line(v(0.0, 1.0), v(2.0, 1.0));
        assert!(!contacts(&t, &u, tol).is_empty());
        // Apart.
        let w = line(v(0.0, 3.0), v(2.0, 3.0));
        assert!(contacts(&u, &w, tol).is_empty());
        // Collinear overlap.
        let o = line(v(1.5, 1.0), v(3.0, 1.0));
        assert!(contacts(&u, &o, tol).len() >= 2);
    }

    #[test]
    fn arcs_meet_lines_only_within_their_sweep() {
        let tol = 1e-9;
        // The upper half of the unit circle.
        let arc = Piece::Arc {
            c: v(0.0, 0.0),
            a: v(1.0, 0.0),
            b: v(-1.0, 0.0),
        };
        let above = line(v(0.0, 0.0), v(0.0, 2.0));
        let below = line(v(0.0, 0.0), v(0.0, -2.0));
        assert_eq!(contacts(&arc, &above, tol).len(), 1);
        assert!(contacts(&arc, &below, tol).is_empty());
        // A tangent line touches.
        let tangent = line(v(-1.0, 1.0), v(1.0, 1.0));
        assert_eq!(contacts(&arc, &tangent, tol).len(), 1);
    }

    #[test]
    fn a_square_offset_outwards_is_valid_and_too_far_inwards_is_not() {
        // The unit square, counter-clockwise; outwards is to the right.
        let sq = [v(0.0, 0.0), v(1.0, 0.0), v(1.0, 1.0), v(0.0, 1.0)];
        let offset = |d: f64, left: bool| -> Vec<OffsetPiece> {
            let s = if left { d } else { -d };
            let corners = [v(s, s), v(1.0 - s, s), v(1.0 - s, 1.0 - s), v(s, 1.0 - s)];
            (0..4)
                .map(|k| OffsetPiece {
                    src: line(sq[k], sq[(k + 1) % 4]),
                    img: line(corners[k], corners[(k + 1) % 4]),
                    left,
                    ends: Some((k, (k + 1) % 4)),
                    held: (false, false),
                })
                .collect()
        };
        assert_eq!(check(&offset(0.25, false), 0.25, 1e-9), Ok(()));
        assert_eq!(check(&offset(0.25, true), 0.25, 1e-9), Ok(()));
        // Inwards by more than half the side: every piece turns back.
        assert_eq!(
            check(&offset(0.75, true), 0.75, 1e-9),
            Err(Invalid::TurnedBack(0))
        );
        // A piece off its curve.
        let mut p = offset(0.25, false);
        p[1].img = line(v(1.25, -0.25), v(1.3, 1.25));
        assert_eq!(check(&p, 0.25, 1e-9), Err(Invalid::OffCurve(1, true)));
    }

    #[test]
    fn pieces_that_cross_away_from_a_shared_joint_are_refused() {
        // Two crossing lines offset to their left by 0.5: their images
        // cross at (1.5, 0.5).
        let src = [
            line(v(0.0, 0.0), v(4.0, 0.0)),
            line(v(2.0, -2.0), v(2.0, 2.0)),
        ];
        let img = [
            line(v(0.0, 0.5), v(4.0, 0.5)),
            line(v(1.5, -2.0), v(1.5, 2.0)),
        ];
        let mut pieces: Vec<OffsetPiece> = (0..2)
            .map(|k| OffsetPiece {
                src: src[k],
                img: img[k],
                left: true,
                ends: Some((2 * k, 2 * k + 1)),
                held: (true, true),
            })
            .collect();
        assert_eq!(check(&pieces, 0.5, 1e-9), Err(Invalid::Cross(0, 1)));
        // Sharing an end point elsewhere does not excuse a crossing.
        pieces[1].ends = Some((1, 3));
        assert_eq!(check(&pieces, 0.5, 1e-9), Err(Invalid::Cross(0, 1)));
        // Consecutive pieces meeting only at their joint are fine: an L.
        let l = [
            OffsetPiece {
                src: line(v(0.0, 0.0), v(4.0, 0.0)),
                img: line(v(0.0, 0.5), v(3.5, 0.5)),
                left: true,
                ends: Some((0, 1)),
                held: (true, false),
            },
            OffsetPiece {
                src: line(v(4.0, 0.0), v(4.0, 4.0)),
                img: line(v(3.5, 0.5), v(3.5, 4.0)),
                left: true,
                ends: Some((1, 2)),
                held: (false, true),
            },
        ];
        assert_eq!(check(&l, 0.5, 1e-9), Ok(()));
    }

    #[test]
    fn an_arc_swallowed_by_the_distance_turns_back() {
        // A quarter arc of radius 1 offset outwards stays valid; an image
        // whose ends have swapped order is refused.
        let src = Piece::Arc {
            c: v(0.0, 0.0),
            a: v(1.0, 0.0),
            b: v(0.0, 1.0),
        };
        let good = OffsetPiece {
            src,
            img: Piece::Arc {
                c: v(0.0, 0.0),
                a: v(2.0, 0.0),
                b: v(0.0, 2.0),
            },
            left: false,
            ends: Some((0, 1)),
            held: (true, true),
        };
        assert_eq!(check(&[good], 1.0, 1e-9), Ok(()));
        let a = 50f64.to_radians();
        let b = 40f64.to_radians();
        let bad = OffsetPiece {
            img: Piece::Arc {
                c: v(0.0, 0.0),
                a: DVec2::from_angle(a) * 2.0,
                b: DVec2::from_angle(b) * 2.0,
            },
            ..good
        };
        assert_eq!(check(&[bad], 1.0, 1e-9), Err(Invalid::TurnedBack(0)));
        // An end held on the normal but across the centre: the normal's
        // other crossing of the offset circle.
        let flipped = OffsetPiece {
            img: Piece::Arc {
                c: v(0.0, 0.0),
                a: v(2.0, 0.0),
                b: v(0.0, -2.0),
            },
            ..good
        };
        assert_eq!(check(&[flipped], 1.0, 1e-9), Err(Invalid::TurnedBack(0)));
    }
}
