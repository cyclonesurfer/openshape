//! Shared helpers for the solver's integration tests: a sketch builder and
//! an independent checker that measures every constraint's violation from
//! the geometry with its own formulas (not the solver's equations).

#![allow(dead_code)]

pub mod gen;

use openshape_sketch_solver::{
    resolve_curves, Constraint, ConstraintEntry, DVec2, DistanceKind, Entity, Geometry,
    ResolvedCurve, Sketch, TangentKind,
};

/// Builds sketches with `u32` keys; entities and constraints share one
/// counter so keys never collide.
#[derive(Default, Clone)]
pub struct B {
    pub sketch: Sketch<u32>,
    next: u32,
}

impl B {
    pub fn new() -> Self {
        Self::default()
    }

    fn key(&mut self) -> u32 {
        self.next += 1;
        self.next
    }

    pub fn entity(&mut self, geometry: Geometry<u32>) -> u32 {
        let id = self.key();
        self.sketch.entities.push(Entity {
            id,
            geometry,
            construction: false,
        });
        id
    }

    pub fn point(&mut self, x: f64, y: f64) -> u32 {
        self.entity(Geometry::Point {
            at: DVec2::new(x, y),
        })
    }

    /// A line with its own two end points: (line, start, end).
    pub fn line(&mut self, x0: f64, y0: f64, x1: f64, y1: f64) -> (u32, u32, u32) {
        let s = self.point(x0, y0);
        let e = self.point(x1, y1);
        let l = self.entity(Geometry::Line { start: s, end: e });
        (l, s, e)
    }

    pub fn line_between(&mut self, s: u32, e: u32) -> u32 {
        self.entity(Geometry::Line { start: s, end: e })
    }

    /// A circle with its own centre: (circle, centre).
    pub fn circle(&mut self, x: f64, y: f64, r: f64) -> (u32, u32) {
        let c = self.point(x, y);
        let k = self.entity(Geometry::Circle {
            center: c,
            radius: r,
        });
        (k, c)
    }

    /// An arc with its own points: (arc, centre, start, end).
    pub fn arc(&mut self, c: DVec2, r: f64, a0: f64, a1: f64) -> (u32, u32, u32, u32) {
        let cc = self.point(c.x, c.y);
        let s = self.point(c.x + r * a0.cos(), c.y + r * a0.sin());
        let e = self.point(c.x + r * a1.cos(), c.y + r * a1.sin());
        let a = self.entity(Geometry::Arc {
            center: cc,
            start: s,
            end: e,
        });
        (a, cc, s, e)
    }

    pub fn c(&mut self, constraint: Constraint<u32>) -> u32 {
        let id = self.key();
        self.sketch
            .constraints
            .push(ConstraintEntry { id, constraint });
        id
    }

    pub fn coincident(&mut self, a: u32, b: u32) -> u32 {
        self.c(Constraint::Coincident { a, b })
    }

    pub fn fix(&mut self, p: u32) -> u32 {
        let at = self.sketch.point(p).unwrap();
        self.c(Constraint::Fix { point: p, at })
    }

    pub fn build(&self) -> Sketch<u32> {
        self.sketch.clone()
    }
}

pub fn pos(s: &Sketch<u32>, p: u32) -> DVec2 {
    s.point(p).expect("point")
}

/// (centre, radius) of a circle or arc.
pub fn circ(s: &Sketch<u32>, k: u32) -> (DVec2, f64) {
    match s.entity(k).unwrap().geometry {
        Geometry::Circle { center, radius } => (pos(s, center), radius),
        Geometry::Arc { center, start, .. } => {
            let c = pos(s, center);
            (c, (pos(s, start) - c).length())
        }
        _ => panic!("not circular"),
    }
}

/// (centre, first semi-axis, second semi-axis) of an ellipse or
/// elliptical arc.
pub fn ell(s: &Sketch<u32>, k: u32) -> (DVec2, DVec2, DVec2) {
    match s.entity(k).unwrap().geometry {
        Geometry::Ellipse {
            center,
            major,
            minor,
        }
        | Geometry::EllipticalArc {
            center,
            major,
            minor,
            ..
        } => {
            let c = pos(s, center);
            (c, pos(s, major) - c, pos(s, minor) - c)
        }
        _ => panic!("not an ellipse"),
    }
}

/// `q` in the ellipse's frame scaled to the unit circle: `(x/a, y/b)` with
/// `x`, `y` the coordinates along the (perpendicular) axes.
fn unit_frame(s: &Sketch<u32>, k: u32, q: DVec2) -> DVec2 {
    let (c, a, b) = ell(s, k);
    let w = q - c;
    DVec2::new(w.dot(a) / a.length_squared(), w.dot(b) / b.length_squared())
}

/// Distance of `q` from an ellipse to first order (Sampson's distance
/// `|F|/‖∇F‖` of `F = x²/a² + y²/b² − 1` in the axis frame): exact enough
/// for violations near zero, independent of the solver's parameter.
pub fn dist_to_ellipse(s: &Sketch<u32>, k: u32, q: DVec2) -> f64 {
    let (c, a, b) = ell(s, k);
    let (la, lb) = (a.length(), b.length());
    let w = q - c;
    let (x, y) = (w.dot(a) / la, w.dot(b) / lb);
    let f = x * x / (la * la) + y * y / (lb * lb) - 1.0;
    let g = DVec2::new(2.0 * x / (la * la), 2.0 * y / (lb * lb)).length();
    f.abs() / g
}

/// (start, control, end, ρ) of a conic.
pub fn cnc(s: &Sketch<u32>, k: u32) -> ([DVec2; 3], f64) {
    match s.entity(k).unwrap().geometry {
        Geometry::Conic {
            start,
            control,
            end,
            rho,
        } => ([pos(s, start), pos(s, control), pos(s, end)], rho),
        _ => panic!("not a conic"),
    }
}

/// Distance of `q` from a conic arc, from Lee's implicit equation in the
/// barycentric coordinates `τ` of the control triangle, `F = τ₁² −
/// 4w²τ₀τ₂ = 0` with `w = ρ/(1 − ρ)` (independent of the solver's
/// parametric form): Sampson's first-order distance `|F|/‖∇F‖`, plus how
/// far `q` lies outside the triangle (the arc is the part of the conic
/// with every `τ ≥ 0`).
pub fn dist_to_conic(s: &Sketch<u32>, k: u32, q: DVec2) -> f64 {
    let ([p0, p1, p2], rho) = cnc(s, k);
    let w = rho / (1.0 - rho);
    let area = (p1 - p0).perp_dot(p2 - p0);
    // τᵢ(x) = cross(Pⱼ − x, Pₖ − x)/area, affine in x with gradient
    // (Pⱼ − Pₖ)ᵖᵉʳᵖ/area.
    let tau = |a: DVec2, b: DVec2| (a - q).perp_dot(b - q) / area;
    let grad = |a: DVec2, b: DVec2| (a - b).perp() / area;
    let (t0, t1, t2) = (tau(p1, p2), tau(p2, p0), tau(p0, p1));
    let (g0, g1, g2) = (grad(p1, p2), grad(p2, p0), grad(p0, p1));
    let f = t1 * t1 - 4.0 * w * w * t0 * t2;
    let gf = g1 * (2.0 * t1) - (g0 * t2 + g2 * t0) * (4.0 * w * w);
    let size = (p1 - p0).length().max((p2 - p1).length());
    let outside = [t0, t1, t2]
        .iter()
        .map(|t| (-t).max(0.0))
        .fold(0.0, f64::max)
        * size;
    f.abs() / gf.length() + outside
}

/// (start, end) of a line.
pub fn line(s: &Sketch<u32>, k: u32) -> (DVec2, DVec2) {
    match s.entity(k).unwrap().geometry {
        Geometry::Line { start, end } => (pos(s, start), pos(s, end)),
        _ => panic!("not a line"),
    }
}

fn wrap(a: f64) -> f64 {
    let mut a = a % std::f64::consts::TAU;
    if a > std::f64::consts::PI {
        a -= std::f64::consts::TAU;
    } else if a <= -std::f64::consts::PI {
        a += std::f64::consts::TAU;
    }
    a
}

fn dist_to_line(p: DVec2, (s, e): (DVec2, DVec2)) -> f64 {
    let d = e - s;
    d.perp_dot(p - s).abs() / d.length()
}

/// Evaluates a clamped B-spline (de Boor's algorithm, P&T A3.1 in its
/// recursive form) — written independently of the solver's basis code.
pub fn de_boor(degree: usize, knots: &[f64], pts: &[DVec2], t: f64) -> DVec2 {
    de_boor_in(degree, knots, pts, t, f64::NEG_INFINITY)
}

/// [`de_boor`] using the polynomial piece of the knot span that starts at
/// `span_start` when `t` lies in it or at its end (the span is found from
/// `span_start` if that is a knot, otherwise from `t`).
pub fn de_boor_in(degree: usize, knots: &[f64], pts: &[DVec2], t: f64, span_start: f64) -> DVec2 {
    let n = pts.len() - 1;
    let t = t.clamp(knots[degree], knots[n + 1]);
    let mut k = degree;
    if span_start > f64::NEG_INFINITY {
        while k < n && knots[k + 1] <= span_start {
            k += 1;
        }
    } else {
        while k < n && t >= knots[k + 1] {
            k += 1;
        }
    }
    let mut d: Vec<DVec2> = (0..=degree).map(|j| pts[j + k - degree]).collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = j + k - degree;
            let den = knots[i + degree + 1 - r] - knots[i];
            let alpha = if den > 0.0 { (t - knots[i]) / den } else { 0.0 };
            d[j] = d[j - 1] * (1.0 - alpha) + d[j] * alpha;
        }
    }
    d[degree]
}

/// The explicit geometry of curve `k` alone (other curves may be
/// degenerate without affecting it).
pub fn resolved(s: &Sketch<u32>, k: u32) -> ResolvedCurve {
    let e = s.entity(k).expect("curve").clone();
    let mut own = Sketch::default();
    for p in e.geometry.referenced_points() {
        if own.entity(p).is_none() {
            own.entities.push(s.entity(p).expect("point").clone());
        }
    }
    own.entities.push(e);
    resolve_curves(&own)
        .expect("resolvable")
        .into_iter()
        .find(|(id, _, _)| *id == k)
        .expect("curve")
        .2
}

/// Distance from `p` to a spline curve. Each knot span (a polynomial
/// piece) is sampled densely, and golden-section search refines every
/// local minimum of the samples within the piece (a looping curve can pass
/// near `p` more than once; at a corner the two pieces are separate).
pub fn dist_to_spline(s: &Sketch<u32>, k: u32, p: DVec2) -> f64 {
    let ResolvedCurve::Spline {
        degree,
        knots,
        points,
    } = resolved(s, k)
    else {
        panic!("not a spline")
    };
    let mut best = f64::INFINITY;
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for span in knots.windows(2) {
        let (lo, hi) = (span[0], span[1]);
        if hi <= lo {
            continue;
        }
        // Evaluate inside the span, approaching its ends from within (the
        // ends themselves are evaluated by the neighbouring spans too).
        let f = |t: f64| (de_boor_in(degree, &knots, &points, t, lo) - p).length();
        let n = 4000;
        let ts: Vec<f64> = (0..=n)
            .map(|i| lo + (hi - lo) * i as f64 / n as f64)
            .collect();
        let d: Vec<f64> = ts.iter().map(|&t| f(t)).collect();
        for i in 0..=n {
            let left = if i > 0 { d[i - 1] } else { f64::INFINITY };
            let right = if i < n { d[i + 1] } else { f64::INFINITY };
            best = best.min(d[i]);
            if !(d[i] <= left && d[i] <= right) {
                continue;
            }
            let (mut a, mut b) = (ts[i.saturating_sub(1)], ts[(i + 1).min(n)]);
            for _ in 0..200 {
                let c = b - g * (b - a);
                let e = a + g * (b - a);
                if f(c) < f(e) {
                    b = e;
                } else {
                    a = c;
                }
            }
            best = best.min(f(0.5 * (a + b)));
        }
    }
    best
}

/// Unit tangent direction of a curve at its end `key` (lines: direction;
/// arcs: perpendicular to the radius there; splines: the first or last
/// control leg — the end that is `key` itself, or else the nearer one).
fn end_tangent(s: &Sketch<u32>, k: u32, key: u32) -> DVec2 {
    let p = pos(s, key);
    match s.entity(k).unwrap().geometry.clone() {
        Geometry::Line { .. } => {
            let (a, b) = line(s, k);
            (b - a).normalize()
        }
        Geometry::Arc { center, .. } => (p - pos(s, center)).perp().normalize(),
        Geometry::EllipticalArc { .. } => {
            // Perpendicular to the implicit equation's gradient
            // (x/a², y/b²) in the axis frame.
            let (_, a, b) = ell(s, k);
            let u = unit_frame(s, k, p);
            let grad = a.normalize() * (u.x / a.length()) + b.normalize() * (u.y / b.length());
            grad.perp().normalize()
        }
        Geometry::Conic { start, end, .. } => {
            // Along the leg from the end to the control point (the end
            // that is the point or joined to it; else the nearer one).
            let ([p0, p1, p2], _) = cnc(s, k);
            let joined = |q: u32| coincident_class(s, key).contains(&q);
            let at_start = match (joined(start), joined(end)) {
                (true, false) => true,
                (false, true) => false,
                _ => (p0 - p).length() <= (p2 - p).length(),
            };
            if at_start {
                (p1 - p0).normalize()
            } else {
                (p2 - p1).normalize()
            }
        }
        Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. } => {
            let ResolvedCurve::Spline { points, .. } = resolved(s, k) else {
                unreachable!()
            };
            let n = points.len();
            let refs = s.entity(k).unwrap().geometry.referenced_points();
            // The end that is the point, or is joined to it by coincident
            // constraints (as the solver decides); else the nearer end.
            let joined = |q: u32| coincident_class(s, key).contains(&q);
            let at_start = match (joined(refs[0]), joined(refs[refs.len() - 1])) {
                (true, false) => true,
                (false, true) => false,
                _ => (points[0] - p).length() <= (points[n - 1] - p).length(),
            };
            if at_start {
                (points[1] - points[0]).normalize()
            } else {
                (points[n - 1] - points[n - 2]).normalize()
            }
        }
        _ => panic!("no tangent"),
    }
}

/// Which end of curve `k` is joined to point `key`: true for its start
/// (the end that is `key` or joined to it by coincident constraints; else
/// the nearer end).
fn joined_at_start(s: &Sketch<u32>, key: u32, first: u32, last: u32) -> bool {
    let joined = |q: u32| coincident_class(s, key).contains(&q);
    match (joined(first), joined(last)) {
        (true, false) => true,
        (false, true) => false,
        _ => (pos(s, first) - pos(s, key)).length() <= (pos(s, last) - pos(s, key)).length(),
    }
}

/// The unit tangent pointing away from the joint at point `key` along curve
/// `k`, and the signed curvature for travel that way (positive turning
/// left), each from its own closed form (design note 0010 §6): lines 0;
/// arcs ±1/r; elliptical arcs `|A×B|/‖E′(θ)‖³` at the point's θ; conics
/// from the quotient rule on the rational form; splines from the
/// derivative control points of Piegl & Tiller eq. 3.8 applied twice.
pub fn away_curvature(s: &Sketch<u32>, k: u32, key: u32) -> (DVec2, f64) {
    let curv = |d1: DVec2, d2: DVec2| (d1.normalize(), d1.perp_dot(d2) / d1.length().powi(3));
    match s.entity(k).unwrap().geometry.clone() {
        Geometry::Line { start, end } => {
            let (a, b) = line(s, k);
            let d = if joined_at_start(s, key, start, end) {
                b - a
            } else {
                a - b
            };
            (d.normalize(), 0.0)
        }
        Geometry::Arc { center, start, end } => {
            let c = pos(s, center);
            let r = (pos(s, start) - c).length();
            let q = pos(s, key);
            if joined_at_start(s, key, start, end) {
                ((q - c).perp().normalize(), 1.0 / r)
            } else {
                (-(q - c).perp().normalize(), -1.0 / r)
            }
        }
        Geometry::EllipticalArc { start, end, .. } => {
            let (_, a, b) = ell(s, k);
            let u = unit_frame(s, k, pos(s, key)).normalize();
            // E′(θ) with cos θ = u.x, sin θ = u.y; the arc runs
            // counter-clockwise, which is increasing θ when A × B > 0.
            let mut t = b * u.x - a * u.y;
            if a.perp_dot(b) < 0.0 {
                t = -t;
            }
            let k_ccw = a.perp_dot(b).abs() / t.length().powi(3);
            if joined_at_start(s, key, start, end) {
                (t.normalize(), k_ccw)
            } else {
                (-t.normalize(), -k_ccw)
            }
        }
        Geometry::Conic { start, end, .. } => {
            let ([p0, p1, p2], rho) = cnc(s, k);
            let w = rho / (1.0 - rho);
            // From the joint's end: C = N/D, N = (1−t)²P₀ + 2t(1−t)wP₁ +
            // t²P₂, D = (1−t)² + 2t(1−t)w + t², at t = 0 where D = 1:
            // C′ = N′ − P₀D′, C″ = N″ − 2C′D′ − P₀D″.
            let (q0, q2) = if joined_at_start(s, key, start, end) {
                (p0, p2)
            } else {
                (p2, p0)
            };
            let (n1, n2) = (-2.0 * q0 + 2.0 * w * p1, 2.0 * q0 - 4.0 * w * p1 + 2.0 * q2);
            let (dd1, dd2) = (2.0 * w - 2.0, 4.0 - 4.0 * w);
            let c1 = n1 - q0 * dd1;
            let c2 = n2 - 2.0 * c1 * dd1 - q0 * dd2;
            curv(c1, c2)
        }
        Geometry::ControlSpline { .. } | Geometry::InterpolatedSpline { .. } => {
            let ResolvedCurve::Spline {
                degree,
                knots,
                points,
            } = resolved(s, k)
            else {
                unreachable!()
            };
            let refs = s.entity(k).unwrap().geometry.referenced_points();
            let (mut u, mut p) = (knots.clone(), points.clone());
            if !joined_at_start(s, key, refs[0], refs[refs.len() - 1]) {
                // The reversed curve, t ↦ 1 − t.
                u = knots.iter().rev().map(|x| 1.0 - x).collect();
                p.reverse();
            }
            // Q_i = p(P_{i+1} − P_i)/(u_{i+p+1} − u_{i+1}); C′(0) = Q₀ and
            // C″(0) = (p − 1)(Q₁ − Q₀)/(u_{p+1} − u₂).
            let q = |i: usize| (p[i + 1] - p[i]) * (degree as f64 / (u[i + degree + 1] - u[i + 1]));
            let c1 = q(0);
            let c2 = (q(1) - q(0)) * ((degree as f64 - 1.0) / (u[degree + 1] - u[2]));
            curv(c1, c2)
        }
        _ => panic!("no end curvature"),
    }
}

/// The points joined to `key` by chains of coincident constraints
/// (including `key`).
fn coincident_class(s: &Sketch<u32>, key: u32) -> Vec<u32> {
    let mut class = vec![key];
    loop {
        let before = class.len();
        for c in &s.constraints {
            if let Constraint::Coincident { a, b } = c.constraint {
                if class.contains(&a) && !class.contains(&b) {
                    class.push(b);
                } else if class.contains(&b) && !class.contains(&a) {
                    class.push(a);
                }
            }
        }
        if class.len() == before {
            return class;
        }
    }
}

/// The violation of one constraint, measured from the geometry:
/// `(value, is_angle)`, metres or radians.
pub fn violation(s: &Sketch<u32>, c: &Constraint<u32>) -> (f64, bool) {
    let p = |k: u32| pos(s, k);
    match c {
        Constraint::Coincident { a, b } => ((p(*a) - p(*b)).length(), false),
        Constraint::PointOnCurve { point, curve } => {
            let q = p(*point);
            match s.entity(*curve).unwrap().geometry {
                Geometry::Line { .. } => (dist_to_line(q, line(s, *curve)), false),
                Geometry::Circle { .. } | Geometry::Arc { .. } => {
                    let (cc, r) = circ(s, *curve);
                    (((q - cc).length() - r).abs(), false)
                }
                Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } => {
                    (dist_to_ellipse(s, *curve, q), false)
                }
                Geometry::Conic { .. } => (dist_to_conic(s, *curve, q), false),
                _ => (dist_to_spline(s, *curve, q), false),
            }
        }
        Constraint::Horizontal { line: l } => {
            let (a, b) = line(s, *l);
            ((a.y - b.y).abs(), false)
        }
        Constraint::Vertical { line: l } => {
            let (a, b) = line(s, *l);
            ((a.x - b.x).abs(), false)
        }
        Constraint::HorizontalPoints { a, b } => ((p(*a).y - p(*b).y).abs(), false),
        Constraint::VerticalPoints { a, b } => ((p(*a).x - p(*b).x).abs(), false),
        Constraint::Parallel { a, b } => {
            let (a0, a1) = line(s, *a);
            let (b0, b1) = line(s, *b);
            let (da, db) = ((a1 - a0).normalize(), (b1 - b0).normalize());
            (da.perp_dot(db).abs(), true)
        }
        Constraint::Perpendicular { a, b } => {
            let (a0, a1) = line(s, *a);
            let (b0, b1) = line(s, *b);
            let (da, db) = ((a1 - a0).normalize(), (b1 - b0).normalize());
            (da.dot(db).abs(), true)
        }
        Constraint::Equal { a, b } => match s.entity(*a).unwrap().geometry {
            Geometry::Line { .. } => {
                let (a0, a1) = line(s, *a);
                let (b0, b1) = line(s, *b);
                (((a1 - a0).length() - (b1 - b0).length()).abs(), false)
            }
            _ => ((circ(s, *a).1 - circ(s, *b).1).abs(), false),
        },
        Constraint::Midpoint { point, line: l } => {
            let (a, b) = line(s, *l);
            ((p(*point) - (a + b) * 0.5).length(), false)
        }
        Constraint::Concentric { a, b } => ((circ(s, *a).0 - circ(s, *b).0).length(), false),
        Constraint::Symmetric { a, b, axis } => {
            let (s0, s1) = line(s, *axis);
            let d = (s1 - s0).normalize();
            let q = p(*a) - s0;
            let refl = s0 + d * (2.0 * q.dot(d)) - q;
            ((refl - p(*b)).length(), false)
        }
        Constraint::Curvature { a, b, point } => {
            // As a length: how far the tangents away from the joint are
            // from opposite (a smooth joint; a cusp counts as 2) times the
            // sketch size, and the curvature mismatch for travel through
            // the joint times the size squared, the gap it opens over the
            // sketch.
            let (ta, ka) = away_curvature(s, *a, *point);
            let (tb, kb) = away_curvature(s, *b, *point);
            let l = size(s);
            let turn = (ta + tb).length() * l;
            let gap = (ka + kb).abs() * l * l;
            (turn.max(gap), false)
        }
        Constraint::Tangent { a, b, kind } => match kind {
            TangentKind::AtPoint { point } => {
                let ta = end_tangent(s, *a, *point);
                let tb = end_tangent(s, *b, *point);
                (ta.perp_dot(tb).abs(), true)
            }
            TangentKind::LineCircle { center_left } => {
                let (l, k) = if matches!(s.entity(*a).unwrap().geometry, Geometry::Line { .. }) {
                    (*a, *b)
                } else {
                    (*b, *a)
                };
                let (s0, s1) = line(s, l);
                if let Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } =
                    s.entity(k).unwrap().geometry
                {
                    // In the frame where the ellipse is the unit circle the
                    // line must be at distance 1 from the centre, on the
                    // recorded side (an affine map keeps tangency and sides).
                    // With p = c + [a b]·q, the line n̂·p = k is
                    // ([a b]ᵀn̂)·q = k − n̂·c, so a unit-frame distance D is
                    // D·‖[a b]ᵀn̂‖ in the plane: the gap in metres is the
                    // unit-frame gap times that norm (not times the longer
                    // axis, which overstates it by up to the axis ratio).
                    let (c, a, b) = ell(s, k);
                    let (q0, q1) = (unit_frame(s, k, s0), unit_frame(s, k, s1));
                    let d = q1 - q0;
                    let dist = d.perp_dot(-q0).abs() / d.length();
                    let side = (s1 - s0).perp_dot(c - s0) > 0.0;
                    let gap = if side == *center_left {
                        (dist - 1.0).abs()
                    } else {
                        dist + 1.0
                    };
                    let n = (s1 - s0).perp().normalize();
                    return (gap * DVec2::new(a.dot(n), b.dot(n)).length(), false);
                }
                let (cc, r) = circ(s, k);
                let d = s1 - s0;
                let signed = d.perp_dot(cc - s0) / d.length();
                let want = if *center_left { r } else { -r };
                ((signed - want).abs(), false)
            }
            TangentKind::CircleCircle { internal } => {
                let (ca, ra) = circ(s, *a);
                let (cb, rb) = circ(s, *b);
                let dist = (ca - cb).length();
                let want = if *internal { (ra - rb).abs() } else { ra + rb };
                ((dist - want).abs(), false)
            }
        },
        Constraint::Normal { line: l, curve } => {
            (dist_to_line(circ(s, *curve).0, line(s, *l)), false)
        }
        Constraint::Fix { point, at } => ((p(*point) - *at).length(), false),
        Constraint::Distance { a, b, kind, value } => {
            let ea = &s.entity(*a).unwrap().geometry;
            let eb = &s.entity(*b).unwrap().geometry;
            let m = match (ea, eb) {
                (Geometry::Point { .. }, Geometry::Point { .. }) => {
                    let d = p(*b) - p(*a);
                    match kind {
                        DistanceKind::Aligned => d.length(),
                        DistanceKind::Horizontal => d.x.abs(),
                        DistanceKind::Vertical => d.y.abs(),
                    }
                }
                (Geometry::Point { .. }, Geometry::Line { .. }) => dist_to_line(p(*a), line(s, *b)),
                (Geometry::Line { .. }, Geometry::Point { .. }) => dist_to_line(p(*b), line(s, *a)),
                _ => dist_to_line(line(s, *a).0, line(s, *b)),
            };
            ((m - value).abs(), false)
        }
        Constraint::Length { line: l, value } => {
            let (a, b) = line(s, *l);
            (((b - a).length() - value).abs(), false)
        }
        Constraint::Angle { a, b, value } => {
            let (a0, a1) = line(s, *a);
            let (b0, b1) = line(s, *b);
            let (da, db) = (a1 - a0, b1 - b0);
            let th = da.perp_dot(db).atan2(da.dot(db));
            (wrap(th - value).abs(), true)
        }
        Constraint::Radius { curve, value } => ((circ(s, *curve).1 - value).abs(), false),
        Constraint::Diameter { curve, value } => ((2.0 * circ(s, *curve).1 - value).abs(), false),
        // A plain number: compared like an angle, unscaled.
        Constraint::Rho { curve, value } => ((cnc(s, *curve).1 - value).abs(), true),
        Constraint::Mirror { axis, pairs } => {
            let (s0, s1) = line(s, *axis);
            let d = (s1 - s0).normalize();
            let reflect = |q: DVec2| {
                let v = q - s0;
                s0 + d * (2.0 * v.dot(d)) - v
            };
            let mut worst = 0.0f64;
            for (src, img) in pairs {
                // Image point k against the reflection of source point k,
                // an arc's ends crossed over; circles also by radius.
                let (gs, gi) = (
                    &s.entity(*src).unwrap().geometry,
                    &s.entity(*img).unwrap().geometry,
                );
                let (ps, pi) = (gs.referenced_points(), gi.referenced_points());
                let ps = match gs {
                    Geometry::Arc { center, start, end } => vec![*center, *end, *start],
                    Geometry::EllipticalArc {
                        center,
                        major,
                        minor,
                        start,
                        end,
                    } => vec![*center, *major, *minor, *end, *start],
                    _ => ps,
                };
                assert_eq!(ps.len(), pi.len(), "mirror pair of different kinds");
                for (a, b) in ps.iter().zip(&pi) {
                    worst = worst.max((reflect(p(*a)) - p(*b)).length());
                }
                if let (Geometry::Circle { .. }, Geometry::Circle { .. }) = (gs, gi) {
                    worst = worst.max((circ(s, *src).1 - circ(s, *img).1).abs());
                }
                if let (Geometry::Conic { .. }, Geometry::Conic { .. }) = (gs, gi) {
                    worst = worst.max((cnc(s, *src).1 - cnc(s, *img).1).abs() * size(s));
                }
            }
            (worst, false)
        }
        Constraint::Offset {
            distance, pieces, ..
        } => (offset_violation(s, *distance, pieces), false),
        Constraint::LinearPattern { direction, copies } => {
            let (s0, s1) = line(s, *direction);
            (
                pattern_violation(s, copies, |q, k| q + (s1 - s0) * f64::from(k)),
                false,
            )
        }
        Constraint::CircularPattern {
            centre,
            angle,
            copies,
        } => {
            let c = p(*centre);
            (
                pattern_violation(s, copies, |q, k| {
                    c + DVec2::from_angle(f64::from(k) * angle).rotate(q - c)
                }),
                false,
            )
        }
    }
}

/// How far a pattern's copies are from their sources moved by `map` (a
/// point and the copy number): every image point against its source
/// point's, circles also by radius.
pub fn pattern_violation(
    s: &Sketch<u32>,
    copies: &[(u32, u32, u32)],
    map: impl Fn(DVec2, u32) -> DVec2,
) -> f64 {
    let mut worst = 0.0f64;
    for &(src, img, k) in copies {
        let (gs, gi) = (
            &s.entity(src).unwrap().geometry,
            &s.entity(img).unwrap().geometry,
        );
        let (ps, pi) = (gs.referenced_points(), gi.referenced_points());
        assert_eq!(ps.len(), pi.len(), "pattern copy of a different kind");
        for (a, b) in ps.iter().zip(&pi) {
            worst = worst.max((map(pos(s, *a), k) - pos(s, *b)).length());
        }
        if let (Geometry::Circle { .. }, Geometry::Circle { .. }) = (gs, gi) {
            worst = worst.max((circ(s, src).1 - circ(s, img).1).abs());
        }
        if let (Geometry::Conic { .. }, Geometry::Conic { .. }) = (gs, gi) {
            worst = worst.max((cnc(s, src).1 - cnc(s, img).1).abs() * size(s));
        }
    }
    worst
}

/// Bounding-box diagonal of the sketch's points and circles.
pub fn size(s: &Sketch<u32>) -> f64 {
    let mut lo = DVec2::splat(f64::INFINITY);
    let mut hi = DVec2::splat(f64::NEG_INFINITY);
    for e in &s.entities {
        match e.geometry {
            Geometry::Point { at } => {
                lo = lo.min(at);
                hi = hi.max(at);
            }
            Geometry::Circle { center, radius } => {
                let c = pos(s, center);
                lo = lo.min(c - DVec2::splat(radius));
                hi = hi.max(c + DVec2::splat(radius));
            }
            _ => {}
        }
    }
    if lo.x > hi.x {
        0.0
    } else {
        (hi - lo).length()
    }
}

/// The largest violation, scaled: lengths relative to `scale`, angles in
/// radians. Arcs' implicit equal-radius condition is included.
pub fn worst(s: &Sketch<u32>, scale: f64) -> (f64, String) {
    let mut w = (0.0f64, String::new());
    for c in &s.constraints {
        let (v, ang) = violation(s, &c.constraint);
        let rel = if ang { v } else { v / scale };
        // A NaN violation counts as the worst.
        if rel.is_nan() || rel > w.0 {
            w = (rel, format!("{:?} {:?}", c.id, c.constraint));
        }
    }
    for e in &s.entities {
        if let Geometry::Arc { center, start, end } = e.geometry {
            let c = pos(s, center);
            let v = ((pos(s, end) - c).length() - (pos(s, start) - c).length()).abs() / scale;
            if v.is_nan() || v > w.0 {
                w = (v, format!("arc {:?} radius mismatch", e.id));
            }
        }
        if let Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. } = e.geometry {
            let (_, a, b) = ell(s, e.id);
            let v = a.normalize().dot(b.normalize()).abs();
            if v.is_nan() || v > w.0 {
                w = (v, format!("ellipse {:?} axes not perpendicular", e.id));
            }
        }
        if let Geometry::EllipticalArc { start, end, .. } = e.geometry {
            for q in [start, end] {
                let v = dist_to_ellipse(s, e.id, pos(s, q)) / scale;
                if v.is_nan() || v > w.0 {
                    w = (v, format!("elliptical arc {:?} end off its ellipse", e.id));
                }
            }
        }
    }
    w
}

/// Asserts every constraint holds: lengths within `1e-9` of the sketch
/// size, angles within `1e-9` rad.
pub fn assert_satisfied(s: &Sketch<u32>) {
    let scale = size(s).max(1e-300);
    let (w, what) = worst(s, scale);
    assert!(w <= 1e-9, "violation {w:e} (relative) at {what}");
}

/// Bits of every coordinate and radius, for determinism checks.
pub fn bits(s: &Sketch<u32>) -> Vec<u64> {
    let mut out = Vec::new();
    for e in &s.entities {
        match e.geometry {
            Geometry::Point { at } => {
                out.push(at.x.to_bits());
                out.push(at.y.to_bits());
            }
            Geometry::Circle { radius, .. } => out.push(radius.to_bits()),
            Geometry::Conic { rho, .. } => out.push(rho.to_bits()),
            _ => {}
        }
    }
    out
}

/// How far an offset's images are from where design note 0007 §10 puts
/// them, measured from the geometry: every image end on its source's
/// offset curve (a circle's radius by the distance), and an image end that
/// no other image shares on the normal at its source's end.
pub fn offset_violation(s: &Sketch<u32>, d: f64, pieces: &[(u32, u32, bool)]) -> f64 {
    let ends = |k: u32| match s.entity(k).unwrap().geometry {
        Geometry::Line { start, end } | Geometry::Arc { start, end, .. } => Some((start, end)),
        _ => None,
    };
    let mut uses: std::collections::BTreeMap<u32, usize> = Default::default();
    for (_, img, _) in pieces {
        if let Some((a, b)) = ends(*img) {
            *uses.entry(a).or_default() += 1;
            *uses.entry(b).or_default() += 1;
        }
    }
    let mut worst = 0.0f64;
    for &(src, img, left) in pieces {
        let side = if left { 1.0 } else { -1.0 };
        match (&s.entity(src).unwrap().geometry, ends(img)) {
            (Geometry::Line { .. }, Some((ia, ib))) => {
                let (a, b) = line(s, src);
                let u = (b - a).normalize();
                for (x, v) in [(ia, a), (ib, b)] {
                    let q = pos(s, x);
                    worst = worst.max((u.perp_dot(q - a) - side * d).abs());
                    if uses[&x] == 1 {
                        worst = worst.max((q - v).dot(u).abs());
                    }
                }
            }
            (Geometry::Arc { start, end, .. }, Some((ia, ib))) => {
                let (c, r) = circ(s, src);
                for (x, v) in [(ia, *start), (ib, *end)] {
                    let q = pos(s, x);
                    worst = worst.max(((q - c).length() - (r - side * d)).abs());
                    if uses[&x] == 1 {
                        let n = (pos(s, v) - c).normalize();
                        worst = worst.max(n.perp_dot(q - c).abs());
                    }
                }
            }
            (Geometry::Circle { .. }, None) => {
                let ((c, r), (ic, ir)) = (circ(s, src), circ(s, img));
                worst = worst
                    .max((c - ic).length())
                    .max((ir - (r - side * d)).abs());
            }
            _ => panic!("offset piece of different kinds"),
        }
    }
    worst
}
