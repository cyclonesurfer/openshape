//! Random sketches whose geometry satisfies their constraints by
//! construction (design note 0005 §8, "constructive property tests"):
//! motifs place entities in relation to each other (coincident ends,
//! parallel lines, tangent circles, points on curves…), constraints are
//! added that hold, and dimensions are measured from the geometry.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use openshape_sketch_solver::{
    measure, resolve_curves, Constraint, DVec2, DistanceKind, Geometry, ResolvedCurve, Sketch,
    TangentKind,
};

use super::{de_boor, B};

/// SplitMix64: a small deterministic generator for test geometry.
#[derive(Clone)]
pub struct Rng(pub u64);

impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    pub fn f(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.f()
    }
    pub fn int(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    pub fn chance(&mut self, p: f64) -> bool {
        self.f() < p
    }
    pub fn pt(&mut self) -> DVec2 {
        DVec2::new(self.range(-1.0, 1.0), self.range(-1.0, 1.0))
    }
    pub fn dir(&mut self) -> DVec2 {
        let a = self.range(-PI, PI);
        DVec2::new(a.cos(), a.sin())
    }
}

/// What the generator built, for tests that need handles.
pub struct Generated {
    pub sketch: Sketch<u32>,
    /// Point entities.
    pub points: Vec<u32>,
    /// Line entities.
    pub lines: Vec<u32>,
    /// Circle and arc entities.
    pub circulars: Vec<u32>,
}

struct G {
    b: B,
    rng: Rng,
    points: Vec<u32>,
    lines: Vec<u32>,
    circulars: Vec<u32>,
    circles: Vec<u32>,
}

impl G {
    fn point(&mut self, p: DVec2) -> u32 {
        let k = self.b.point(p.x, p.y);
        self.points.push(k);
        k
    }
    fn line(&mut self, a: DVec2, c: DVec2) -> (u32, u32, u32) {
        let s = self.point(a);
        let e = self.point(c);
        let l = self.b.line_between(s, e);
        self.lines.push(l);
        (l, s, e)
    }
    fn circle(&mut self, c: DVec2, r: f64) -> (u32, u32) {
        let cc = self.point(c);
        let k = self.b.entity(Geometry::Circle {
            center: cc,
            radius: r,
        });
        self.circulars.push(k);
        self.circles.push(k);
        (k, cc)
    }
    fn arc(&mut self, c: DVec2, r: f64, a0: f64, sweep: f64) -> (u32, u32, u32, u32) {
        let cc = self.point(c);
        let s = self.point(c + r * DVec2::new(a0.cos(), a0.sin()));
        let a1 = a0 + sweep;
        let e = self.point(c + r * DVec2::new(a1.cos(), a1.sin()));
        let k = self.b.entity(Geometry::Arc {
            center: cc,
            start: s,
            end: e,
        });
        self.circulars.push(k);
        (k, cc, s, e)
    }
    fn con(&mut self, c: Constraint<u32>) {
        self.b.c(c);
    }
    fn len(&mut self) -> f64 {
        self.rng.range(0.25, 1.0)
    }

    fn polygon(&mut self) {
        let k = 3 + self.rng.int(3);
        let center = self.rng.pt() * 0.5;
        let rad = self.rng.range(0.4, 0.9);
        let base = self.rng.range(-PI, PI);
        let verts: Vec<DVec2> = (0..k)
            .map(|i| {
                let a = base + TAU * (i as f64 + self.rng.range(-0.25, 0.25)) / k as f64;
                center + rad * self.rng.range(0.7, 1.0) * DVec2::new(a.cos(), a.sin())
            })
            .collect();
        let closed = self.rng.chance(0.6);
        let n_edges = if closed { k } else { k - 1 };
        let mut ends = Vec::new();
        for i in 0..n_edges {
            let (_, s, e) = self.line(verts[i], verts[(i + 1) % k]);
            ends.push((s, e));
        }
        for i in 0..n_edges {
            if i + 1 < n_edges || closed {
                let next = ends[(i + 1) % n_edges].0;
                self.con(Constraint::Coincident { a: ends[i].1, b: next });
            }
        }
    }

    fn rectangle(&mut self) {
        let c = self.rng.pt() * 0.5;
        let (w, h) = (self.rng.range(0.3, 1.0), self.rng.range(0.3, 1.0));
        let p = [
            c + DVec2::new(-w, -h) * 0.5,
            c + DVec2::new(w, -h) * 0.5,
            c + DVec2::new(w, h) * 0.5,
            c + DVec2::new(-w, h) * 0.5,
        ];
        let mut ends = Vec::new();
        let mut ls = Vec::new();
        for i in 0..4 {
            let (l, s, e) = self.line(p[i], p[(i + 1) % 4]);
            ends.push((s, e));
            ls.push(l);
        }
        for i in 0..4 {
            self.con(Constraint::Coincident {
                a: ends[i].1,
                b: ends[(i + 1) % 4].0,
            });
            if i % 2 == 0 {
                self.con(Constraint::Horizontal { line: ls[i] });
            } else {
                self.con(Constraint::Vertical { line: ls[i] });
            }
        }
        if self.rng.chance(0.5) {
            self.con(Constraint::Length {
                line: ls[0],
                value: 1.0,
            });
        }
        if self.rng.chance(0.5) {
            // Diagonal corners: horizontal and vertical distance dimensions.
            self.con(Constraint::Distance {
                a: ends[0].0,
                b: ends[2].0,
                kind: DistanceKind::Vertical,
                value: 1.0,
            });
        }
    }

    fn line_pair(&mut self) {
        let a0 = self.rng.pt();
        let d = self.rng.dir();
        let la = self.len();
        let (l1, _, _) = self.line(a0, a0 + d * la);
        let b0 = self.rng.pt();
        let kind = self.rng.int(3);
        let db = match kind {
            0 => d,
            1 => -d,
            _ => d.perp(),
        };
        let equal = self.rng.chance(0.4);
        let lb = if equal { la } else { self.len() };
        // Keep parallel lines apart so the line-line distance is not tiny.
        let b0 = if kind < 2 && d.perp_dot(b0 - a0).abs() < 0.1 {
            b0 + d.perp() * 0.3
        } else {
            b0
        };
        let (l2, _, _) = self.line(b0, b0 + db * lb);
        if kind < 2 {
            self.con(Constraint::Parallel { a: l1, b: l2 });
            if self.rng.chance(0.5) {
                self.con(Constraint::Distance {
                    a: l1,
                    b: l2,
                    kind: DistanceKind::Aligned,
                    value: 1.0,
                });
            }
        } else {
            self.con(Constraint::Perpendicular { a: l1, b: l2 });
        }
        if equal {
            self.con(Constraint::Equal { a: l1, b: l2 });
        }
    }

    fn point_on_line(&mut self) {
        let a = self.rng.pt();
        let d = self.rng.dir();
        let l = self.len();
        let (ln, s, e) = self.line(a, a + d * l);
        let t = self.rng.range(-0.5, 1.5);
        let p = self.point(a + d * l * t);
        self.con(Constraint::PointOnCurve { point: p, curve: ln });
        if self.rng.chance(0.4) {
            let m = self.point(a + d * l * 0.5);
            self.con(Constraint::Midpoint { point: m, line: ln });
        }
        if self.rng.chance(0.3) {
            // A point off the line with a point–line distance dimension.
            let (u, h) = (self.rng.range(0.0, 1.0), self.rng.range(0.1, 0.5));
            let q = self.point(a + d * l * u + d.perp() * h);
            self.con(Constraint::Distance {
                a: q,
                b: ln,
                kind: DistanceKind::Aligned,
                value: 1.0,
            });
        }
        let _ = (s, e);
    }

    fn symmetric(&mut self) {
        let a = self.rng.pt();
        let d = self.rng.dir();
        let l = self.len();
        let (ax, _, _) = self.line(a, a + d * l);
        let n = d.perp();
        let along = self.rng.range(-0.5, 1.5) * l;
        let off = self.rng.range(0.1, 0.6) * if self.rng.chance(0.5) { 1.0 } else { -1.0 };
        let p = a + d * along + n * off;
        let q = a + d * along - n * off;
        let (pp, qq) = (self.point(p), self.point(q));
        self.con(Constraint::Symmetric {
            a: pp,
            b: qq,
            axis: ax,
        });
    }

    fn circles(&mut self) {
        let c = self.rng.pt() * 0.6;
        let r = self.rng.range(0.15, 0.5);
        let (k, cc) = self.circle(c, r);
        if self.rng.chance(0.5) {
            let a = self.rng.range(-PI, PI);
            let p = self.point(c + r * DVec2::new(a.cos(), a.sin()));
            self.con(Constraint::PointOnCurve { point: p, curve: k });
        }
        if self.rng.chance(0.4) {
            let r2 = r + self.rng.range(0.1, 0.4);
            let (k2, _) = self.circle(c, r2);
            self.con(Constraint::Concentric { a: k, b: k2 });
        }
        if self.rng.chance(0.4) {
            let c3 = self.rng.pt();
            let (k3, _) = self.circle(c3, r);
            self.con(Constraint::Equal { a: k, b: k3 });
        }
        if self.rng.chance(0.5) {
            let internal = self.rng.chance(0.4);
            let r4 = if internal {
                r * self.rng.range(0.2, 0.7)
            } else {
                self.rng.range(0.1, 0.4)
            };
            let dist = if internal { r - r4 } else { r + r4 };
            let c4 = c + self.rng.dir() * dist;
            let (k4, _) = self.circle(c4, r4);
            self.con(Constraint::Tangent {
                a: k,
                b: k4,
                kind: TangentKind::CircleCircle { internal },
            });
        }
        if self.rng.chance(0.5) {
            let d = self.rng.dir();
            let left = self.rng.chance(0.5);
            let sigma = if left { 1.0 } else { -1.0 };
            let q = c - d.perp() * (sigma * r);
            let t0 = self.rng.range(-0.8, 0.0);
            let ll = self.len();
            let (l, _, _) = self.line(q + d * t0, q + d * (t0 + ll));
            let (a, b) = if self.rng.chance(0.5) { (l, k) } else { (k, l) };
            self.con(Constraint::Tangent {
                a,
                b,
                kind: TangentKind::LineCircle { center_left: left },
            });
        }
        if self.rng.chance(0.3) {
            let d = self.rng.dir();
            let t0 = self.rng.range(-0.6, 0.3);
            let ll = self.len();
            let (l, _, _) = self.line(c + d * t0, c + d * (t0 + ll));
            self.con(Constraint::Normal { line: l, curve: k });
        }
        if self.rng.chance(0.4) {
            self.con(Constraint::Radius {
                curve: k,
                value: 1.0,
            });
        } else if self.rng.chance(0.3) {
            self.con(Constraint::Diameter {
                curve: k,
                value: 1.0,
            });
        }
        let _ = cc;
    }

    fn arcs(&mut self) {
        let c = self.rng.pt() * 0.5;
        let r = self.rng.range(0.2, 0.6);
        let a0 = self.rng.range(-PI, PI);
        let sweep = self.rng.range(0.6, 4.5);
        let (k, cc, s, e) = self.arc(c, r, a0, sweep);
        let a1 = a0 + sweep;
        let pe = c + r * DVec2::new(a1.cos(), a1.sin());
        let tangent = DVec2::new(-a1.sin(), a1.cos());
        match self.rng.int(3) {
            0 => {
                // A line leaving the arc's end along its tangent.
                let ll = self.len();
                let (l, ls, _) = self.line(pe, pe + tangent * ll);
                self.con(Constraint::Coincident { a: e, b: ls });
                self.con(Constraint::Tangent {
                    a: k,
                    b: l,
                    kind: TangentKind::AtPoint { point: ls },
                });
            }
            1 => {
                // A second arc continuing smoothly (centre on the radius line).
                let r2 = self.rng.range(0.2, 0.6);
                let toward = (c - pe).normalize();
                let same_side = self.rng.chance(0.5);
                let c2 = if same_side { pe + toward * r2 } else { pe - toward * r2 };
                let b0 = (pe - c2).y.atan2((pe - c2).x);
                let sweep2 = self.rng.range(0.6, 3.0);
                let (k2, _, s2, _) = self.arc(c2, r2, b0, sweep2);
                self.con(Constraint::Coincident { a: e, b: s2 });
                self.con(Constraint::Tangent {
                    a: k,
                    b: k2,
                    kind: TangentKind::AtPoint { point: e },
                });
            }
            _ => {
                let a = self.rng.range(-PI, PI);
                let p = self.point(c + r * DVec2::new(a.cos(), a.sin()));
                self.con(Constraint::PointOnCurve { point: p, curve: k });
            }
        }
        if self.rng.chance(0.4) {
            self.con(Constraint::Radius {
                curve: k,
                value: 1.0,
            });
        }
        let _ = (cc, s);
    }

    fn splines(&mut self) {
        let n = 3 + self.rng.int(4);
        let x0 = self.rng.range(-1.0, 0.0);
        let y0 = self.rng.range(-1.0, 1.0);
        let pts: Vec<DVec2> = (0..n)
            .map(|i| DVec2::new(x0 + 0.3 * i as f64, y0 + self.rng.range(-0.3, 0.3)))
            .collect();
        let keys: Vec<u32> = pts.iter().map(|&p| self.point(p)).collect();
        if self.rng.chance(0.5) {
            let degree = 1 + self.rng.int((n - 1).min(5));
            let sp = self.b.entity(Geometry::ControlSpline {
                control: keys.clone(),
                degree: degree as u32,
            });
            let knots = control_knots(n, degree);
            let t = self.rng.range(0.05, 0.95);
            let q = self.point(de_boor(degree, &knots, &pts, t));
            self.con(Constraint::PointOnCurve { point: q, curve: sp });
            if self.rng.chance(0.5) {
                self.spline_tangent_line(sp, &keys, true);
            }
        } else {
            // Chord-length parameters.
            let mut params = vec![0.0];
            for w in pts.windows(2) {
                let last = *params.last().unwrap();
                params.push(last + (w[1] - w[0]).length());
            }
            let total = *params.last().unwrap();
            for p in params.iter_mut() {
                *p /= total;
            }
            *params.last_mut().unwrap() = 1.0;
            let sp = self.b.entity(Geometry::InterpolatedSpline {
                through: keys.clone(),
                params,
            });
            if self.rng.chance(0.6) {
                self.spline_tangent_line(sp, &keys, false);
            }
            if self.rng.chance(0.5) {
                // A point on the curve, from its resolved form.
                if let Some((_, _, ResolvedCurve::Spline { degree, knots, points })) =
                    resolve_curves(&self.b.sketch)
                        .ok()
                        .and_then(|v| v.into_iter().find(|(k, _, _)| *k == sp))
                {
                    let t = self.rng.range(0.05, 0.95);
                    let q = self.point(de_boor(degree, &knots, &points, t));
                    self.con(Constraint::PointOnCurve { point: q, curve: sp });
                }
            }
        }
    }

    /// A line leaving the spline's last point along its end tangent.
    fn spline_tangent_line(&mut self, sp: u32, keys: &[u32], _control: bool) {
        let Some((_, _, ResolvedCurve::Spline { points, .. })) = resolve_curves(&self.b.sketch)
            .ok()
            .and_then(|v| v.into_iter().find(|(k, _, _)| *k == sp))
        else {
            return;
        };
        let m = points.len();
        let leg = points[m - 1] - points[m - 2];
        if leg.length() < 1e-3 {
            return;
        }
        let end = points[m - 1];
        let ll = self.len();
        let (l, ls, _) = self.line(end, end + leg.normalize() * ll);
        self.con(Constraint::Coincident {
            a: *keys.last().unwrap(),
            b: ls,
        });
        self.con(Constraint::Tangent {
            a: sp,
            b: l,
            kind: TangentKind::AtPoint { point: ls },
        });
    }

    fn aligned_points(&mut self) {
        if self.points.is_empty() {
            return;
        }
        let k = self.points[self.rng.int(self.points.len())];
        let p = self.b.sketch.point(k).unwrap();
        let horizontal = self.rng.chance(0.5);
        let q = if horizontal {
            DVec2::new(p.x + self.rng.range(0.2, 0.8), p.y)
        } else {
            DVec2::new(p.x, p.y + self.rng.range(0.2, 0.8))
        };
        let kq = self.point(q);
        self.con(if horizontal {
            Constraint::HorizontalPoints { a: k, b: kq }
        } else {
            Constraint::VerticalPoints { a: k, b: kq }
        });
    }

    fn dimensions(&mut self) {
        let n = self.rng.int(4);
        for _ in 0..n {
            if self.points.len() < 2 {
                return;
            }
            let a = self.points[self.rng.int(self.points.len())];
            let b = self.points[self.rng.int(self.points.len())];
            if a == b {
                continue;
            }
            let d = self.b.sketch.point(b).unwrap() - self.b.sketch.point(a).unwrap();
            let kind = match self.rng.int(3) {
                0 if d.length() > 0.1 => DistanceKind::Aligned,
                1 if d.x.abs() > 0.1 => DistanceKind::Horizontal,
                2 if d.y.abs() > 0.1 => DistanceKind::Vertical,
                _ => continue,
            };
            self.con(Constraint::Distance {
                a,
                b,
                kind,
                value: 1.0,
            });
        }
        if self.lines.len() >= 2 && self.rng.chance(0.5) {
            let a = self.lines[self.rng.int(self.lines.len())];
            let b = self.lines[self.rng.int(self.lines.len())];
            if a != b {
                self.con(Constraint::Angle { a, b, value: 0.0 });
            }
        }
        if !self.lines.is_empty() && self.rng.chance(0.5) {
            let l = self.lines[self.rng.int(self.lines.len())];
            self.con(Constraint::Length { line: l, value: 1.0 });
        }
    }
}

/// Clamped uniform knots, computed independently of the solver.
pub fn control_knots(n: usize, p: usize) -> Vec<f64> {
    let segs = n - p;
    let mut u = vec![0.0; p + 1];
    for j in 1..segs {
        u.push(j as f64 / segs as f64);
    }
    u.extend(std::iter::repeat_n(1.0, p + 1));
    u
}

/// Generates a sketch from `seed` with about `motifs` motifs, scaled by
/// `scale` and shifted by `offset · scale`. Dimension values and fixed
/// positions are measured from the final geometry.
pub fn generate(seed: u64, motifs: usize, scale: f64, offset: DVec2) -> Generated {
    let mut g = G {
        b: B::new(),
        rng: Rng(seed),
        points: Vec::new(),
        lines: Vec::new(),
        circulars: Vec::new(),
        circles: Vec::new(),
    };
    for _ in 0..motifs {
        match g.rng.int(9) {
            0 => g.polygon(),
            1 => g.rectangle(),
            2 => g.line_pair(),
            3 => g.point_on_line(),
            4 => g.symmetric(),
            5 => g.circles(),
            6 => g.arcs(),
            7 => g.splines(),
            _ => g.aligned_points(),
        }
    }
    g.dimensions();
    // Fix some points: none, a few, many or most.
    let fix_chance = [0.0, 0.12, 0.35, 0.7][g.rng.int(4)];
    let pts = g.points.clone();
    for &p in &pts {
        if g.rng.chance(fix_chance) {
            g.con(Constraint::Fix {
                point: p,
                at: DVec2::ZERO,
            });
        }
    }
    // Transform to the requested scale and place.
    let mut s = g.b.sketch.clone();
    for e in s.entities.iter_mut() {
        match &mut e.geometry {
            Geometry::Point { at } => *at = *at * scale + offset * scale,
            Geometry::Circle { radius, .. } => *radius *= scale,
            _ => {}
        }
    }
    // Fixed positions and dimension values from the geometry.
    let snapshot = s.clone();
    for c in s.constraints.iter_mut() {
        match &mut c.constraint {
            Constraint::Fix { point, at } => *at = snapshot.point(*point).unwrap(),
            other if other.value().is_some() => {
                let v = measure(&snapshot, other).expect("measurable dimension");
                other.set_value(v);
            }
            _ => {}
        }
    }
    // Angles at exactly ±π are stored as π (the (−π, π] convention).
    for c in s.constraints.iter_mut() {
        if let Constraint::Angle { value, .. } = &mut c.constraint {
            if *value <= -PI {
                *value = PI;
            }
            let _ = FRAC_PI_2;
        }
    }
    Generated {
        sketch: s,
        points: g.points,
        lines: g.lines,
        circulars: g.circulars,
    }
}

/// Moves every point by up to `amount · scale` in each coordinate and
/// scales circle radii by up to `1 ± amount`.
pub fn perturb(s: &Sketch<u32>, seed: u64, amount: f64, scale: f64) -> Sketch<u32> {
    let mut rng = Rng(seed ^ 0x5DEE_CE66_D1CE_4E5B);
    let mut out = s.clone();
    for e in out.entities.iter_mut() {
        match &mut e.geometry {
            Geometry::Point { at } => {
                *at += DVec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)) * (amount * scale);
            }
            Geometry::Circle { radius, .. } => *radius *= 1.0 + amount * rng.range(-1.0, 1.0),
            _ => {}
        }
    }
    out
}
