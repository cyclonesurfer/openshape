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
    /// A line length: long enough that perturbing its ends (by up to 0.1
    /// in each coordinate, see [`perturb`]) cannot swap them.
    fn len(&mut self) -> f64 {
        self.rng.range(0.4, 1.0)
    }

    fn polygon(&mut self) {
        let k = 3 + self.rng.int(3);
        let center = self.rng.pt() * 0.5;
        let rad = self.rng.range(0.6, 1.0);
        let base = self.rng.range(-PI, PI);
        // Edges of at least ~0.35.
        let verts: Vec<DVec2> = (0..k)
            .map(|i| {
                let a = base + TAU * (i as f64 + self.rng.range(-0.15, 0.15)) / k as f64;
                center + rad * self.rng.range(0.85, 1.0) * DVec2::new(a.cos(), a.sin())
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
                self.con(Constraint::Coincident {
                    a: ends[i].1,
                    b: next,
                });
            }
        }
    }

    fn rectangle(&mut self) {
        let c = self.rng.pt() * 0.5;
        let (w, h) = (self.rng.range(0.45, 1.0), self.rng.range(0.45, 1.0));
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
        // Keep parallel lines well apart: the line–line distance is
        // unsigned, and the perturbation must not move one line across the
        // other (a different branch).
        let off = d.perp_dot(b0 - a0);
        let b0 = if kind < 2 && off.abs() < 0.35 {
            b0 + d.perp() * (0.5 * if off < 0.0 { -1.0 } else { 1.0 })
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
        let (ln, _, _) = self.line(a, a + d * l);
        // Away from the midpoint (which may be constrained too): two
        // constrained points close together would pin the line's direction
        // only weakly.
        let side = if self.rng.chance(0.5) { 1.0 } else { -1.0 };
        let t = 0.5 + side * self.rng.range(0.3, 1.0);
        let p = self.point(a + d * l * t);
        self.con(Constraint::PointOnCurve {
            point: p,
            curve: ln,
        });
        if self.rng.chance(0.4) {
            let m = self.point(a + d * l * 0.5);
            self.con(Constraint::Midpoint { point: m, line: ln });
        }
        if self.rng.chance(0.3) {
            // A point off the line with a point–line distance dimension.
            let (u, h) = (self.rng.range(0.0, 1.0), self.rng.range(0.35, 0.6));
            let q = self.point(a + d * l * u + d.perp() * h);
            self.con(Constraint::Distance {
                a: q,
                b: ln,
                kind: DistanceKind::Aligned,
                value: 1.0,
            });
        }
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
        let (k, _) = self.circle(c, r);
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
    }

    fn arcs(&mut self) {
        let c = self.rng.pt() * 0.5;
        let r = self.rng.range(0.2, 0.6);
        let a0 = self.rng.range(-PI, PI);
        let sweep = self.rng.range(0.6, 4.5);
        let (k, _, _, e) = self.arc(c, r, a0, sweep);
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
                let c2 = if same_side {
                    pe + toward * r2
                } else {
                    pe - toward * r2
                };
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
            // Inside a knot span, away from its ends: a degree-1 spline has
            // corners there, and a point near a corner may be closer to the
            // other segment after the perturbation (a different branch).
            let segs = n - degree;
            let j = self.rng.int(segs);
            let t = (j as f64 + self.rng.range(0.2, 0.8)) / segs as f64;
            let at = de_boor(degree, &knots, &pts, t);
            // On a polyline, a point close to the line of another segment
            // is ambiguous at the perturbation's scale (it could be on
            // either): keep it at least 0.3 from those lines.
            let clear = degree > 1
                || (0..n - 1).filter(|&i| i != j).all(|i| {
                    let d = pts[i + 1] - pts[i];
                    d.length() > 0.0 && d.perp_dot(at - pts[i]).abs() / d.length() >= 0.3
                });
            if clear {
                let q = self.point(at);
                self.con(Constraint::PointOnCurve {
                    point: q,
                    curve: sp,
                });
            }
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
                if let Some((
                    _,
                    _,
                    ResolvedCurve::Spline {
                        degree,
                        knots,
                        points,
                    },
                )) = resolve_curves(&self.b.sketch)
                    .ok()
                    .and_then(|v| v.into_iter().find(|(k, _, _)| *k == sp))
                {
                    let t = self.rng.range(0.05, 0.95);
                    let q = self.point(de_boor(degree, &knots, &points, t));
                    self.con(Constraint::PointOnCurve {
                        point: q,
                        curve: sp,
                    });
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

    /// An ellipse or elliptical arc of either handedness, with a point on
    /// it, a line tangent to it, or a line leaving an arc's end along its
    /// tangent; sometimes construction lines along its semi-axes (which
    /// the dimensions motif may then dimension).
    /// An eccentric anomaly in [0, 2π) at least 0.3 from a multiple of π/2.
    fn off_vertex(&mut self) -> f64 {
        self.rng.int(4) as f64 * FRAC_PI_2 + self.rng.range(0.3, FRAC_PI_2 - 0.3)
    }

    fn ellipses(&mut self) {
        let c = self.rng.pt() * 0.5;
        let u = self.rng.dir();
        let la = self.rng.range(0.3, 0.7);
        let lb = la * self.rng.range(0.3, 0.9);
        let a = u * la;
        let b = if self.rng.chance(0.3) {
            -u.perp() * lb
        } else {
            u.perp() * lb
        };
        let at = |t: f64| c + a * t.cos() + b * t.sin();
        let tangent = |t: f64| b * t.cos() - a * t.sin();
        let (kc, km, kn) = (self.point(c), self.point(c + a), self.point(c + b));
        // Counter-clockwise in the plane is increasing θ when b is left
        // of a.
        let ccw = a.perp_dot(b) > 0.0;
        let arc = self.rng.chance(0.6);
        let (k, end) = if arc {
            // Ends at least 0.3 rad from a vertex, where the tangent's
            // direction does not depend on the other axis's length (a fold
            // when an angle or a tangency holds that direction).
            let t0 = self.off_vertex();
            let t1 = loop {
                let t1 = self.off_vertex();
                let sweep = if ccw { t1 - t0 } else { t0 - t1 }.rem_euclid(TAU);
                if (0.6..=4.5).contains(&sweep) {
                    break if ccw { t0 + sweep } else { t0 - sweep };
                }
            };
            let (ks, ke) = (self.point(at(t0)), self.point(at(t1)));
            let k = self.b.entity(Geometry::EllipticalArc {
                center: kc,
                major: km,
                minor: kn,
                start: ks,
                end: ke,
            });
            (k, Some((ke, t1)))
        } else {
            let k = self.b.entity(Geometry::Ellipse {
                center: kc,
                major: km,
                minor: kn,
            });
            (k, None)
        };
        match (self.rng.int(3), end) {
            (0, Some((ke, t1))) => {
                let d = tangent(t1).normalize();
                let ll = self.len();
                let (l, ls, _) = self.line(at(t1), at(t1) + d * ll);
                self.con(Constraint::Coincident { a: ke, b: ls });
                self.con(Constraint::Tangent {
                    a: k,
                    b: l,
                    kind: TangentKind::AtPoint { point: ls },
                });
            }
            (1, _) => {
                // A line touching the full ellipse at θ, centred on it, at
                // least 0.3 rad from a vertex: touching at a vertex, the
                // tangency does not depend on the other axis's length (its
                // derivative vanishes), a fold where that axis is otherwise
                // free.
                let t = self.off_vertex();
                let d = tangent(t).normalize();
                let (h0, h1) = (self.rng.range(0.3, 0.8), self.rng.range(0.3, 0.8));
                let (l, _, _) = self.line(at(t) - d * h0, at(t) + d * h1);
                self.con(Constraint::Tangent {
                    a: l,
                    b: k,
                    kind: TangentKind::LineCircle {
                        center_left: d.perp_dot(c - at(t)) > 0.0,
                    },
                });
            }
            _ => {
                let t = self.rng.range(-PI, PI);
                let q = self.point(at(t));
                self.con(Constraint::PointOnCurve { point: q, curve: k });
            }
        }
        if self.rng.chance(0.5) {
            let l = self.b.line_between(kc, km);
            self.lines.push(l);
            let l2 = self.b.line_between(kc, kn);
            self.lines.push(l2);
        }
    }

    /// A conic arc (design note 0009) of any type, sometimes with its ρ
    /// dimensioned, with a point on it, a line leaving an end along its
    /// tangent, or an arc continuing tangentially from its start.
    fn conics(&mut self) {
        let p0 = self.rng.pt() * 0.5;
        let chord = self.rng.dir() * self.rng.range(0.5, 1.0);
        let p2 = p0 + chord;
        // The control point above the chord, anywhere along it (the legs
        // well apart from the chord, so that a perturbation cannot flatten
        // the triangle).
        let along = self.rng.range(0.1, 0.9);
        let height = self.rng.range(0.4, 0.8) * if self.rng.chance(0.5) { 1.0 } else { -1.0 };
        let p1 = p0 + chord * along + chord.perp().normalize() * height;
        let rho = self.rng.range(0.15, 0.85);
        let w = rho / (1.0 - rho);
        // The curve by its homogeneous de Casteljau form, written here
        // independently of the solver.
        let at = |t: f64| {
            let (a, b) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t) * w);
            (p0 * a + p1 * b + p2 * (t * t)) / (a + b + t * t)
        };
        let (k0, k1, k2) = (self.point(p0), self.point(p1), self.point(p2));
        let k = self.b.entity(Geometry::Conic {
            start: k0,
            control: k1,
            end: k2,
            rho,
        });
        if self.rng.chance(0.4) {
            self.con(Constraint::Rho {
                curve: k,
                value: 0.5,
            });
        }
        match self.rng.int(3) {
            0 => {
                let t = self.rng.range(0.15, 0.85);
                let q = self.point(at(t));
                self.con(Constraint::PointOnCurve { point: q, curve: k });
            }
            1 => {
                // A line on from the end, along the leg into it.
                let d = (p2 - p1).normalize();
                let ll = self.len();
                let (l, ls, _) = self.line(p2, p2 + d * ll);
                self.con(Constraint::Coincident { a: k2, b: ls });
                self.con(Constraint::Tangent {
                    a: k,
                    b: l,
                    kind: TangentKind::AtPoint { point: ls },
                });
            }
            _ => {
                // An arc ending at the start, tangent to the first leg
                // there: its centre on the normal to the leg, away from
                // the control point.
                let d = (p1 - p0).normalize();
                let side = if d.perp().dot(p2 - p0) > 0.0 {
                    -1.0
                } else {
                    1.0
                };
                let r = self.rng.range(0.3, 0.6);
                let c = p0 + d.perp() * (side * r);
                // The arc (counter-clockwise) ends or starts there.
                let a1 = (p0 - c).to_angle();
                let sweep = self.rng.range(0.6, 1.8);
                let ends_there = self.rng.chance(0.5);
                let a0 = if ends_there { a1 - sweep } else { a1 };
                let (ka, _, ks, ke) = self.arc(c, r, a0, sweep);
                let join = if ends_there { ke } else { ks };
                self.con(Constraint::Coincident { a: k0, b: join });
                self.con(Constraint::Tangent {
                    a: k,
                    b: ka,
                    kind: TangentKind::AtPoint { point: join },
                });
            }
        }
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
            // Unsigned horizontal and vertical distances only where the
            // perturbation cannot change the sign of the difference. The
            // same for aligned distances in both coordinates: other
            // constraints may pin one coordinate of a point, and then the
            // sign of the other difference picks one of two solutions (a
            // perturbation across it starts the solve beyond a fold).
            let kind = match self.rng.int(3) {
                0 if d.x.abs() > 0.3 && d.y.abs() > 0.3 => DistanceKind::Aligned,
                1 if d.x.abs() > 0.3 => DistanceKind::Horizontal,
                2 if d.y.abs() > 0.3 => DistanceKind::Vertical,
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
            self.con(Constraint::Length {
                line: l,
                value: 1.0,
            });
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
    generate_with(seed, motifs, scale, offset, false)
}

/// As [`generate`], with ellipses and elliptical arcs among the motifs
/// when `ellipses` (without, the sketches are those of [`generate`], so
/// existing seeds keep their sketches).
pub fn generate_with(
    seed: u64,
    motifs: usize,
    scale: f64,
    offset: DVec2,
    ellipses: bool,
) -> Generated {
    generate_kinds(seed, motifs, scale, offset, if ellipses { 10 } else { 9 })
}

/// As [`generate`], with ellipses, elliptical arcs and conics among the
/// motifs.
pub fn generate_with_conics(seed: u64, motifs: usize, scale: f64, offset: DVec2) -> Generated {
    generate_kinds(seed, motifs, scale, offset, 11)
}

/// The generator over the first `kinds` motifs.
fn generate_kinds(seed: u64, motifs: usize, scale: f64, offset: DVec2, kinds: usize) -> Generated {
    let mut g = G {
        b: B::new(),
        rng: Rng(seed),
        points: Vec::new(),
        lines: Vec::new(),
        circulars: Vec::new(),
        circles: Vec::new(),
    };
    for _ in 0..motifs {
        match g.rng.int(kinds) {
            9 => g.ellipses(),
            10 => g.conics(),
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
        }
    }
    Generated {
        sketch: s,
        points: g.points,
        lines: g.lines,
        circulars: g.circulars,
    }
}

/// Moves every point by up to `amount · scale` in each coordinate, scales
/// circle radii by up to `1 ± amount` and moves conics' ρ by up to
/// `± amount / 2`.
pub fn perturb(s: &Sketch<u32>, seed: u64, amount: f64, scale: f64) -> Sketch<u32> {
    let mut rng = Rng(seed ^ 0x5DEE_CE66_D1CE_4E5B);
    let mut out = s.clone();
    for e in out.entities.iter_mut() {
        match &mut e.geometry {
            Geometry::Point { at } => {
                *at += DVec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)) * (amount * scale);
            }
            Geometry::Circle { radius, .. } => *radius *= 1.0 + amount * rng.range(-1.0, 1.0),
            Geometry::Conic { rho, .. } => {
                *rho = (*rho + 0.5 * amount * rng.range(-1.0, 1.0)).clamp(0.01, 0.99);
            }
            _ => {}
        }
    }
    out
}
