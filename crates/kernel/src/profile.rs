//! Planar profiles: the kernel-neutral input of sweep operations, and their
//! normalization into oriented loops and regions.

use std::f64::consts::PI;

use glam::{DVec2, DVec3};
use openshape_nurbs::NurbsCurve;
use serde::{Deserialize, Serialize};

use crate::geom::Curve3;
use crate::math::Frame;
use crate::{KernelError, Tolerance};

/// One curve of a planar profile, in the profile plane's 2D coordinates.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ProfileCurve {
    /// Straight segment.
    Line {
        /// Start point.
        start: DVec2,
        /// End point.
        end: DVec2,
    },
    /// Circular arc. Positive `sweep_angle` runs counter-clockwise.
    Arc {
        /// Centre.
        center: DVec2,
        /// Radius.
        radius: f64,
        /// Angle of the start point (radians, from the +x axis).
        start_angle: f64,
        /// Signed angular extent (radians), `0 < |sweep| ≤ 2π`.
        sweep_angle: f64,
    },
    /// Full circle.
    Circle {
        /// Centre.
        center: DVec2,
        /// Radius.
        radius: f64,
    },
    /// Non-uniform rational B-spline.
    Spline {
        /// Degree (≥ 1).
        degree: usize,
        /// Knot vector, length `points.len() + degree + 1`.
        knots: Vec<f64>,
        /// Control points.
        points: Vec<DVec2>,
        /// Optional positive weights (all 1 when absent).
        weights: Option<Vec<f64>>,
    },
}

impl ProfileCurve {
    /// `n + 1` points at evenly spaced parameters along the curve (angle for
    /// arcs and circles, knot parameter for splines), starting and ending at
    /// its ends: for drawing and picking. Fails for an invalid spline.
    pub fn sample(&self, n: usize) -> Result<Vec<DVec2>, KernelError> {
        let n = n.max(1);
        let f = |i: usize| i as f64 / n as f64;
        Ok(match self {
            ProfileCurve::Line { start, end } => (0..=n).map(|i| start.lerp(*end, f(i))).collect(),
            ProfileCurve::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => (0..=n)
                .map(|i| *center + DVec2::from_angle(start_angle + sweep_angle * f(i)) * *radius)
                .collect(),
            ProfileCurve::Circle { center, radius } => (0..=n)
                .map(|i| *center + DVec2::from_angle(2.0 * PI * f(i)) * *radius)
                .collect(),
            ProfileCurve::Spline {
                degree,
                knots,
                points,
                weights,
            } => {
                let p3: Vec<DVec3> = points.iter().map(|p| p.extend(0.0)).collect();
                let c = NurbsCurve::new(*degree, knots.clone(), p3, weights.clone())
                    .map_err(|e| KernelError::InvalidInput(format!("spline: {e}")))?;
                let (a, b) = c.domain();
                (0..=n)
                    .map(|i| c.point(a + (b - a) * f(i)).truncate())
                    .collect()
            }
        })
    }

    /// The elliptical arc `center + cos θ·major + sin θ·minor`, θ from
    /// `start` to `start + sweep` (`0 < |sweep| ≤ 2π`), as an exact rational
    /// quadratic B-spline (design note 0008 §6).
    ///
    /// `major` and `minor` are conjugate semi-diameters (the semi-axes when
    /// perpendicular). The θ range is split into `k = ⌈|sweep|/(π/2)⌉`
    /// equal segments of `δ = sweep/k`; each is the rational quadratic Bézier
    /// arc of Piegl & Tiller A7.1 (end points `E(θᵢ)`, `E(θᵢ₊₁)`, middle
    /// point `center + (cos θₘ·major + sin θₘ·minor)/cos(δ/2)` of weight
    /// `cos(δ/2)`, `θₘ` the segment's middle), mapped by the affine map that
    /// takes the unit circle to the ellipse; rational curves are affine
    /// invariant (P&T §4.2), so the result is the ellipse exactly. Knots are
    /// `0,0,0, 1/k,1/k, …, 1,1,1`. A sweep of a full turn ends exactly on its
    /// start point.
    ///
    /// Fails for non-finite input, parallel (or zero) semi-diameters, and a
    /// sweep that is zero or beyond a full turn.
    pub fn elliptical_arc(
        center: DVec2,
        major: DVec2,
        minor: DVec2,
        start: f64,
        sweep: f64,
    ) -> Result<ProfileCurve, KernelError> {
        if !(center.is_finite() && major.is_finite() && minor.is_finite() && start.is_finite()) {
            return Err(KernelError::InvalidInput(
                "elliptical arc: non-finite input".into(),
            ));
        }
        let turn = major.perp_dot(minor);
        if turn == 0.0 || turn.is_nan() {
            return Err(KernelError::InvalidInput(
                "elliptical arc: the semi-diameters are parallel or zero".into(),
            ));
        }
        if !(sweep.abs() > 0.0 && sweep.abs() <= 2.0 * PI) {
            return Err(KernelError::InvalidInput(format!(
                "elliptical arc: sweep {sweep} is not within (0, 2π]"
            )));
        }
        let at = |t: f64| {
            let (s, c) = t.sin_cos();
            center + major * c + minor * s
        };
        let full = sweep.abs() == 2.0 * PI;
        let k = ((sweep.abs() / (0.5 * PI)).ceil() as usize).max(1);
        let delta = sweep / k as f64;
        let w = (0.5 * delta).cos();
        let mut points = Vec::with_capacity(2 * k + 1);
        let mut weights = Vec::with_capacity(2 * k + 1);
        for i in 0..k {
            let t0 = start + delta * i as f64;
            let tm = t0 + 0.5 * delta;
            points.push(at(t0));
            weights.push(1.0);
            let (s, c) = tm.sin_cos();
            points.push(center + (major * c + minor * s) / w);
            weights.push(w);
        }
        points.push(if full { points[0] } else { at(start + sweep) });
        weights.push(1.0);
        let mut knots = vec![0.0; 3];
        for i in 1..k {
            let u = i as f64 / k as f64;
            knots.extend([u, u]);
        }
        knots.extend([1.0; 3]);
        Ok(ProfileCurve::Spline {
            degree: 2,
            knots,
            points,
            weights: Some(weights),
        })
    }

    /// A cubic (or lower, for few points) spline interpolating `points`.
    pub fn interpolate(points: &[DVec2]) -> Result<ProfileCurve, KernelError> {
        let degree = 3.min(points.len().saturating_sub(1)).max(1);
        let p3: Vec<DVec3> = points.iter().map(|p| p.extend(0.0)).collect();
        let c = NurbsCurve::interpolate(&p3, degree)
            .map_err(|e| KernelError::InvalidInput(format!("spline: {e}")))?;
        Ok(ProfileCurve::Spline {
            degree: c.degree(),
            knots: c.knots().as_slice().to_vec(),
            points: c.control_points().iter().map(|p| p.truncate()).collect(),
            weights: None,
        })
    }
}

/// A planar profile: a set of curves forming one or more closed loops in a
/// plane. Loops may be nested (holes, islands) but must not intersect.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Profile {
    /// The profile plane. Curve coordinates are in its x/y axes.
    pub plane: Frame,
    /// The curves, in any order and orientation.
    pub curves: Vec<ProfileCurve>,
}

impl Profile {
    /// Axis-aligned rectangle with corners `min` and `max`.
    pub fn rectangle(plane: Frame, min: DVec2, max: DVec2) -> Profile {
        Self::polygon(
            plane,
            &[min, DVec2::new(max.x, min.y), max, DVec2::new(min.x, max.y)],
        )
    }

    /// Closed polygon through `points`.
    pub fn polygon(plane: Frame, points: &[DVec2]) -> Profile {
        let n = points.len();
        Profile {
            plane,
            curves: (0..n)
                .map(|i| ProfileCurve::Line {
                    start: points[i],
                    end: points[(i + 1) % n],
                })
                .collect(),
        }
    }

    /// A circle.
    pub fn circle(plane: Frame, center: DVec2, radius: f64) -> Profile {
        Profile {
            plane,
            curves: vec![ProfileCurve::Circle { center, radius }],
        }
    }

    /// Adds the curves of another profile in the same plane (for holes).
    pub fn with(mut self, other: Profile) -> Profile {
        self.curves.extend(other.curves);
        self
    }
}

/// An oriented 2D segment in profile-plane coordinates (`z = 0`), traversed
/// from `t0` to `t1` with `t0 < t1`.
#[derive(Clone, Debug)]
pub(crate) struct Seg {
    pub curve: Curve3,
    pub t0: f64,
    pub t1: f64,
}

impl Seg {
    pub fn start(&self) -> DVec2 {
        self.curve.point(self.t0).truncate()
    }
    pub fn end(&self) -> DVec2 {
        self.curve.point(self.t1).truncate()
    }
    pub fn point(&self, f: f64) -> DVec2 {
        self.curve
            .point(self.t0 + (self.t1 - self.t0) * f)
            .truncate()
    }
    pub fn reversed(&self) -> Seg {
        match &self.curve {
            Curve3::Line { origin, dir } => {
                let a = *origin + *dir * self.t0;
                let b = *origin + *dir * self.t1;
                Seg {
                    curve: Curve3::Line {
                        origin: b,
                        dir: a - b,
                    },
                    t0: 0.0,
                    t1: 1.0,
                }
            }
            Curve3::Circle {
                center,
                x,
                y,
                radius,
            } => {
                // Start at the old end and run the other way.
                let (s, c) = self.t1.sin_cos();
                let nx = *x * c + *y * s;
                let ny = -(*y * c - *x * s);
                Seg {
                    curve: Curve3::Circle {
                        center: *center,
                        x: nx,
                        y: ny,
                        radius: *radius,
                    },
                    t0: 0.0,
                    t1: self.t1 - self.t0,
                }
            }
            Curve3::Nurbs(n) => {
                let (a, b) = n.domain();
                let r = n.reversed();
                // reversed() maps t to a + b - t.
                Seg {
                    curve: Curve3::Nurbs(r),
                    t0: a + b - self.t1,
                    t1: a + b - self.t0,
                }
            }
        }
    }

    /// Polyline approximation (including both ends).
    pub fn sample(&self, n: usize) -> Vec<DVec2> {
        (0..=n).map(|i| self.point(i as f64 / n as f64)).collect()
    }
}

/// A closed, oriented loop of segments.
pub(crate) type SegLoop = Vec<Seg>;

/// A connected planar region: an outer loop (counter-clockwise about the
/// plane normal) followed by hole loops (clockwise).
#[derive(Clone, Debug)]
pub(crate) struct Region {
    pub loops: Vec<SegLoop>,
}

fn arc_segs(center: DVec2, radius: f64, start: f64, sweep: f64) -> Vec<Seg> {
    let pieces = ((sweep.abs() / PI) - 1e-9).ceil().max(1.0) as usize;
    let step = sweep / pieces as f64;
    (0..pieces)
        .map(|i| {
            let a = start + step * i as f64;
            let x = DVec3::new(a.cos(), a.sin(), 0.0);
            let y = DVec3::new(-a.sin(), a.cos(), 0.0) * step.signum();
            Seg {
                curve: Curve3::Circle {
                    center: center.extend(0.0),
                    x,
                    y,
                    radius,
                },
                t0: 0.0,
                t1: step.abs(),
            }
        })
        .collect()
}

fn to_segs(c: &ProfileCurve, tol: &Tolerance) -> Result<Vec<Seg>, KernelError> {
    let bad = |m: &str| Err(KernelError::InvalidInput(m.to_string()));
    match c {
        ProfileCurve::Line { start, end } => {
            if !(start.is_finite() && end.is_finite()) {
                return bad("line has non-finite coordinates");
            }
            if (*end - *start).length() <= tol.linear {
                return bad("line is shorter than the model tolerance");
            }
            Ok(vec![Seg {
                curve: Curve3::Line {
                    origin: start.extend(0.0),
                    dir: (*end - *start).extend(0.0),
                },
                t0: 0.0,
                t1: 1.0,
            }])
        }
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            if !(center.is_finite()
                && radius.is_finite()
                && start_angle.is_finite()
                && sweep_angle.is_finite())
            {
                return bad("arc has non-finite parameters");
            }
            if *radius <= tol.linear {
                return bad("arc radius must be positive");
            }
            if sweep_angle.abs() * radius <= tol.linear || sweep_angle.abs() > 2.0 * PI + 1e-12 {
                return bad("arc sweep must be in (0, 2π]");
            }
            Ok(arc_segs(*center, *radius, *start_angle, *sweep_angle))
        }
        ProfileCurve::Circle { center, radius } => {
            if !(center.is_finite() && radius.is_finite()) || *radius <= tol.linear {
                return bad("circle radius must be positive");
            }
            Ok(arc_segs(*center, *radius, 0.0, 2.0 * PI))
        }
        ProfileCurve::Spline {
            degree,
            knots,
            points,
            weights,
        } => {
            let p3: Vec<DVec3> = points.iter().map(|p| p.extend(0.0)).collect();
            let n = NurbsCurve::new(*degree, knots.clone(), p3, weights.clone())
                .map_err(|e| KernelError::InvalidInput(format!("spline: {e}")))?;
            let (a, b) = n.domain();
            let curve = Curve3::Nurbs(n);
            let closed = (curve.point(a) - curve.point(b)).length() <= tol.linear;
            if closed {
                // Split closed splines so no face is periodic in this direction.
                let m = 0.5 * (a + b);
                Ok(vec![
                    Seg {
                        curve: curve.clone(),
                        t0: a,
                        t1: m,
                    },
                    Seg {
                        curve,
                        t0: m,
                        t1: b,
                    },
                ])
            } else {
                Ok(vec![Seg {
                    curve,
                    t0: a,
                    t1: b,
                }])
            }
        }
    }
}

/// Signed area of a closed polyline (positive = counter-clockwise).
fn signed_area(poly: &[DVec2]) -> f64 {
    let n = poly.len();
    (0..n)
        .map(|i| poly[i].perp_dot(poly[(i + 1) % n]))
        .sum::<f64>()
        * 0.5
}

fn loop_polyline(l: &SegLoop) -> Vec<DVec2> {
    let mut out = Vec::new();
    for s in l {
        let pts = s.sample(48);
        out.extend_from_slice(&pts[..pts.len() - 1]);
    }
    out
}

fn point_in_polygon(p: DVec2, poly: &[DVec2]) -> bool {
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if p.x < x {
                inside = !inside;
            }
        }
    }
    inside
}

fn orient(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    )
}

/// True if the closed segments `ab` and `cd` properly intersect or touch.
fn segments_touch(a: DVec2, b: DVec2, c: DVec2, d: DVec2) -> bool {
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    if ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0))
        && ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))
    {
        return true;
    }
    let on = |p: DVec2, q: DVec2, r: DVec2| {
        r.x >= p.x.min(q.x) && r.x <= p.x.max(q.x) && r.y >= p.y.min(q.y) && r.y <= p.y.max(q.y)
    };
    (o1 == 0.0 && on(a, b, c))
        || (o2 == 0.0 && on(a, b, d))
        || (o3 == 0.0 && on(c, d, a))
        || (o4 == 0.0 && on(c, d, b))
}

/// Checks loops (as polylines) for self- and mutual intersections.
fn check_intersections(polys: &[Vec<DVec2>]) -> Result<(), KernelError> {
    let mut segs: Vec<(usize, usize, DVec2, DVec2)> = Vec::new();
    for (li, p) in polys.iter().enumerate() {
        for i in 0..p.len() {
            segs.push((li, i, p[i], p[(i + 1) % p.len()]));
        }
    }
    for i in 0..segs.len() {
        let (la, ia, a, b) = segs[i];
        let (bmin, bmax) = (a.min(b), a.max(b));
        for &(lb, ib, c, d) in &segs[i + 1..] {
            if la == lb {
                let n = polys[la].len();
                if ib == (ia + 1) % n || ia == (ib + 1) % n {
                    continue; // adjacent segments share an end point
                }
            }
            if c.max(d).cmplt(bmin).any() || c.min(d).cmpgt(bmax).any() {
                continue;
            }
            if segments_touch(a, b, c, d) {
                return Err(KernelError::InvalidInput(
                    "profile curves intersect or touch each other".into(),
                ));
            }
        }
    }
    Ok(())
}

/// Chains curves into closed loops, orients them and groups them into regions.
pub(crate) fn regions(profile: &Profile, tol: &Tolerance) -> Result<Vec<Region>, KernelError> {
    if !profile.plane.is_valid(1e-9) {
        return Err(KernelError::InvalidInput(
            "profile plane is not orthonormal".into(),
        ));
    }
    let mut segs: Vec<Seg> = Vec::new();
    for c in &profile.curves {
        segs.extend(to_segs(c, tol)?);
    }
    if segs.is_empty() {
        return Err(KernelError::InvalidInput("profile is empty".into()));
    }
    // Cluster end points into vertices.
    let snap = tol.linear * 10.0;
    let mut verts: Vec<DVec2> = Vec::new();
    let vid = |p: DVec2, verts: &mut Vec<DVec2>| -> usize {
        if let Some(i) = verts.iter().position(|v| (*v - p).length() <= snap) {
            i
        } else {
            verts.push(p);
            verts.len() - 1
        }
    };
    let ends: Vec<(usize, usize)> = segs
        .iter()
        .map(|s| (vid(s.start(), &mut verts), vid(s.end(), &mut verts)))
        .collect();
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); verts.len()];
    for (i, (a, b)) in ends.iter().enumerate() {
        if a == b {
            return Err(KernelError::InvalidInput("degenerate closed curve".into()));
        }
        incident[*a].push(i);
        incident[*b].push(i);
    }
    if let Some(v) = incident.iter().position(|inc| inc.len() != 2) {
        let n = incident[v].len();
        return Err(KernelError::InvalidInput(if n < 2 {
            format!("profile is open at ({:.6}, {:.6})", verts[v].x, verts[v].y)
        } else {
            format!("profile branches at ({:.6}, {:.6})", verts[v].x, verts[v].y)
        }));
    }
    // Walk loops.
    let mut used = vec![false; segs.len()];
    let mut loops: Vec<SegLoop> = Vec::new();
    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        let mut lp: SegLoop = Vec::new();
        let mut cur = start;
        let mut at = ends[start].0; // vertex we are leaving from
        loop {
            used[cur] = true;
            let (a, b) = ends[cur];
            let (seg, next_at) = if a == at {
                (segs[cur].clone(), b)
            } else {
                (segs[cur].reversed(), a)
            };
            lp.push(seg);
            at = next_at;
            let next = incident[at].iter().copied().find(|&s| s != cur && !used[s]);
            match next {
                Some(n) => cur = n,
                None => break,
            }
        }
        loops.push(lp);
    }
    // Snap segment ends exactly onto shared vertices is left to the B-rep
    // builder (vertices are shared); here we only orient and nest.
    let polys: Vec<Vec<DVec2>> = loops.iter().map(loop_polyline).collect();
    check_intersections(&polys)?;
    let areas: Vec<f64> = polys.iter().map(|p| signed_area(p)).collect();
    if areas.iter().any(|a| a.abs() <= tol.linear * tol.linear) {
        return Err(KernelError::InvalidInput(
            "profile loop has zero area".into(),
        ));
    }
    // Depth = number of loops containing this loop.
    let n = loops.len();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut depth = vec![0usize; n];
    for i in 0..n {
        let probe = polys[i][0];
        let mut best: Option<usize> = None;
        for j in 0..n {
            if i != j && point_in_polygon(probe, &polys[j]) {
                depth[i] += 1;
                if best.is_none_or(|b| areas[j].abs() < areas[b].abs()) {
                    best = Some(j);
                }
            }
        }
        parent[i] = best;
    }
    let mut out_regions: Vec<Region> = Vec::new();
    let mut region_of: Vec<Option<usize>> = vec![None; n];
    let orient_loop = |l: &SegLoop, ccw: bool, area: f64| -> SegLoop {
        if (area > 0.0) == ccw {
            l.clone()
        } else {
            l.iter().rev().map(|s| s.reversed()).collect()
        }
    };
    for i in 0..n {
        if depth[i] % 2 == 0 {
            region_of[i] = Some(out_regions.len());
            out_regions.push(Region {
                loops: vec![orient_loop(&loops[i], true, areas[i])],
            });
        }
    }
    for i in 0..n {
        if depth[i] % 2 == 1 {
            let p = parent[i].expect("odd depth implies a parent");
            let r = region_of[p].expect("parent of a hole is an outer loop");
            out_regions[r]
                .loops
                .push(orient_loop(&loops[i], false, areas[i]));
        }
    }
    Ok(out_regions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangle_with_hole_gives_one_region() {
        let p = Profile::rectangle(Frame::WORLD, DVec2::ZERO, DVec2::new(4.0, 3.0))
            .with(Profile::circle(Frame::WORLD, DVec2::new(2.0, 1.5), 0.5));
        let r = regions(&p, &Tolerance::default()).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].loops.len(), 2);
        assert!(signed_area(&loop_polyline(&r[0].loops[0])) > 0.0);
        assert!(signed_area(&loop_polyline(&r[0].loops[1])) < 0.0);
    }

    #[test]
    fn reversed_and_shuffled_curves_chain() {
        let p = Profile {
            plane: Frame::WORLD,
            curves: vec![
                ProfileCurve::Line {
                    start: DVec2::new(1.0, 0.0),
                    end: DVec2::ZERO,
                },
                ProfileCurve::Line {
                    start: DVec2::new(0.0, 1.0),
                    end: DVec2::new(1.0, 1.0),
                },
                ProfileCurve::Line {
                    start: DVec2::new(1.0, 0.0),
                    end: DVec2::new(1.0, 1.0),
                },
                ProfileCurve::Line {
                    start: DVec2::ZERO,
                    end: DVec2::new(0.0, 1.0),
                },
            ],
        };
        let r = regions(&p, &Tolerance::default()).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].loops[0].len(), 4);
        let l = &r[0].loops[0];
        for i in 0..4 {
            assert!((l[i].end() - l[(i + 1) % 4].start()).length() < 1e-12);
        }
    }

    #[test]
    fn open_and_crossing_profiles_are_rejected() {
        let open = Profile {
            plane: Frame::WORLD,
            curves: vec![ProfileCurve::Line {
                start: DVec2::ZERO,
                end: DVec2::X,
            }],
        };
        assert!(regions(&open, &Tolerance::default()).is_err());
        // Bow-tie.
        let bow = Profile::polygon(
            Frame::WORLD,
            &[
                DVec2::ZERO,
                DVec2::new(1.0, 1.0),
                DVec2::new(1.0, 0.0),
                DVec2::new(0.0, 1.0),
            ],
        );
        assert!(regions(&bow, &Tolerance::default()).is_err());
        // Two overlapping circles.
        let two = Profile::circle(Frame::WORLD, DVec2::ZERO, 1.0).with(Profile::circle(
            Frame::WORLD,
            DVec2::new(1.0, 0.0),
            1.0,
        ));
        assert!(regions(&two, &Tolerance::default()).is_err());
    }

    #[test]
    fn arc_reversal_preserves_geometry() {
        let s = &arc_segs(DVec2::new(1.0, 2.0), 3.0, 0.3, 1.2)[0];
        let r = s.reversed();
        assert!((r.start() - s.end()).length() < 1e-12);
        assert!((r.end() - s.start()).length() < 1e-12);
        assert!((r.point(0.25) - s.point(0.75)).length() < 1e-12);
    }
}
