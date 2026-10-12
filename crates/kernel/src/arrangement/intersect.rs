//! Contacts between curves: transversal crossings, tangencies, end points
//! lying on other curves, and the ends of overlaps (design note 0006, §4.1–2).
//!
//! * Line–line, line–circle and circle–circle pairs use closed forms
//!   (arcs as their circles, filtered by parameter range).
//! * Pairs involving a spline use recursive subdivision of rational Bézier
//!   pieces with separating-axis exclusion, in the spirit of Bézier clipping
//!   (Sederberg & Nishita, "Curve intersection using Bézier clipping",
//!   *CAD* 22(9), 1990; Patrikalakis & Maekawa, *Shape Interrogation for
//!   CAD/CAM*, 2002, §5.6–5.8): once both pieces are flat to an eighth of the
//!   tolerance, their chords give starting points that are refined on the
//!   exact curves by Levenberg–Marquardt iteration.
//! * Self-intersections of a spline: the spline is cut into pieces that are
//!   monotone along some direction (which cannot self-intersect), and the
//!   pieces are intersected pairwise.
//!
//! Every reported contact `(ta, tb)` satisfies `|A(ta) − B(tb)| ≤ tol`.
//! Collinear or co-circular overlaps are reported by their two ends only;
//! the duplicated piece between them is removed later, when pieces joining
//! the same two vertices are found to lie within the tolerance of each other.

use std::collections::HashMap;

use glam::DVec2;

use super::bezier;
use super::curve::{BPiece, Crv, Geom};
use crate::KernelError;

/// Shared state: the tolerance and a work budget that bounds the
/// subdivision on pathological inputs.
pub(super) struct Ctx {
    /// Linear tolerance.
    pub tol: f64,
    /// Remaining subdivision steps.
    pub budget: u64,
}

impl Ctx {
    fn spend(&mut self) -> Result<(), KernelError> {
        if self.budget == 0 {
            return Err(KernelError::InvalidInput(
                "curves are too complex to intersect reliably".into(),
            ));
        }
        self.budget -= 1;
        Ok(())
    }
}

/// Contacts between two different curves.
pub(super) fn pair_contacts(
    a: &Crv,
    b: &Crv,
    ctx: &mut Ctx,
) -> Result<Vec<(f64, f64)>, KernelError> {
    let tol = ctx.tol;
    let mut out = Vec::new();
    endpoint_contacts(a, b, tol, &mut out);
    match (&a.geom, &b.geom) {
        (Geom::Line { .. }, Geom::Line { .. }) => line_line(a, b, tol, &mut out),
        (Geom::Line { .. }, Geom::Arc { c, r, .. }) => {
            for x in line_circle(a, *c, *r, tol) {
                let s = a.closest(x, 0.0, 1.0).0;
                let t = b.closest(a.point(s), b.lo, b.hi).0;
                out.push((s, t));
            }
        }
        (Geom::Arc { c, r, .. }, Geom::Line { .. }) => {
            for x in line_circle(b, *c, *r, tol) {
                let t = b.closest(x, 0.0, 1.0).0;
                let s = a.closest(b.point(t), a.lo, a.hi).0;
                out.push((s, t));
            }
        }
        (Geom::Arc { c: ca, r: ra, .. }, Geom::Arc { c: cb, r: rb, .. }) => {
            if (*ca - *cb).length() + (ra - rb).abs() <= tol {
                // Co-circular: the end points (already recorded) bound the
                // overlap; two full circles share a vertex at angle 0.
                if a.periodic && b.periodic {
                    let p = a.point(0.0);
                    out.push((0.0, b.closest(p, b.lo, b.hi).0));
                }
            } else {
                for x in circle_circle(*ca, *ra, *cb, *rb, tol) {
                    out.push((a.closest(x, a.lo, a.hi).0, b.closest(x, b.lo, b.hi).0));
                }
            }
        }
        _ => bez_pair(a, b, ctx, &mut out)?,
    }
    out.retain(|&(s, t)| {
        s.is_finite() && t.is_finite() && (a.point(s) - b.point(t)).length() <= tol
    });
    Ok(out)
}

/// End points of either curve that lie within the tolerance of the other.
fn endpoint_contacts(a: &Crv, b: &Crv, tol: f64, out: &mut Vec<(f64, f64)>) {
    if !a.periodic {
        for ta in [a.lo, a.hi] {
            let (tb, d) = b.closest(a.point(ta), b.lo, b.hi);
            if d <= tol {
                out.push((ta, tb));
            }
        }
    }
    if !b.periodic {
        for tb in [b.lo, b.hi] {
            let (ta, d) = a.closest(b.point(tb), a.lo, a.hi);
            if d <= tol {
                out.push((ta, tb));
            }
        }
    }
}

/// Line–line: the crossing of the two segments; collinear overlaps are
/// reduced to their two extreme contacts (by convexity of the distance
/// between two segments, everything between two contacts is within the
/// tolerance).
fn line_line(a: &Crv, b: &Crv, tol: f64, out: &mut Vec<(f64, f64)>) {
    let (a0, a1) = (a.point(0.0), a.point(1.0));
    let (b0, b1) = (b.point(0.0), b.point(1.0));
    let da = a1 - a0;
    let db = b1 - b0;
    let den = da.perp_dot(db);
    if den != 0.0 {
        let w = b0 - a0;
        let s = w.perp_dot(db) / den;
        let t = w.perp_dot(da) / den;
        let ea = tol / da.length();
        let eb = tol / db.length();
        if s.is_finite()
            && t.is_finite()
            && (-ea..=1.0 + ea).contains(&s)
            && (-eb..=1.0 + eb).contains(&t)
        {
            out.push((s.clamp(0.0, 1.0), t.clamp(0.0, 1.0)));
        }
    }
    let valid: Vec<(f64, f64)> = out
        .iter()
        .copied()
        .filter(|&(s, t)| (a.point(s) - b.point(t)).length() <= tol)
        .collect();
    if valid.len() >= 2 {
        let lo = valid
            .iter()
            .copied()
            .min_by(|x, y| x.1.total_cmp(&y.1))
            .unwrap_or(valid[0]);
        let hi = valid
            .iter()
            .copied()
            .max_by(|x, y| x.1.total_cmp(&y.1))
            .unwrap_or(valid[0]);
        if (b.point(hi.1) - b.point(lo.1)).length() > tol {
            *out = vec![lo, hi];
        }
    }
}

/// Points of the infinite line through segment `l` that lie on (or, for a
/// tangency, within the tolerance of) the circle `(c, r)`.
fn line_circle(l: &Crv, c: DVec2, r: f64, tol: f64) -> Vec<DVec2> {
    let a0 = l.point(0.0);
    let da = l.point(1.0) - a0;
    let len = da.length();
    let u = da / len;
    let f = a0 + u * (c - a0).dot(u);
    let h = (c - a0).perp_dot(u).abs();
    if h > r + tol {
        Vec::new()
    } else if h >= r {
        vec![f]
    } else {
        let q = ((r - h) * (r + h)).sqrt();
        vec![f - u * q, f + u * q]
    }
}

/// Points common to two circles; for a tangency within the tolerance, the
/// point of the first circle on the centre line.
fn circle_circle(ca: DVec2, ra: f64, cb: DVec2, rb: f64, tol: f64) -> Vec<DVec2> {
    let d = (cb - ca).length();
    if d == 0.0 || d > ra + rb + tol || d < (ra - rb).abs() - tol {
        return Vec::new();
    }
    let u = (cb - ca) / d;
    let a = (d * d + ra * ra - rb * rb) / (2.0 * d);
    let h2 = (ra - a) * (ra + a);
    if h2 <= 0.0 {
        vec![ca + u * (ra * a.signum())]
    } else {
        let h = h2.sqrt();
        vec![ca + u * a - u.perp() * h, ca + u * a + u.perp() * h]
    }
}

/// A node of the subdivision: a rational Bézier piece with cached bounds.
struct Node {
    piece: BPiece,
    hull: Vec<DVec2>,
    lo: DVec2,
    hi: DVec2,
    flat: f64,
    diam: f64,
}

impl Node {
    fn new(piece: BPiece) -> Node {
        let hull = piece.bez.hull();
        let (lo, hi) = bezier::bounds(&hull);
        let flat = bezier::flatness(&hull);
        Node {
            piece,
            hull,
            lo,
            hi,
            flat,
            diam: (hi - lo).length(),
        }
    }

    fn split(&self, crv: &Crv) -> (Node, Node) {
        let (l, r) = self.piece.bez.split(0.5);
        let tm = self.piece.param_at(crv, 0.5);
        (
            Node::new(BPiece {
                bez: l,
                t0: self.piece.t0,
                t1: tm,
                linear: self.piece.linear,
            }),
            Node::new(BPiece {
                bez: r,
                t0: tm,
                t1: self.piece.t1,
                linear: self.piece.linear,
            }),
        )
    }
}

/// True if the convex hulls are more than `tol` apart along one of the
/// coordinate axes or the chord normals (separating axis theorem).
fn separated(p: &Node, q: &Node, tol: f64) -> bool {
    if p.lo.x > q.hi.x + tol || q.lo.x > p.hi.x + tol {
        return true;
    }
    if p.lo.y > q.hi.y + tol || q.lo.y > p.hi.y + tol {
        return true;
    }
    for n in [chord_normal(&p.hull), chord_normal(&q.hull)]
        .into_iter()
        .flatten()
    {
        let (pl, ph) = extent(&p.hull, n);
        let (ql, qh) = extent(&q.hull, n);
        if pl > qh + tol || ql > ph + tol {
            return true;
        }
    }
    false
}

fn chord_normal(h: &[DVec2]) -> Option<DVec2> {
    (h[h.len() - 1] - h[0]).try_normalize().map(|d| d.perp())
}

fn extent(h: &[DVec2], n: DVec2) -> (f64, f64) {
    h.iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let x = p.dot(n);
            (lo.min(x), hi.max(x))
        })
}

/// Contacts between two curves of which at least one is a spline.
fn bez_pair(a: &Crv, b: &Crv, ctx: &mut Ctx, out: &mut Vec<(f64, f64)>) -> Result<(), KernelError> {
    let na: Vec<Node> = a.pieces(a.lo, a.hi).into_iter().map(Node::new).collect();
    let nb: Vec<Node> = b.pieces(b.lo, b.hi).into_iter().map(Node::new).collect();
    let mut found = Found::new(ctx.tol);
    for p in &na {
        for q in &nb {
            let mut s = Subdiv {
                a,
                b,
                ctx,
                bounds: None,
                out: &mut found,
            };
            s.rec(p, q, 0)?;
        }
    }
    out.extend(found.list.iter().map(|f| (f.0, f.1)));
    Ok(())
}

/// Contacts found so far, with a grid over their points so that duplicates
/// (both points within a quarter of the tolerance of an earlier contact's)
/// are found in constant time; along an overlap there are many contacts.
struct Found {
    list: Vec<(f64, f64, DVec2, DVec2)>,
    grid: HashMap<(i64, i64), Vec<usize>>,
    cell: f64,
}

impl Found {
    fn new(tol: f64) -> Found {
        Found {
            list: Vec::new(),
            grid: HashMap::new(),
            cell: 0.25 * tol,
        }
    }

    fn key(&self, p: DVec2) -> (i64, i64) {
        (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
        )
    }

    /// Records `(s, t)` with points `pa`, `pb` unless it duplicates an
    /// earlier contact.
    fn insert(&mut self, s: f64, t: f64, pa: DVec2, pb: DVec2) {
        let (kx, ky) = self.key(pa);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(ids) = self.grid.get(&(kx + dx, ky + dy)) {
                    for &i in ids {
                        let (_, _, qa, qb) = self.list[i];
                        if (qa - pa).length() <= self.cell && (qb - pb).length() <= self.cell {
                            return;
                        }
                    }
                }
            }
        }
        self.grid.entry((kx, ky)).or_default().push(self.list.len());
        self.list.push((s, t, pa, pb));
    }
}

/// Recursive subdivision between pieces of `a` and `b`.
struct Subdiv<'a, 'c> {
    a: &'a Crv,
    b: &'a Crv,
    ctx: &'c mut Ctx,
    /// Parameter bounds for refinement (`a`, then `b`), when restricted.
    bounds: Option<[(f64, f64); 2]>,
    out: &'c mut Found,
}

impl Subdiv<'_, '_> {
    fn rec(&mut self, p: &Node, q: &Node, depth: u32) -> Result<(), KernelError> {
        self.ctx.spend()?;
        let tol = self.ctx.tol;
        if separated(p, q, tol) {
            return Ok(());
        }
        let eps = tol / 8.0;
        let pf = p.flat <= eps;
        let qf = q.flat <= eps;
        if (pf && qf) || depth >= 64 || (p.diam <= tol / 4.0 && q.diam <= tol / 4.0) {
            self.leaf(p, q);
            return Ok(());
        }
        let split_p = !pf && (qf || p.diam >= q.diam);
        if split_p {
            let (l, r) = p.split(self.a);
            self.rec(&l, q, depth + 1)?;
            self.rec(&r, q, depth + 1)
        } else {
            let (l, r) = q.split(self.b);
            self.rec(p, &l, depth + 1)?;
            self.rec(p, &r, depth + 1)
        }
    }

    /// Both pieces are (nearly) straight: candidate contacts from their
    /// chords, refined on the exact curves.
    fn leaf(&mut self, p: &Node, q: &Node) {
        let tol = self.ctx.tol;
        let (a0, a1) = (p.piece.bez.start(), p.piece.bez.end());
        let (b0, b1) = (q.piece.bez.start(), q.piece.bez.end());
        let reach = tol + p.flat + q.flat;
        let mut cand: Vec<(f64, f64)> = Vec::new();
        let da = a1 - a0;
        let db = b1 - b0;
        let den = da.perp_dot(db);
        if den != 0.0 {
            let w = b0 - a0;
            let s = w.perp_dot(db) / den;
            let t = w.perp_dot(da) / den;
            if s.is_finite()
                && t.is_finite()
                && (-0.01..=1.01).contains(&s)
                && (-0.01..=1.01).contains(&t)
            {
                cand.push((s.clamp(0.0, 1.0), t.clamp(0.0, 1.0)));
            }
        }
        let proj = |x: DVec2, o: DVec2, d: DVec2| {
            let l2 = d.length_squared();
            if l2 == 0.0 {
                0.0
            } else {
                ((x - o).dot(d) / l2).clamp(0.0, 1.0)
            }
        };
        for (ua, x) in [(0.0, a0), (1.0, a1)] {
            let ub = proj(x, b0, db);
            if (b0 + db * ub - x).length() <= reach {
                cand.push((ua, ub));
            }
        }
        for (ub, x) in [(0.0, b0), (1.0, b1)] {
            let ua = proj(x, a0, da);
            if (a0 + da * ua - x).length() <= reach {
                cand.push((ua, ub));
            }
        }
        for (ua, ub) in cand {
            let s = p.piece.param_at(self.a, ua);
            let t = q.piece.param_at(self.b, ub);
            let (s, t) = refine(self.a, self.b, s, t, self.bounds);
            let (s, t) = touch(self.a, self.b, s, t, self.bounds, tol);
            let pa = self.a.point(s);
            let pb = self.b.point(t);
            if (pa - pb).length() > tol {
                continue;
            }
            self.out.insert(s, t, pa, pb);
        }
    }
}

/// Levenberg–Marquardt refinement of `A(s) = B(t)` from `(s, t)`: converges
/// quadratically to transversal crossings and to the closest approach at
/// tangencies (Nocedal & Wright, *Numerical Optimization*, §10.3).
pub(super) fn refine(
    a: &Crv,
    b: &Crv,
    s0: f64,
    t0: f64,
    bounds: Option<[(f64, f64); 2]>,
) -> (f64, f64) {
    let clamp = |c: &Crv, t: f64, k: usize| match bounds {
        Some(bd) => t.clamp(bd[k].0, bd[k].1),
        None => c.clamp(t),
    };
    let mut s = clamp(a, s0, 0);
    let mut t = clamp(b, t0, 1);
    let mut f = a.point(s) - b.point(t);
    let mut f2 = f.length_squared();
    let mut lam = 1e-3;
    for _ in 0..64 {
        if f2 == 0.0 {
            break;
        }
        let da = a.eval(s, true)[1];
        let db = b.eval(t, true)[1];
        let m00 = da.dot(da);
        let m11 = db.dot(db);
        let m01 = -da.dot(db);
        let g0 = da.dot(f);
        let g1 = -db.dot(f);
        let floor = 1e-30 * (m00 + m11).max(f64::MIN_POSITIVE);
        let mut accepted = None;
        for _ in 0..24 {
            let a00 = m00 + lam * m00.max(floor);
            let a11 = m11 + lam * m11.max(floor);
            let det = a00 * a11 - m01 * m01;
            if !(det.is_finite() && det > 0.0) {
                lam *= 10.0;
                continue;
            }
            let ds = -(a11 * g0 - m01 * g1) / det;
            let dt = -(a00 * g1 - m01 * g0) / det;
            if !(ds.is_finite() && dt.is_finite()) {
                lam *= 10.0;
                continue;
            }
            let sn = clamp(a, s + ds, 0);
            let tn = clamp(b, t + dt, 1);
            let fnew = a.point(sn) - b.point(tn);
            let f2n = fnew.length_squared();
            if f2n < f2 {
                accepted = Some((sn, tn, fnew, f2n));
                lam = (lam * 0.3).max(1e-15);
                break;
            }
            lam *= 10.0;
        }
        let Some((sn, tn, fnew, f2n)) = accepted else {
            break;
        };
        let moved = (sn - s).abs() + (tn - t).abs();
        s = sn;
        t = tn;
        f = fnew;
        f2 = f2n;
        if moved <= 1e-16 * (1.0 + s.abs() + t.abs()) {
            break;
        }
    }
    (s, t)
}

/// Where [`refine`] stopped short of a tangency, the point of contact
/// itself.
///
/// Where two curves touch, the Jacobian of `A(s) − B(t)` has rank one and
/// the residual left is normal to both curves, which first-order steps
/// cannot reduce: refinement stops wherever the gap is too small to shrink
/// further, which can be microns along the curves from the touch (the
/// curves stay within ε of each other for about √(2Rε) either side, R the
/// radius of curvature). The point of contact, or of closest approach, is a
/// root of
///
///   G(s, t) = ((A − B)·B′, A′ × B′)
///
/// (A(s)'s foot on B, where the tangents are parallel), and a simple one:
/// at a touch with A′ = αT and B′ = βT, det ∂G/∂(s, t) = α²β³(κ_B − κ_A)
/// for the signed curvatures in the frame (T, N), non-zero unless the
/// curves osculate. Newton's method on G therefore converges to it
/// quadratically. Its closest iterate replaces refinement's result when it
/// is a closer approach and every iterate stays within the tolerance; a
/// converged crossing (residual at rounding level) or a miss (beyond the
/// tolerance) is returned unchanged.
fn touch(
    a: &Crv,
    b: &Crv,
    s0: f64,
    t0: f64,
    bounds: Option<[(f64, f64); 2]>,
    tol: f64,
) -> (f64, f64) {
    let clamp = |c: &Crv, t: f64, k: usize| match bounds {
        Some(bd) => t.clamp(bd[k].0, bd[k].1),
        None => c.clamp(t),
    };
    let gap = |s: f64, t: f64| (a.point(s) - b.point(t)).length();
    let (pa, pb) = (a.point(s0), b.point(t0));
    let noise = 64.0 * f64::EPSILON * pa.abs().max_element().max(pb.abs().max_element());
    let g0 = (pa - pb).length();
    if !(g0 > noise && g0 <= tol) {
        return (s0, t0);
    }
    let mut best = (s0, t0, g0);
    let (mut s, mut t) = (s0, t0);
    for _ in 0..32 {
        let [pa, da, dda] = a.eval(s, true);
        let [pb, db, ddb] = b.eval(t, true);
        let f = pa - pb;
        let g = [f.dot(db), da.perp_dot(db)];
        let j = [
            [da.dot(db), f.dot(ddb) - db.dot(db)],
            [dda.perp_dot(db), da.perp_dot(ddb)],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let ds = -(j[1][1] * g[0] - j[0][1] * g[1]) / det;
        let dt = -(j[0][0] * g[1] - j[1][0] * g[0]) / det;
        if !(ds.is_finite() && dt.is_finite()) {
            break;
        }
        let (sn, tn) = (clamp(a, s + ds, 0), clamp(b, t + dt, 1));
        let d = gap(sn, tn);
        if d > tol {
            return (s0, t0);
        }
        if d < best.2 {
            best = (sn, tn, d);
        }
        let moved = (sn - s).abs() + (tn - t).abs();
        (s, t) = (sn, tn);
        if moved <= 1e-15 * (1.0 + s.abs() + t.abs()) {
            break;
        }
    }
    (best.0, best.1)
}

/// Self-contacts of a spline: pairs `(s, t)`, `s < t`, of distinct points of
/// the curve that coincide within the tolerance and are joined by a genuine
/// loop of the curve.
pub(super) fn self_contacts(a: &Crv, ctx: &mut Ctx) -> Result<Vec<(f64, f64)>, KernelError> {
    if !matches!(a.geom, Geom::Spline { .. }) {
        return Ok(Vec::new());
    }
    let tol = ctx.tol;
    let nodes: Vec<Node> = a.monotone_pieces(tol).into_iter().map(Node::new).collect();
    let mut out: Vec<(f64, f64)> = Vec::new();
    for i in 0..nodes.len() {
        for j in i + 1..nodes.len() {
            let (p, q) = (&nodes[i], &nodes[j]);
            if separated(p, q, tol) {
                continue;
            }
            let adjacent = p.piece.t1 == q.piece.t0;
            let mut found = Found::new(tol);
            let mut s = Subdiv {
                a,
                b: a,
                ctx,
                bounds: Some([(p.piece.t0, p.piece.t1), (q.piece.t0, q.piece.t1)]),
                out: &mut found,
            };
            s.rec(p, q, 0)?;
            let junction = a.point(p.piece.t1);
            for (s, t, ps, pt) in found.list {
                if adjacent
                    && (ps - junction).length() <= 2.0 * tol
                    && (pt - junction).length() <= 2.0 * tol
                {
                    continue;
                }
                let (s, t) = if s <= t { (s, t) } else { (t, s) };
                if t <= s || !leaves(a, s, t, ps, 2.0 * tol) {
                    continue;
                }
                out.push((s, t));
            }
        }
    }
    Ok(out)
}

/// True if the curve between `s` and `t` gets farther than `r` from `p`.
fn leaves(a: &Crv, s: f64, t: f64, p: DVec2, r: f64) -> bool {
    a.pieces(s, t)
        .iter()
        .any(|pc| (0..=16).any(|k| (pc.bez.point(k as f64 / 16.0) - p).length() > r))
}
