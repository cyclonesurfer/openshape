//! Building the arrangement from curve contacts: vertices, edges, faces and
//! regions (design note 0006, §4.3–4.8). The data structure is a half-edge
//! map in the sense of de Berg, Cheong, van Kreveld & Overmars,
//! *Computational Geometry*, 3rd ed., ch. 2 (doubly-connected edge list,
//! face tracing, assignment of hole boundaries to faces).
//!
//! Robustness strategy, in order:
//! 1. Contacts are computed on the exact curves (see `intersect`).
//! 2. Contact and end points closer than the tolerance are merged into one
//!    vertex (union–find), and every curve passing within the tolerance of a
//!    vertex is routed through it (the "hot pixel" rule of snap rounding,
//!    Hobby 1999; Guibas & Marimont 1998), repeated to a fixed point. The
//!    vertex is placed at its most precise members (curve end points, then
//!    transversal crossings, then other contacts), and each curve's pieces
//!    end at its most precise incidence there, so pieces meeting at a
//!    vertex join to rounding error wherever the input curves do.
//! 3. Pieces of different curves (or of one curve) that join the same two
//!    vertices and lie within the tolerance of each other are overlaps; the
//!    piece of the lower curve index (then lower parameter) is kept. Vertices
//!    that only split one curve inside such an overlap are removed again.
//! 4. The result is checked: Euler's formula per connected component, exactly
//!    one outer boundary per component, positive region areas, an interior
//!    point at least the tolerance from each region's boundary, and no
//!    region's interior point inside another region. Any failure is an
//!    error, never a silently wrong region.

use std::collections::{BTreeMap, HashMap};
use std::f64::consts::{PI, TAU};

use glam::DVec2;

use super::curve::{Crv, Geom};
use super::intersect::{self, Ctx};
use super::measure::{self, RegionShape, Side};
use super::{ArrangementEdge, ArrangementVertex, Region, RegionLoop, SketchArrangement};
use crate::{KernelError, ProfileCurve};

/// Half-edges whose tangents differ by less than this angle (radians) are
/// ordered by their actual separation at a distance rather than by angle.
const TIE_ANGLE: f64 = 1e-3;
/// Total subdivision steps allowed for one arrangement.
pub(super) const BUDGET: u64 = 40_000_000;
/// Contacts whose tangents make at least this sine are transversal: their
/// position is known to about the rounding error divided by the sine.
const TRANSVERSAL_SINE: f64 = 1e-3;

/// How precisely an incidence locates its vertex, best first.
mod rank {
    /// A curve's own end point or corner: exact data.
    pub const EXACT: u8 = 0;
    /// A transversal crossing of two curves.
    pub const CROSSING: u8 = 1;
    /// A tangency, an overlap or a near-miss: placed along the curves only
    /// to about the tolerance.
    pub const CONTACT: u8 = 2;
    /// A curve routed through a vertex it passes within the tolerance of.
    pub const SNAPPED: u8 = 3;
}

fn inconsistent(m: impl std::fmt::Display) -> KernelError {
    KernelError::InvalidResult(format!("sketch arrangement is inconsistent: {m}"))
}

fn near_degenerate(p: DVec2, m: impl std::fmt::Display) -> KernelError {
    KernelError::InvalidInput(format!(
        "{m} near ({:.9}, {:.9}); curves nearly coincide or touch within the model tolerance",
        p.x, p.y
    ))
}

/// Union–find with path halving.
struct Uf {
    p: Vec<usize>,
}

impl Uf {
    fn add(&mut self) -> usize {
        self.p.push(self.p.len());
        self.p.len() - 1
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.p[x] != x {
            self.p[x] = self.p[self.p[x]];
            x = self.p[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            let (lo, hi) = (ra.min(rb), ra.max(rb));
            self.p[hi] = lo;
        }
    }
}

/// A point of a curve that is a vertex of the arrangement.
#[derive(Clone, Debug)]
struct Inc {
    curve: usize,
    t: f64,
    p: DVec2,
    /// True for a domain end of an open curve.
    end: bool,
    /// How precisely it locates its vertex (see [`rank`]).
    rank: u8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum State {
    Alive,
    /// Overlaps a kept piece.
    Duplicate,
    /// Absorbed into the preceding piece of the same curve.
    Merged,
    /// Dangling.
    Pruned,
}

#[derive(Clone, Debug)]
struct Piece {
    curve: usize,
    t0: f64,
    t1: f64,
    v0: usize,
    v1: usize,
    state: State,
    /// The next piece along the same curve (`usize::MAX` for none); its
    /// `t0` is this piece's `t1` (modulo 2π for a periodic curve).
    next: usize,
}

struct Vtx {
    /// Mean of the members: the centre used for routing curves and
    /// ordering edges.
    pos: DVec2,
    /// Mean of the most precise members (see [`rank`]): the output point,
    /// and where the pieces of each curve through the vertex end.
    point: DVec2,
    members: Vec<usize>,
    end: bool,
    removed: bool,
    /// The parameter at which each curve's pieces end at this vertex (the
    /// smallest where the curve passes more than once).
    params: BTreeMap<usize, f64>,
}

/// A tangent-tied half-edge being ordered around a vertex.
#[derive(Clone, Copy)]
struct Tied {
    he: usize,
    dir: DVec2,
    kappa: f64,
}

struct Builder<'a> {
    crvs: &'a [Crv],
    tol: f64,
    incs: Vec<Inc>,
    uf: Uf,
    /// Monotone pieces of each spline: `(t0, t1, lo, hi)`.
    mono: Vec<Vec<(f64, f64, DVec2, DVec2)>>,
}

/// Builds the arrangement of `curves` (design note 0006).
pub(super) fn build(curves: &[ProfileCurve], tol: f64) -> Result<SketchArrangement, KernelError> {
    let crvs: Vec<Crv> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| Crv::new(i, c, tol))
        .collect::<Result<_, _>>()?;
    if crvs.is_empty() {
        return Ok(SketchArrangement::default());
    }
    let mut b = Builder {
        crvs: &crvs,
        tol,
        incs: Vec::new(),
        uf: Uf { p: Vec::new() },
        mono: Vec::new(),
    };
    b.collect()?;
    b.cluster();
    b.snap()?;
    let (mut verts, vert_of) = b.vertices()?;
    let mut pieces = b.split(&mut verts, &vert_of);
    b.eliminate(&mut pieces);
    b.cleanup(&mut pieces, &mut verts);
    let used: Vec<bool> = {
        let mut u = vec![false; verts.len()];
        for p in pieces.iter().filter(|p| p.state == State::Alive) {
            u[p.v0] = true;
            u[p.v1] = true;
        }
        u
    };
    prune(&mut pieces, verts.len());
    b.finish(&verts, &used, &pieces)
}

impl Builder<'_> {
    fn add(&mut self, curve: usize, t: f64, end: bool, rank: u8) -> usize {
        let c = &self.crvs[curve];
        let t = if c.periodic {
            let w = t.rem_euclid(TAU);
            if w >= TAU {
                0.0
            } else {
                w
            }
        } else {
            t.clamp(c.lo, c.hi)
        };
        self.incs.push(Inc {
            curve,
            t,
            p: c.point(t),
            end,
            rank,
        });
        self.uf.add()
    }

    /// [`rank::CROSSING`] if curves `a` and `b` cross transversally at
    /// parameters `ta` and `tb`, else [`rank::CONTACT`].
    fn contact_rank(&self, a: usize, ta: f64, b: usize, tb: f64) -> u8 {
        let da = self.crvs[a].eval(self.crvs[a].clamp(ta), true)[1];
        let db = self.crvs[b].eval(self.crvs[b].clamp(tb), true)[1];
        let sine = da.perp_dot(db).abs() / (da.length() * db.length());
        if sine.is_finite() && sine >= TRANSVERSAL_SINE {
            rank::CROSSING
        } else {
            rank::CONTACT
        }
    }

    /// End points, sharp spline corners, self-contacts and pairwise contacts.
    fn collect(&mut self) -> Result<(), KernelError> {
        let crvs = self.crvs;
        let tol = self.tol;
        let mut ctx = Ctx {
            tol,
            budget: BUDGET,
        };
        for (i, c) in crvs.iter().enumerate() {
            if !c.periodic {
                self.add(i, c.lo, true, rank::EXACT);
                self.add(i, c.hi, true, rank::EXACT);
            }
            for t in c.corners() {
                self.add(i, t, false, rank::EXACT);
            }
            let mono = if matches!(c.geom, Geom::Spline { .. }) {
                c.monotone_pieces(tol)
                    .into_iter()
                    .map(|p| {
                        let (lo, hi) = super::bezier::bounds(&p.bez.hull());
                        (p.t0, p.t1, lo, hi)
                    })
                    .collect()
            } else {
                Vec::new()
            };
            self.mono.push(mono);
            for (s, t) in intersect::self_contacts(c, &mut ctx)? {
                let r = self.contact_rank(i, s, i, t);
                let a = self.add(i, s, false, r);
                let b = self.add(i, t, false, r);
                self.uf.union(a, b);
            }
        }
        let mut order: Vec<usize> = (0..crvs.len()).collect();
        order.sort_by(|&a, &b| crvs[a].bmin.x.total_cmp(&crvs[b].bmin.x).then(a.cmp(&b)));
        for (k, &i) in order.iter().enumerate() {
            for &j in &order[k + 1..] {
                if crvs[j].bmin.x > crvs[i].bmax.x + tol {
                    break;
                }
                if crvs[j].bmin.y > crvs[i].bmax.y + tol || crvs[i].bmin.y > crvs[j].bmax.y + tol {
                    continue;
                }
                let (a, b) = (i.min(j), i.max(j));
                for (ta, tb) in intersect::pair_contacts(&crvs[a], &crvs[b], &mut ctx)? {
                    let r = self.contact_rank(a, ta, b, tb);
                    let x = self.add(a, ta, false, r);
                    let y = self.add(b, tb, false, r);
                    self.uf.union(x, y);
                }
            }
        }
        Ok(())
    }

    /// Unites all incidences closer than the tolerance, using a grid of
    /// cells of half the tolerance (any two points in one cell are within
    /// the tolerance).
    fn cluster(&mut self) {
        let cell = 0.5 * self.tol;
        let key = |p: DVec2| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (k, inc) in self.incs.iter().enumerate() {
            grid.entry(key(inc.p)).or_default().push(k);
        }
        let mut keys: Vec<(i64, i64)> = grid.keys().copied().collect();
        keys.sort_unstable();
        for k in &keys {
            let list = &grid[k];
            for w in list.windows(2) {
                self.uf.union(w[0], w[1]);
            }
        }
        for &(cx, cy) in &keys {
            let a = &grid[&(cx, cy)];
            for dx in -2i64..=2 {
                for dy in -2i64..=2 {
                    if (dx, dy) <= (0, 0) {
                        continue;
                    }
                    let Some(b) = grid.get(&(cx + dx, cy + dy)) else {
                        continue;
                    };
                    if self.uf.find(a[0]) == self.uf.find(b[0]) {
                        continue;
                    }
                    'pairs: for &i in a {
                        for &j in b {
                            if (self.incs[i].p - self.incs[j].p).length() <= self.tol {
                                self.uf.union(i, j);
                                break 'pairs;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Clusters as sorted member lists, ordered by their first member.
    fn groups(&mut self) -> Vec<Vec<usize>> {
        let mut by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for k in 0..self.incs.len() {
            let r = self.uf.find(k);
            by_root.entry(r).or_default().push(k);
        }
        let mut g: Vec<Vec<usize>> = by_root.into_values().collect();
        g.sort_by_key(|m| m[0]);
        g
    }

    fn mean(&self, members: &[usize]) -> DVec2 {
        members.iter().map(|&m| self.incs[m].p).sum::<DVec2>() / members.len() as f64
    }

    /// Mean of the best-ranked members.
    fn precise_mean(&self, members: &[usize]) -> DVec2 {
        let best = members.iter().map(|&m| self.incs[m].rank).min();
        let top: Vec<usize> = members
            .iter()
            .copied()
            .filter(|&m| Some(self.incs[m].rank) == best)
            .collect();
        self.mean(&top)
    }

    /// Parameters of the points of curve `ci` nearest to `p` (one per
    /// monotone piece for splines), with their distances.
    fn near_params(&self, ci: usize, p: DVec2) -> Vec<(f64, f64)> {
        let c = &self.crvs[ci];
        if self.mono[ci].is_empty() {
            return vec![c.closest(p, c.lo, c.hi)];
        }
        let tol = self.tol;
        self.mono[ci]
            .iter()
            .filter(|(_, _, lo, hi)| {
                p.x >= lo.x - tol && p.y >= lo.y - tol && p.x <= hi.x + tol && p.y <= hi.y + tol
            })
            .map(|&(t0, t1, _, _)| c.closest(p, t0, t1))
            .collect()
    }

    /// Routes every curve passing within the tolerance of a vertex through
    /// it, until nothing changes.
    fn snap(&mut self) -> Result<(), KernelError> {
        let tol = self.tol;
        for _ in 0..12 {
            let groups = self.groups();
            let mut new: Vec<(usize, f64, usize)> = Vec::new();
            for g in &groups {
                let pos = self.mean(g);
                // This group's additions so far (a suffix of `new`).
                let first_new = new.len();
                for (ci, c) in self.crvs.iter().enumerate() {
                    if pos.x < c.bmin.x - tol
                        || pos.y < c.bmin.y - tol
                        || pos.x > c.bmax.x + tol
                        || pos.y > c.bmax.y + tol
                    {
                        continue;
                    }
                    for (t, d) in self.near_params(ci, pos) {
                        if d > tol {
                            continue;
                        }
                        let pt = c.point(t);
                        let has = g.iter().any(|&m| {
                            self.incs[m].curve == ci && (self.incs[m].p - pt).length() <= 2.0 * tol
                        }) || new[first_new..].iter().any(|&(c2, t2, _)| {
                            c2 == ci && (c.point(t2) - pt).length() <= 2.0 * tol
                        });
                        if !has {
                            new.push((ci, t, g[0]));
                        }
                    }
                }
            }
            if new.is_empty() {
                return Ok(());
            }
            for (ci, t, m) in new {
                let k = self.add(ci, t, false, rank::SNAPPED);
                self.uf.union(k, m);
            }
            self.cluster();
        }
        Err(KernelError::InvalidInput(
            "curve contacts within the model tolerance do not settle; curves nearly coincide"
                .into(),
        ))
    }

    /// Final vertices and the vertex of each incidence.
    fn vertices(&mut self) -> Result<(Vec<Vtx>, Vec<usize>), KernelError> {
        let groups = self.groups();
        let mut vert_of = vec![0; self.incs.len()];
        let mut verts = Vec::with_capacity(groups.len());
        for (v, g) in groups.into_iter().enumerate() {
            let pos = self.mean(&g);
            let spread = g
                .iter()
                .map(|&m| (self.incs[m].p - pos).length())
                .fold(0.0, f64::max);
            if spread > 4.0 * self.tol {
                return Err(near_degenerate(
                    pos,
                    "too many curve features within the tolerance",
                ));
            }
            for &m in &g {
                vert_of[m] = v;
            }
            verts.push(Vtx {
                pos,
                point: self.precise_mean(&g),
                end: g.iter().any(|&m| self.incs[m].end),
                members: g,
                removed: false,
                params: BTreeMap::new(),
            });
        }
        Ok((verts, vert_of))
    }

    /// True if the piece `[t0, t1]` of curve `ci` gets farther than `r` from `p`.
    fn leaves(&self, ci: usize, t0: f64, t1: f64, p: DVec2, r: f64) -> bool {
        self.crvs[ci]
            .pieces(t0, t1)
            .iter()
            .any(|pc| (0..=16).any(|k| (pc.bez.point(k as f64 / 16.0) - p).length() > r))
    }

    /// Splits every curve at its vertices. Where a curve's incidences at one
    /// vertex follow each other without the curve leaving the vertex, they
    /// are one stop, at the parameter of the best-ranked of them (then the
    /// one nearest the vertex's point), so that consecutive pieces of a
    /// curve meet exactly and pieces of different curves meet where their
    /// most precise contact is.
    fn split(&mut self, verts: &mut Vec<Vtx>, vert_of: &[usize]) -> Vec<Piece> {
        let n = self.crvs.len();
        // (parameter, vertex, rank) of each incidence, by curve.
        let mut by_curve: Vec<Vec<(f64, usize, u8)>> = vec![Vec::new(); n];
        for (k, inc) in self.incs.iter().enumerate() {
            by_curve[inc.curve].push((inc.t, vert_of[k], inc.rank));
        }
        let mut pieces = Vec::new();
        for (ci, list) in by_curve.iter_mut().enumerate() {
            let c = &self.crvs[ci];
            list.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
            list.dedup_by(|b, a| a.0 == b.0 && a.1 == b.1);
            if c.periodic && list.is_empty() {
                // An isolated circle gets one vertex, at parameter 0.
                let k = self.add(ci, 0.0, false, rank::EXACT);
                let v = verts.len();
                verts.push(Vtx {
                    pos: self.incs[k].p,
                    point: self.incs[k].p,
                    members: vec![k],
                    end: false,
                    removed: false,
                    params: BTreeMap::new(),
                });
                list.push((0.0, v, rank::EXACT));
            }
            let m = list.len();
            if m == 0 {
                continue;
            }
            // Stop `i` of the (cyclically extended) list.
            let stop = |i: usize| {
                let (t, v, r) = list[i % m];
                (t + TAU * (i / m) as f64, v, r)
            };
            // True if the curve stays at one vertex from stop `i` to `i + 1`.
            let stays = |i: usize| {
                let (a, b) = (stop(i), stop(i + 1));
                a.1 == b.1 && !self.leaves(ci, a.0, b.0, verts[a.1].pos, 2.0 * self.tol)
            };
            // Runs of stops at one vertex, starting after a place where the
            // curve leaves its vertex (for a closed curve).
            let (first, count) = if c.periodic {
                match (1..=m).find(|&b| !stays(b - 1)) {
                    Some(b) => (b, m),
                    // The curve never leaves its vertex: nothing to split.
                    None => continue,
                }
            } else {
                (0, m)
            };
            // (parameter in the domain, vertex) of each run.
            let mut runs: Vec<(f64, usize)> = Vec::new();
            let mut i = first;
            while i < first + count {
                let mut j = i;
                while j + 1 < first + count && stays(j) {
                    j += 1;
                }
                let v = stop(i).1;
                let at = verts[v].point;
                let mut best = list[i % m];
                for k in i + 1..=j {
                    let s = list[k % m];
                    let closer = (c.point(s.0) - at).length() < (c.point(best.0) - at).length();
                    if s.2 < best.2 || (s.2 == best.2 && closer) {
                        best = s;
                    }
                }
                runs.push((best.0, v));
                i = j + 1;
            }
            if c.periodic {
                runs.push(runs[0]);
            }
            let first_piece = pieces.len();
            for w in runs.windows(2) {
                let ((t0, v0), (t, v1)) = (w[0], w[1]);
                // A closed curve's piece may run on past 2π.
                let t1 = if t > t0 { t } else { t + TAU };
                pieces.push(Piece {
                    curve: ci,
                    t0,
                    t1,
                    v0,
                    v1,
                    state: State::Alive,
                    next: usize::MAX,
                });
            }
            for &(t, v) in &runs {
                let e = verts[v].params.entry(ci).or_insert(t);
                if t < *e {
                    *e = t;
                }
            }
            let last_piece = pieces.len();
            for (k, p) in pieces.iter_mut().enumerate().skip(first_piece) {
                if k + 1 < last_piece {
                    p.next = k + 1;
                } else if c.periodic && k > first_piece {
                    p.next = first_piece;
                }
            }
        }
        pieces
    }

    /// True if every point of `p` is within the tolerance of `q` (end points
    /// within twice the tolerance, since both ends may be off their vertex).
    fn one_sided(&self, p: &Piece, q: &Piece) -> bool {
        let cp = &self.crvs[p.curve];
        let cq = &self.crvs[q.curve];
        let n = match &cp.geom {
            Geom::Line { .. } => 8,
            Geom::Arc { .. } => 16,
            Geom::Spline { spans } => {
                let k = spans.iter().filter(|s| s.t1 > p.t0 && s.t0 < p.t1).count();
                (8 * k.max(2)).min(512)
            }
        };
        let order = std::iter::once(n / 2).chain((0..=n).filter(|&k| k != n / 2));
        for k in order {
            let t = if k == n {
                p.t1
            } else {
                p.t0 + (p.t1 - p.t0) * k as f64 / n as f64
            };
            let lim = if k == 0 || k == n {
                2.0 * self.tol
            } else {
                self.tol
            };
            if cq.closest(cp.point(t), q.t0, q.t1).1 > lim {
                return false;
            }
        }
        true
    }

    /// Marks pieces that duplicate a kept piece between the same vertices.
    fn eliminate(&self, pieces: &mut [Piece]) {
        let mut groups: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        for (k, p) in pieces.iter().enumerate() {
            groups
                .entry((p.v0.min(p.v1), p.v0.max(p.v1)))
                .or_default()
                .push(k);
        }
        for g in groups.values() {
            if g.len() < 2 {
                continue;
            }
            let mut kept: Vec<usize> = Vec::new();
            for &k in g {
                let dup = kept.iter().any(|&q| {
                    self.one_sided(&pieces[k], &pieces[q]) && self.one_sided(&pieces[q], &pieces[k])
                });
                if dup {
                    pieces[k].state = State::Duplicate;
                } else {
                    kept.push(k);
                }
            }
        }
    }

    /// Removes vertices that only split one curve inside an overlap: two
    /// live pieces of the same curve, consecutive in parameter, and nothing
    /// else live there.
    fn cleanup(&self, pieces: &mut [Piece], verts: &mut [Vtx]) {
        loop {
            let mut adj: Vec<Vec<(usize, bool)>> = vec![Vec::new(); verts.len()];
            for (k, p) in pieces.iter().enumerate() {
                if p.state == State::Alive {
                    adj[p.v0].push((k, true));
                    adj[p.v1].push((k, false));
                }
            }
            let mut touched = vec![false; pieces.len()];
            let mut changed = false;
            for (v, vx) in verts.iter_mut().enumerate() {
                if vx.end || vx.removed || adj[v].len() != 2 {
                    continue;
                }
                let ((k1, s1), (k2, s2)) = (adj[v][0], adj[v][1]);
                if k1 == k2 || touched[k1] || touched[k2] {
                    continue;
                }
                let (pe, ps) = match (s1, s2) {
                    (false, true) => (k1, k2),
                    (true, false) => (k2, k1),
                    _ => continue,
                };
                // Consecutive pieces of one curve (possibly across an
                // interval that collapsed into the vertex).
                if pieces[pe].next != ps {
                    continue;
                }
                let c = &self.crvs[pieces[pe].curve];
                if c.periodic {
                    let g = (pieces[ps].t0 - pieces[pe].t1).rem_euclid(TAU);
                    let g = if g > PI { g - TAU } else { g };
                    pieces[pe].t1 += g + (pieces[ps].t1 - pieces[ps].t0);
                } else {
                    pieces[pe].t1 = pieces[ps].t1;
                }
                pieces[pe].v1 = pieces[ps].v1;
                pieces[pe].next = pieces[ps].next;
                pieces[ps].state = State::Merged;
                touched[pe] = true;
                touched[ps] = true;
                vx.removed = true;
                changed = true;
            }
            if !changed {
                break;
            }
        }
    }

    /// Unit direction of travel and signed curvature where a half-edge
    /// leaves its start vertex.
    fn travel(&self, p: &Piece, fwd: bool) -> (DVec2, f64) {
        let c = &self.crvs[p.curve];
        let (t, sgn) = if fwd { (p.t0, 1.0) } else { (p.t1, -1.0) };
        let [pt, d1, d2] = c.eval(t, fwd);
        let d1 = d1 * sgn;
        let chord = c.point(0.5 * (p.t0 + p.t1)) - pt;
        let speed = d1.length();
        let dir = if speed.is_finite() && speed * (p.t1 - p.t0) > 1e-9 * chord.length() {
            d1 / speed
        } else if let Some(d) = d2.try_normalize() {
            d
        } else {
            chord.try_normalize().unwrap_or(DVec2::X)
        };
        let kappa = if speed > 0.0 {
            d1.perp_dot(d2) / (speed * speed * speed)
        } else {
            0.0
        };
        (dir, if kappa.is_finite() { kappa } else { 0.0 })
    }

    /// Lateral offset (left positive) from `v` along `dir` of the point where
    /// the half-edge first reaches the longitudinal distance `x`.
    fn probe(&self, p: &Piece, fwd: bool, v: DVec2, dir: DVec2, x: f64) -> Option<f64> {
        let c = &self.crvs[p.curve];
        let (ts, te) = if fwd { (p.t0, p.t1) } else { (p.t1, p.t0) };
        let f = |t: f64| (c.point(t) - v).dot(dir);
        let lat = |t: f64| dir.perp_dot(c.point(t) - v);
        if f(ts) >= x {
            return Some(lat(ts));
        }
        let sgn = (te - ts).signum();
        let speed = c.eval(ts, fwd)[1].length();
        let mut h = if speed > 0.0 && (x / speed).is_finite() {
            x / speed
        } else {
            (te - ts).abs() * 1e-6
        };
        let mut lo = ts;
        let mut hi = None;
        for _ in 0..200 {
            let mut tt = ts + sgn * h;
            let past = if sgn > 0.0 { tt >= te } else { tt <= te };
            if past {
                tt = te;
            }
            if f(tt) >= x {
                hi = Some(tt);
                break;
            }
            if past {
                return None;
            }
            lo = tt;
            h *= 2.0;
        }
        let mut hi = hi?;
        for _ in 0..100 {
            let m = 0.5 * (lo + hi);
            if m == lo || m == hi {
                break;
            }
            if f(m) >= x {
                hi = m;
            } else {
                lo = m;
            }
        }
        Some(lat(hi))
    }

    /// Orders tangent-tied half-edges (all leaving `v` in nearly the same
    /// direction) counter-clockwise: by their lateral separation at
    /// increasing distances along the common direction (the order in which
    /// they cross a small circle about the vertex), then by signed curvature.
    fn order_tied(
        &self,
        pieces: &[Piece],
        v: DVec2,
        tied: &[Tied],
    ) -> Result<Vec<usize>, KernelError> {
        let dir = tied
            .iter()
            .map(|t| t.dir)
            .sum::<DVec2>()
            .try_normalize()
            .unwrap_or(tied[0].dir);
        let limit = tied
            .iter()
            .map(|t| {
                let p = &pieces[t.he / 2];
                let c = &self.crvs[p.curve];
                let far = if t.he % 2 == 0 { p.t1 } else { p.t0 };
                let mid = c.point(0.5 * (p.t0 + p.t1));
                0.45 * (mid - v).length().max((c.point(far) - v).length())
            })
            .fold(f64::INFINITY, f64::min);
        let scale = v.abs().max_element();
        let mut groups: Vec<Vec<usize>> = vec![(0..tied.len()).collect()];
        let mut x = (32.0 * self.tol).min(limit / 8.0);
        while x <= limit && groups.iter().any(|g| g.len() > 1) {
            let noise = 64.0 * f64::EPSILON * (scale + x);
            let mut next = Vec::new();
            for g in groups {
                if g.len() == 1 {
                    next.push(g);
                    continue;
                }
                let ys: Option<Vec<f64>> = g
                    .iter()
                    .map(|&k| {
                        let he = tied[k].he;
                        self.probe(&pieces[he / 2], he % 2 == 0, v, dir, x)
                    })
                    .collect();
                let Some(ys) = ys else {
                    next.push(g);
                    continue;
                };
                let mut idx: Vec<usize> = (0..g.len()).collect();
                idx.sort_by(|&a, &b| ys[a].total_cmp(&ys[b]).then(g[a].cmp(&g[b])));
                let mut cur = vec![g[idx[0]]];
                for w in idx.windows(2) {
                    if ys[w[1]] - ys[w[0]] > noise {
                        next.push(std::mem::take(&mut cur));
                    }
                    cur.push(g[w[1]]);
                }
                next.push(cur);
            }
            groups = next;
            x *= 2.0;
        }
        let mut out = Vec::with_capacity(tied.len());
        for mut g in groups {
            if g.len() > 1 {
                g.sort_by(|&a, &b| tied[a].kappa.total_cmp(&tied[b].kappa));
                for w in g.windows(2) {
                    let (ka, kb) = (tied[w[0]].kappa, tied[w[1]].kappa);
                    if kb - ka <= 1e-9 * (ka.abs() + kb.abs()) {
                        return Err(near_degenerate(v, "coincident curves at a vertex"));
                    }
                }
            }
            out.extend(g.into_iter().map(|k| tied[k].he));
        }
        Ok(out)
    }

    /// Outgoing half-edges of vertex `v` in counter-clockwise order.
    fn rotation(
        &self,
        pieces: &[Piece],
        v: DVec2,
        hes: &[usize],
    ) -> Result<Vec<usize>, KernelError> {
        if hes.len() <= 2 {
            return Ok(hes.to_vec());
        }
        let mut info: Vec<(Tied, f64)> = hes
            .iter()
            .map(|&he| {
                let (dir, kappa) = self.travel(&pieces[he / 2], he % 2 == 0);
                (Tied { he, dir, kappa }, dir.y.atan2(dir.x))
            })
            .collect();
        info.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.he.cmp(&b.0.he)));
        let mut groups: Vec<Vec<Tied>> = Vec::new();
        for (k, (t, ang)) in info.iter().enumerate() {
            if k == 0 || ang - info[k - 1].1 >= TIE_ANGLE {
                groups.push(Vec::new());
            }
            if let Some(g) = groups.last_mut() {
                g.push(*t);
            }
        }
        if groups.len() > 1 && info[0].1 + TAU - info[info.len() - 1].1 < TIE_ANGLE {
            if let Some(mut last) = groups.pop() {
                last.append(&mut groups[0]);
                groups[0] = last;
            }
        }
        let mut out = Vec::with_capacity(hes.len());
        for g in groups {
            if g.len() == 1 {
                out.push(g[0].he);
            } else {
                out.extend(self.order_tied(pieces, v, &g)?);
            }
        }
        Ok(out)
    }

    /// Faces, regions, checks and output.
    fn finish(
        &self,
        verts: &[Vtx],
        used: &[bool],
        all: &[Piece],
    ) -> Result<SketchArrangement, KernelError> {
        let crvs = self.crvs;
        let tol = self.tol;
        // Live edges, ordered by curve and parameter.
        let mut eidx: Vec<usize> = (0..all.len())
            .filter(|&k| all[k].state == State::Alive)
            .collect();
        eidx.sort_by(|&a, &b| {
            (all[a].curve, all[a].t0)
                .partial_cmp(&(all[b].curve, all[b].t0))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let edges: Vec<Piece> = eidx
            .iter()
            .map(|&k| {
                let mut p = all[k].clone();
                if crvs[p.curve].periodic && p.t0 >= TAU {
                    p.t0 -= TAU;
                    p.t1 -= TAU;
                }
                p
            })
            .collect();
        let ne = edges.len();
        let nv = verts.len();
        let start = |h: usize| {
            if h % 2 == 0 {
                edges[h / 2].v0
            } else {
                edges[h / 2].v1
            }
        };
        let end = |h: usize| {
            if h % 2 == 0 {
                edges[h / 2].v1
            } else {
                edges[h / 2].v0
            }
        };
        let side = |h: usize| Side {
            curve: edges[h / 2].curve,
            t0: edges[h / 2].t0,
            t1: edges[h / 2].t1,
            fwd: h % 2 == 0,
        };
        // Rotation system.
        let mut out_he: Vec<Vec<usize>> = vec![Vec::new(); nv];
        for h in 0..2 * ne {
            out_he[start(h)].push(h);
        }
        let mut rot: Vec<Vec<usize>> = vec![Vec::new(); nv];
        let mut pos_in = vec![0usize; 2 * ne];
        for v in 0..nv {
            if out_he[v].is_empty() {
                continue;
            }
            let r = self.rotation(&edges, verts[v].pos, &out_he[v])?;
            for (i, &h) in r.iter().enumerate() {
                pos_in[h] = i;
            }
            rot[v] = r;
        }
        let next = |h: usize| {
            let r = &rot[end(h)];
            let i = pos_in[h ^ 1];
            r[(i + r.len() - 1) % r.len()]
        };
        // Face cycles.
        let mut cycle_of = vec![usize::MAX; 2 * ne];
        let mut cycles: Vec<Vec<usize>> = Vec::new();
        for h0 in 0..2 * ne {
            if cycle_of[h0] != usize::MAX {
                continue;
            }
            let id = cycles.len();
            let mut cyc = Vec::new();
            let mut h = h0;
            loop {
                if cycle_of[h] != usize::MAX {
                    if h == h0 {
                        break;
                    }
                    return Err(inconsistent("face cycles overlap"));
                }
                cycle_of[h] = id;
                cyc.push(h);
                h = next(h);
                if cyc.len() > 2 * ne {
                    return Err(inconsistent("face cycle does not close"));
                }
            }
            cycles.push(cyc);
        }
        let sides: Vec<Vec<Side>> = cycles
            .iter()
            .map(|c| c.iter().map(|&h| side(h)).collect())
            .collect();
        let areas: Vec<f64> = sides.iter().map(|s| measure::loop_area(crvs, s)).collect();
        // Components and Euler's formula.
        let mut uf = Uf {
            p: (0..nv).collect(),
        };
        for e in &edges {
            uf.union(e.v0, e.v1);
        }
        let mut comp_v: HashMap<usize, (usize, usize, usize, usize)> = HashMap::new();
        for (v, hes) in out_he.iter().enumerate() {
            if !hes.is_empty() {
                comp_v.entry(uf.find(v)).or_default().0 += 1;
            }
        }
        for e in &edges {
            comp_v.entry(uf.find(e.v0)).or_default().1 += 1;
        }
        for (c, s) in cycles.iter().zip(&areas) {
            let k = comp_v.entry(uf.find(start(c[0]))).or_default();
            k.2 += 1;
            if *s < 0.0 {
                k.3 += 1;
            }
        }
        for &(v, e, f, neg) in comp_v.values() {
            if v + f != e + 2 {
                return Err(inconsistent(format!(
                    "Euler's formula fails (V={v}, E={e}, F={f})"
                )));
            }
            if neg != 1 {
                return Err(inconsistent(format!(
                    "{neg} outer boundaries in one component"
                )));
            }
        }
        for (c, &a) in cycles.iter().zip(&areas) {
            if a.abs() <= tol * tol || !a.is_finite() {
                return Err(near_degenerate(verts[start(c[0])].pos, "degenerate face"));
            }
        }
        // Regions: positive cycles; holes: negative cycles inside them.
        let mut region_of: Vec<Option<usize>> = vec![None; cycles.len()];
        let mut outer_cycles = Vec::new();
        for (ci, &a) in areas.iter().enumerate() {
            if a > 0.0 {
                region_of[ci] = Some(outer_cycles.len());
                outer_cycles.push(ci);
            }
        }
        let outer_bounds: Vec<(DVec2, DVec2)> = outer_cycles
            .iter()
            .map(|&c| measure::loop_bounds(crvs, &sides[c]))
            .collect();
        let mut holes_of: Vec<Vec<usize>> = vec![Vec::new(); outer_cycles.len()];
        let mut face_of: Vec<Option<usize>> = region_of.clone();
        for (ci, &a) in areas.iter().enumerate() {
            if a > 0.0 {
                continue;
            }
            let comp = uf.find(start(cycles[ci][0]));
            let owner = self.hole_owner(
                &sides,
                ci,
                &outer_cycles,
                &outer_bounds,
                &areas,
                |c: usize| uf_find_const(&uf, start(cycles[c][0])) == comp,
            )?;
            if let Some(r) = owner {
                holes_of[r].push(ci);
            }
            face_of[ci] = owner;
        }
        // Regions with interior points, checked against each other.
        let mut regions = Vec::with_capacity(outer_cycles.len());
        let mut shapes = Vec::with_capacity(outer_cycles.len());
        for (r, &oc) in outer_cycles.iter().enumerate() {
            let area = areas[oc] + holes_of[r].iter().map(|&h| areas[h]).sum::<f64>();
            let anchor = verts[start(cycles[oc][0])].pos;
            if area.is_nan() || area <= tol * tol {
                return Err(near_degenerate(anchor, "region with no area"));
            }
            let mut loops = vec![sides[oc].clone()];
            loops.extend(holes_of[r].iter().map(|&h| sides[h].clone()));
            let shape = RegionShape::new(crvs, loops);
            let Some((p, clearance)) = shape.interior_point(tol) else {
                return Err(near_degenerate(anchor, "region without an interior point"));
            };
            if clearance < tol {
                return Err(near_degenerate(p, "region thinner than the tolerance"));
            }
            let to_loop = |c: usize| RegionLoop {
                edges: cycles[c].iter().map(|&h| (h / 2, h % 2 == 0)).collect(),
            };
            regions.push(Region {
                outer: to_loop(oc),
                holes: holes_of[r].iter().map(|&h| to_loop(h)).collect(),
                area,
                interior_point: p,
            });
            shapes.push(shape);
        }
        for (r, reg) in regions.iter().enumerate() {
            for (s, shape) in shapes.iter().enumerate() {
                if s != r && shape.winding(reg.interior_point) != Some(0) {
                    return Err(inconsistent(format!("regions {r} and {s} overlap")));
                }
            }
        }
        // Output with compact vertex indices.
        let mut vmap = vec![usize::MAX; nv];
        let mut vertices = Vec::new();
        for (v, vx) in verts.iter().enumerate() {
            if !used[v] || vx.removed {
                continue;
            }
            vmap[v] = vertices.len();
            let mut cp: BTreeMap<usize, f64> = BTreeMap::new();
            for &m in &vx.members {
                let inc = &self.incs[m];
                let e = cp.entry(inc.curve).or_insert(inc.t);
                if inc.t < *e {
                    *e = inc.t;
                }
            }
            // Where the curve's pieces end, when it has any here.
            for (ci, t) in &vx.params {
                cp.insert(*ci, *t);
            }
            vertices.push(ArrangementVertex {
                point: vx.point,
                curves: cp.keys().copied().collect(),
                params: cp.values().copied().collect(),
            });
        }
        let mut out_edges = Vec::with_capacity(ne);
        for (k, e) in edges.iter().enumerate() {
            if vmap[e.v0] == usize::MAX || vmap[e.v1] == usize::MAX {
                return Err(inconsistent("edge on a removed vertex"));
            }
            out_edges.push(ArrangementEdge {
                curve: e.curve,
                t0: e.t0,
                t1: e.t1,
                start: vmap[e.v0],
                end: vmap[e.v1],
                left: face_of[cycle_of[2 * k]],
                right: face_of[cycle_of[2 * k + 1]],
            });
        }
        Ok(SketchArrangement {
            vertices,
            edges: out_edges,
            regions,
        })
    }

    /// The region whose outer loop most tightly contains the outer boundary
    /// `ci` of a connected component (`None` if it is unbounded).
    fn hole_owner(
        &self,
        sides: &[Vec<Side>],
        ci: usize,
        outer_cycles: &[usize],
        outer_bounds: &[(DVec2, DVec2)],
        areas: &[f64],
        same_component: impl Fn(usize) -> bool,
    ) -> Result<Option<usize>, KernelError> {
        let crvs = self.crvs;
        let mut order: Vec<&Side> = sides[ci].iter().collect();
        order.sort_by(|a, b| {
            let la = (a.end(crvs) - a.start(crvs)).length();
            let lb = (b.end(crvs) - b.start(crvs)).length();
            lb.total_cmp(&la)
        });
        'points: for s in order {
            let p = crvs[s.curve].point(0.5 * (s.t0 + s.t1));
            let mut best: Option<usize> = None;
            for (r, &oc) in outer_cycles.iter().enumerate() {
                if same_component(oc) {
                    continue;
                }
                let (lo, hi) = outer_bounds[r];
                if p.x < lo.x || p.y < lo.y || p.x > hi.x || p.y > hi.y {
                    continue;
                }
                match measure::winding(crvs, &sides[oc], p) {
                    None => continue 'points,
                    Some(1) => {
                        if best.is_none_or(|b| areas[oc] < areas[outer_cycles[b]]) {
                            best = Some(r);
                        }
                    }
                    Some(0) => {}
                    Some(_) => return Err(inconsistent("a face boundary winds more than once")),
                }
            }
            return Ok(best);
        }
        Err(inconsistent("cannot place a component inside a region"))
    }
}

/// `find` without path compression, for use inside closures.
fn uf_find_const(uf: &Uf, mut x: usize) -> usize {
    while uf.p[x] != x {
        x = uf.p[x];
    }
    x
}

/// Removes dangling pieces: repeatedly, live pieces with an end at a vertex
/// of degree one.
fn prune(pieces: &mut [Piece], nv: usize) {
    let mut deg = vec![0usize; nv];
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (k, p) in pieces.iter().enumerate() {
        if p.state == State::Alive {
            deg[p.v0] += 1;
            deg[p.v1] += 1;
            adj[p.v0].push(k);
            adj[p.v1].push(k);
        }
    }
    let mut stack: Vec<usize> = (0..nv).filter(|&v| deg[v] == 1).collect();
    while let Some(v) = stack.pop() {
        if deg[v] != 1 {
            continue;
        }
        let Some(&k) = adj[v].iter().find(|&&k| pieces[k].state == State::Alive) else {
            continue;
        };
        pieces[k].state = State::Pruned;
        for w in [pieces[k].v0, pieces[k].v1] {
            deg[w] -= 1;
            if deg[w] == 1 {
                stack.push(w);
            }
        }
    }
}
