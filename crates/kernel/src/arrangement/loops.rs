//! Exact checks and measures for the closed loops of a sweep profile
//! (design note 0001, "Intersection check"): the profile builder chains
//! curves into loops by their end points; here the loops are checked for
//! contacts other than at the vertices they share, with the arrangement's
//! contact finder, and measured with its exact area and winding numbers.

use glam::DVec2;

use super::curve::Crv;
use super::intersect::{self, Ctx};
use super::measure::{self, Side};
use super::topology::BUDGET;
use crate::bvh::Bvh;
use crate::{KernelError, ProfileCurve};

/// What profile normalization needs to know about its loops.
#[derive(Clone, Debug)]
pub(crate) struct LoopFacts {
    /// Signed area of each loop (positive when counter-clockwise).
    pub areas: Vec<f64>,
    /// Per loop, the loops it lies inside, in increasing order.
    pub containers: Vec<Vec<usize>>,
}

/// The profile's curves, their vertices and loops.
pub(crate) struct ProfileLoops<'a> {
    /// The curves.
    pub curves: &'a [ProfileCurve],
    /// Each loop as its curves in order, each with its direction (`true`
    /// along the curve's own direction).
    pub loops: &'a [Vec<(usize, bool)>],
    /// The vertices at the start and the end of each curve (one vertex for
    /// a closed curve).
    pub ends: &'a [(usize, usize)],
    /// Vertex positions.
    pub verts: &'a [DVec2],
    /// The distance within which end points are one vertex.
    pub snap: f64,
    /// The linear tolerance.
    pub tol: f64,
}

/// Samples per curve stretch when checking that two curves leaving a
/// vertex stay together.
const RUN_SAMPLES: usize = 64;

impl ProfileLoops<'_> {
    /// Checks that the loops meet only at their vertices and measures them.
    ///
    /// A contact between two curves (or between the two ends of one closed
    /// curve) is allowed only at a vertex where both end: both stretches
    /// from that vertex to the contact stay within `snap` of it, or the two
    /// stretches stay within twice the tolerance of each other all along
    /// (curves leaving a vertex tangentially). Any other contact is
    /// [`KernelError::InvalidInput`] naming where.
    pub fn check(&self) -> Result<LoopFacts, KernelError> {
        let crvs: Vec<Crv> = self
            .curves
            .iter()
            .enumerate()
            .map(|(i, c)| Crv::new(i, c, self.tol))
            .collect::<Result<_, _>>()?;
        let mut ctx = Ctx {
            tol: self.tol,
            budget: BUDGET,
        };
        let touching = |p: DVec2| {
            KernelError::InvalidInput(format!(
                "profile curves intersect or touch each other near ({:.9}, {:.9})",
                p.x, p.y
            ))
        };
        let tol = self.tol;
        // Candidate pairs come from a tree of the curves' boxes, padded by
        // twice the tolerance so rounding cannot lose a pair the exact test
        // below keeps; pairs are taken in the order of a full scan.
        let pad = DVec2::splat(2.0 * tol);
        let tree = Bvh::new(
            crvs.iter()
                .map(|c| ((c.bmin - pad).extend(0.0), (c.bmax + pad).extend(0.0)))
                .collect(),
        );
        let mut near = Vec::new();
        for i in 0..crvs.len() {
            for (s, t) in intersect::self_contacts(&crvs[i], &mut ctx)? {
                if !self.allowed(&crvs, i, s, i, t) {
                    return Err(touching(crvs[i].point(s)));
                }
            }
            tree.overlapping(
                crvs[i].bmin.extend(0.0),
                crvs[i].bmax.extend(0.0),
                &mut near,
            );
            for &j in near.iter().filter(|&&j| j > i) {
                let (a, b) = (&crvs[i], &crvs[j]);
                if a.bmin.x > b.bmax.x + tol
                    || b.bmin.x > a.bmax.x + tol
                    || a.bmin.y > b.bmax.y + tol
                    || b.bmin.y > a.bmax.y + tol
                {
                    continue;
                }
                for (s, t) in intersect::pair_contacts(a, b, &mut ctx)? {
                    if !self.allowed(&crvs, i, s, j, t) {
                        return Err(touching(a.point(s)));
                    }
                }
            }
        }
        let sides: Vec<Vec<Side>> = self
            .loops
            .iter()
            .map(|lp| {
                lp.iter()
                    .map(|&(c, fwd)| Side {
                        curve: c,
                        t0: crvs[c].lo,
                        t1: crvs[c].hi,
                        fwd,
                    })
                    .collect()
            })
            .collect();
        let areas: Vec<f64> = sides.iter().map(|s| measure::loop_area(&crvs, s)).collect();
        // Loops meet nowhere, so every point of one loop is inside another
        // loop or every point is outside it: one probe decides.
        let probes: Vec<DVec2> = self
            .loops
            .iter()
            .map(|lp| {
                let c = &crvs[lp[0].0];
                c.point(0.5 * (c.lo + c.hi))
            })
            .collect();
        // Only loops whose boxes hold the probe can wind around it, or have
        // it on their boundary (a curve lies in the hull of its Bézier
        // control points, which the box bounds).
        let tree = Bvh::new(
            sides
                .iter()
                .map(|s| {
                    let (lo, hi) = measure::loop_bounds(&crvs, s);
                    ((lo - pad).extend(0.0), (hi + pad).extend(0.0))
                })
                .collect(),
        );
        let mut containers = vec![Vec::new(); sides.len()];
        for (i, &p) in probes.iter().enumerate() {
            tree.overlapping(p.extend(0.0), p.extend(0.0), &mut near);
            for &j in near.iter().filter(|&&j| j != i) {
                let w = measure::winding(&crvs, &sides[j], p).ok_or_else(|| touching(p))?;
                if w != 0 {
                    containers[i].push(j);
                }
            }
        }
        Ok(LoopFacts { areas, containers })
    }

    /// True if the contact of curve `i` at `s` with curve `j` at `t` is
    /// part of a vertex at which both curves end (see [`Self::check`]).
    fn allowed(&self, crvs: &[Crv], i: usize, s: f64, j: usize, t: f64) -> bool {
        let (a, b) = (&crvs[i], &crvs[j]);
        if a.periodic || b.periodic {
            return false;
        }
        let ends = |c: usize, crv: &Crv| [(self.ends[c].0, crv.lo), (self.ends[c].1, crv.hi)];
        for (vi, ei) in ends(i, a) {
            for (vj, ej) in ends(j, b) {
                if vi != vj || (i == j && ei == ej) {
                    continue;
                }
                let v = self.verts[vi];
                let at_vertex =
                    stays_near(a, ei, s, v, self.snap) && stays_near(b, ej, t, v, self.snap);
                if at_vertex || self.run(a, ei, s, b, ej, t) {
                    return true;
                }
            }
        }
        false
    }

    /// True if the stretches of `a` from `ea` to `s` and of `b` from `eb`
    /// to `t` each stay within twice the tolerance of the other (sampled).
    fn run(&self, a: &Crv, ea: f64, s: f64, b: &Crv, eb: f64, t: f64) -> bool {
        let r = 2.0 * self.tol;
        let near = |x: &Crv, x0: f64, x1: f64, y: &Crv, y0: f64, y1: f64| {
            let (lo, hi) = (y0.min(y1), y0.max(y1));
            (0..=RUN_SAMPLES).all(|k| {
                let u = x0 + (x1 - x0) * k as f64 / RUN_SAMPLES as f64;
                y.closest(x.point(u), lo, hi).1 <= r
            })
        };
        near(a, ea, s, b, eb, t) && near(b, eb, t, a, ea, s)
    }
}

/// True if the stretch of `c` between parameters `t0` and `t1` lies within
/// `r` of `p`: every control point of its Bézier pieces does (the curve lies
/// in their convex hull), so this never accepts a stretch that leaves.
fn stays_near(c: &Crv, t0: f64, t1: f64, p: DVec2, r: f64) -> bool {
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    if hi <= lo {
        return (c.point(lo) - p).length() <= r;
    }
    c.pieces(lo, hi)
        .iter()
        .all(|pc| pc.bez.hull().iter().all(|q| (*q - p).length() <= r))
}
