//! The pieces of a face cut along segments (design note 0012 §4.2).
//!
//! The segments join given points (clusters), so nothing is intersected
//! here. Two segments crossing, or a point within the tolerance of a
//! segment it does not end, means the bodies come closer than the
//! tolerance somewhere the earlier stages found no contact (seen in the
//! face's plane rather than in space), and the boolean is refused there.
//! The segments form a planar graph; walking it with the face on the left
//! of every half-edge gives its cycles. Counter-clockwise cycles bound
//! regions, clockwise ones are the outsides of connected components and
//! become holes of the smallest region of another component around them.
//! A region is inside the face when it runs along the face's own edges in
//! their direction (material on the left); regions bounded only by
//! segments imprinted inside the face are inside it too.

use std::collections::HashMap;

use glam::DVec2;

use super::{FaceSeg, Key};

/// Why a face could not be split, at a point (a cluster).
#[derive(Debug)]
pub(super) enum Problem {
    /// A point within the tolerance of a segment it does not end.
    Close(usize),
    /// Two segments crossing.
    Cross(usize),
    /// The earlier stages disagree: a kernel error.
    Inconsistent(&'static str, usize),
}

impl Problem {
    /// What went wrong, and where (a cluster).
    pub fn describe(&self) -> (&'static str, usize) {
        match *self {
            Problem::Close(c) => (
                "a point within the tolerance of a segment it does not end",
                c,
            ),
            Problem::Cross(c) => ("two segments crossing", c),
            Problem::Inconsistent(m, c) => (m, c),
        }
    }
}

/// The regions inside the face, each as cycles of points (outer first,
/// counter-clockwise; holes clockwise).
pub(super) fn regions(
    segs: &[FaceSeg],
    pos: &dyn Fn(usize) -> DVec2,
    eps: f64,
) -> Result<Vec<Vec<Vec<usize>>>, Problem> {
    // Unique segments, in order of first appearance, with their sides.
    let mut keys: Vec<Key> = Vec::new();
    let mut side: HashMap<Key, Option<bool>> = HashMap::new();
    for s in segs {
        match side.get_mut(&s.key) {
            None => {
                keys.push(s.key);
                side.insert(s.key, s.side);
            }
            Some(old) => match (*old, s.side) {
                (_, None) => {}
                (None, Some(x)) => *old = Some(x),
                (Some(x), Some(y)) if x == y => {}
                _ => {
                    return Err(Problem::Inconsistent(
                        "an edge with material on both sides",
                        s.key.0,
                    ))
                }
            },
        }
    }
    if keys.is_empty() {
        return Ok(Vec::new());
    }
    // Nodes.
    let mut node_of: HashMap<usize, usize> = HashMap::new();
    let mut nodes: Vec<usize> = Vec::new();
    for k in &keys {
        for c in [k.0, k.1] {
            node_of.entry(c).or_insert_with(|| {
                nodes.push(c);
                nodes.len() - 1
            });
        }
    }
    let p: Vec<DVec2> = nodes.iter().map(|&c| pos(c)).collect();
    let ends: Vec<(usize, usize)> = keys
        .iter()
        .map(|k| (node_of[&k.0], node_of[&k.1]))
        .collect();

    // Consistency: no crossings, no point on a segment it does not end.
    for (i, &(a, b)) in ends.iter().enumerate() {
        for &(c, d) in ends.iter().skip(i + 1) {
            if a == c || a == d || b == c || b == d {
                continue;
            }
            let o1 = orient(p[a], p[b], p[c]);
            let o2 = orient(p[a], p[b], p[d]);
            let o3 = orient(p[c], p[d], p[a]);
            let o4 = orient(p[c], p[d], p[b]);
            if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
                return Err(Problem::Cross(nodes[a]));
            }
        }
        for (n, q) in p.iter().enumerate() {
            if n != a && n != b && point_segment(*q, p[a], p[b]) <= eps {
                return Err(Problem::Close(nodes[n]));
            }
        }
    }

    // Half-edges: 2s runs a → b, 2s + 1 runs b → a.
    let nh = 2 * ends.len();
    let from = |h: usize| {
        if h % 2 == 0 {
            ends[h / 2].0
        } else {
            ends[h / 2].1
        }
    };
    let to = |h: usize| {
        if h % 2 == 0 {
            ends[h / 2].1
        } else {
            ends[h / 2].0
        }
    };
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for h in 0..nh {
        out[from(h)].push(h);
    }
    let angle = |h: usize| {
        let d = p[to(h)] - p[from(h)];
        d.y.atan2(d.x)
    };
    let mut rank = vec![0usize; nh];
    for (n, list) in out.iter_mut().enumerate() {
        list.sort_by(|&x, &y| angle(x).total_cmp(&angle(y)));
        for w in list.windows(2) {
            if angle(w[0]) == angle(w[1]) {
                return Err(Problem::Close(nodes[n]));
            }
        }
        for (r, &h) in list.iter().enumerate() {
            rank[h] = r;
        }
    }
    let next = |h: usize| {
        let t = h ^ 1;
        let v = from(t);
        let list = &out[v];
        list[(rank[t] + list.len() - 1) % list.len()]
    };

    // Cycles, with antennas (segments with the same cycle on both sides)
    // removed.
    let mut seen = vec![false; nh];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for h0 in 0..nh {
        if seen[h0] {
            continue;
        }
        let mut cyc = Vec::new();
        let mut h = h0;
        while !seen[h] {
            seen[h] = true;
            cyc.push(h);
            h = next(h);
        }
        if h != h0 {
            return Err(Problem::Inconsistent(
                "a half-edge walk did not close",
                nodes[from(h0)],
            ));
        }
        cycles.push(prune(cyc));
    }
    let area = |c: &[usize]| -> f64 {
        c.iter()
            .map(|&h| p[from(h)].perp_dot(p[to(h)]))
            .sum::<f64>()
            * 0.5
    };

    // Components, by union–find over nodes.
    let mut comp: Vec<usize> = (0..nodes.len()).collect();
    fn find(c: &mut [usize], mut i: usize) -> usize {
        while c[i] != i {
            c[i] = c[c[i]];
            i = c[i];
        }
        i
    }
    for &(a, b) in &ends {
        let (ra, rb) = (find(&mut comp, a), find(&mut comp, b));
        if ra != rb {
            comp[ra.max(rb)] = ra.min(rb);
        }
    }

    // Which side of the face each cycle is on, from its boundary half-edges.
    let inside_of = |c: &[usize]| -> Result<Option<bool>, Problem> {
        let mut verdict = None;
        for &h in c {
            if let Some(s) = side[&keys[h / 2]] {
                let along = (h % 2 == 0) == s;
                match verdict {
                    None => verdict = Some(along),
                    Some(v) if v != along => {
                        return Err(Problem::Inconsistent(
                            "a region both inside and outside the face",
                            nodes[from(h)],
                        ))
                    }
                    _ => {}
                }
            }
        }
        Ok(verdict)
    };

    struct Region {
        cycle: Vec<usize>,
        area: f64,
        comp: usize,
        inside: bool,
        holes: Vec<Vec<usize>>,
    }
    let mut regions: Vec<Region> = Vec::new();
    let mut outsides: Vec<(Vec<usize>, f64)> = Vec::new();
    for c in cycles {
        if c.is_empty() {
            continue;
        }
        let a = area(&c);
        if a > 0.0 {
            let inside = inside_of(&c)?.unwrap_or(true);
            let k = find(&mut comp, from(c[0]));
            regions.push(Region {
                cycle: c,
                area: a,
                comp: k,
                inside,
                holes: Vec::new(),
            });
        } else if a < 0.0 {
            outsides.push((c, a));
        } else {
            return Err(Problem::Inconsistent(
                "a cycle of zero area",
                nodes[from(c[0])],
            ));
        }
    }
    for (c, _) in outsides {
        let k = find(&mut comp, from(c[0]));
        let q = p[from(c[0])];
        let mut best: Option<usize> = None;
        for (i, r) in regions.iter().enumerate() {
            if r.comp == k || !contains(&r.cycle, q, &p, &from) {
                continue;
            }
            if best.is_none_or(|b| r.area < regions[b].area) {
                best = Some(i);
            }
        }
        let verdict = inside_of(&c)?;
        match best {
            Some(i) => {
                if !regions[i].inside || verdict == Some(false) {
                    return Err(Problem::Inconsistent(
                        "a hole outside the face",
                        nodes[from(c[0])],
                    ));
                }
                regions[i].holes.push(c);
            }
            None => {
                if verdict != Some(false) {
                    return Err(Problem::Inconsistent(
                        "an outer boundary inside the face",
                        nodes[from(c[0])],
                    ));
                }
            }
        }
    }
    Ok(regions
        .into_iter()
        .filter(|r| r.inside)
        .map(|r| {
            std::iter::once(&r.cycle)
                .chain(&r.holes)
                .map(|c| c.iter().map(|&h| nodes[from(h)]).collect())
                .collect()
        })
        .collect())
}

/// Removes antennas: a half-edge followed by its twin.
fn prune(mut c: Vec<usize>) -> Vec<usize> {
    loop {
        let n = c.len();
        if n < 2 {
            return c;
        }
        let Some(i) = (0..n).find(|&i| c[(i + 1) % n] == c[i] ^ 1) else {
            return c;
        };
        let j = (i + 1) % n;
        if j == 0 {
            c.remove(i);
            c.remove(0);
        } else {
            c.remove(j);
            c.remove(i);
        }
    }
}

fn orient(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    )
}

fn point_segment(q: DVec2, a: DVec2, b: DVec2) -> f64 {
    let d = b - a;
    let len2 = d.length_squared();
    if len2 == 0.0 {
        return (q - a).length();
    }
    let t = ((q - a).dot(d) / len2).clamp(0.0, 1.0);
    (q - (a + d * t)).length()
}

/// True if `q` is inside the polygon of half-edge cycle `c` (`q` is not on
/// it: the caller's points keep the tolerance from every segment).
fn contains(c: &[usize], q: DVec2, p: &[DVec2], from: &dyn Fn(usize) -> usize) -> bool {
    let mut inside = false;
    for (i, &h) in c.iter().enumerate() {
        let a = p[from(h)];
        let b = p[from(c[(i + 1) % c.len()])];
        if (a.y > q.y) != (b.y > q.y) {
            let x = a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if x > q.x {
                inside = !inside;
            }
        }
    }
    inside
}
