//! Independent oracles for the arrangement property tests.
//!
//! * [`Shapes`]: simple closed shapes (counter-clockwise star polygons and
//!   circles) in general position. The boundary of their union is found by
//!   splitting every boundary element at its crossings with the other shapes
//!   and keeping the pieces outside all other shapes; chaining those pieces
//!   gives the union's boundary loops. The bounded faces of the arrangement
//!   fill exactly the loops that no other loop encloses, so their total area
//!   is the sum of those loops' areas (Green's theorem, exact for lines and
//!   arcs).
//! * [`CellGrid`]: axis-aligned segments on an integer grid. Every face of
//!   their arrangement is a union of unit cells, so flooding cells across
//!   cell sides not covered by a segment gives the faces exactly.

use std::collections::{HashSet, VecDeque};
use std::f64::consts::TAU;

use openshape_kernel::{DVec2, ProfileCurve};

use super::{circle, line};

/// A simple closed shape.
#[derive(Clone, Debug)]
pub enum Shape {
    /// Simple polygon, in either direction ([`Shapes::new`] runs every
    /// boundary counter-clockwise).
    Poly(Vec<DVec2>),
    /// Circle.
    Circle(DVec2, f64),
}

#[derive(Clone, Copy, Debug)]
enum Elem {
    Seg(DVec2, DVec2),
    Circ(DVec2, f64),
}

impl Elem {
    fn point(&self, t: f64) -> DVec2 {
        match *self {
            Elem::Seg(a, b) => a + (b - a) * t,
            Elem::Circ(c, r) => c + r * DVec2::new(t.cos(), t.sin()),
        }
    }

    /// `½ ∫ x dy − y dx` from `t0` to `t1`.
    fn green(&self, t0: f64, t1: f64) -> f64 {
        match *self {
            Elem::Seg(..) => 0.5 * self.point(t0).perp_dot(self.point(t1)),
            Elem::Circ(c, r) => {
                0.5 * (c.perp_dot(self.point(t1) - self.point(t0)) + r * r * (t1 - t0))
            }
        }
    }

    /// Parameters on `self` where it meets `o`.
    fn hits(&self, o: &Elem) -> Vec<f64> {
        let angle = |c: DVec2, p: DVec2| {
            let a = (p - c).y.atan2((p - c).x);
            a.rem_euclid(TAU)
        };
        match (*self, *o) {
            (Elem::Seg(a, b), Elem::Seg(c, d)) => {
                let (da, dc) = (b - a, d - c);
                let den = da.perp_dot(dc);
                if den == 0.0 {
                    return vec![];
                }
                let s = (c - a).perp_dot(dc) / den;
                let t = (c - a).perp_dot(da) / den;
                if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
                    vec![s]
                } else {
                    vec![]
                }
            }
            (Elem::Seg(a, b), Elem::Circ(c, r)) => seg_circle(a, b, c, r),
            (Elem::Circ(c, r), Elem::Seg(a, b)) => seg_circle(a, b, c, r)
                .into_iter()
                .map(|s| angle(c, a + (b - a) * s))
                .collect(),
            (Elem::Circ(c1, r1), Elem::Circ(c2, r2)) => {
                let d = (c2 - c1).length();
                if d == 0.0 || d >= r1 + r2 || d <= (r1 - r2).abs() {
                    return vec![];
                }
                let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
                let h = (r1 * r1 - a * a).max(0.0).sqrt();
                let u = (c2 - c1) / d;
                let m = c1 + u * a;
                vec![angle(c1, m + u.perp() * h), angle(c1, m - u.perp() * h)]
            }
        }
    }
}

fn seg_circle(a: DVec2, b: DVec2, c: DVec2, r: f64) -> Vec<f64> {
    let d = b - a;
    let f = a - c;
    let qa = d.dot(d);
    let qb = 2.0 * f.dot(d);
    let qc = f.dot(f) - r * r;
    let disc = qb * qb - 4.0 * qa * qc;
    if disc <= 0.0 {
        return vec![];
    }
    let sq = disc.sqrt();
    [(-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa)]
        .into_iter()
        .filter(|s| (0.0..=1.0).contains(s))
        .collect()
}

/// Crossing-number point-in-polygon test.
pub fn in_polygon(p: DVec2, poly: &[DVec2]) -> bool {
    let mut inside = false;
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if p.x < x {
                inside = !inside;
            }
        }
    }
    inside
}

/// Distance from `p` to a closed polyline.
pub fn polyline_distance(p: DVec2, poly: &[DVec2]) -> f64 {
    let mut d = f64::INFINITY;
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        let ab = b - a;
        let s = if ab.length_squared() > 0.0 {
            ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0)
        } else {
            0.0
        };
        d = d.min((a + ab * s - p).length());
    }
    d
}

impl Shape {
    /// The shape's boundary as sketch curves.
    pub fn curves(&self) -> Vec<ProfileCurve> {
        match self {
            Shape::Poly(p) => (0..p.len())
                .map(|i| line(p[i], p[(i + 1) % p.len()]))
                .collect(),
            Shape::Circle(c, r) => vec![circle(*c, *r)],
        }
    }

    fn elems(&self) -> Vec<Elem> {
        match self {
            Shape::Poly(p) => (0..p.len())
                .map(|i| Elem::Seg(p[i], p[(i + 1) % p.len()]))
                .collect(),
            Shape::Circle(c, r) => vec![Elem::Circ(*c, *r)],
        }
    }

    /// True if `p` is inside.
    pub fn contains(&self, p: DVec2) -> bool {
        match self {
            Shape::Poly(poly) => in_polygon(p, poly),
            Shape::Circle(c, r) => (p - *c).length() < *r,
        }
    }

    /// The same shape with its boundary counter-clockwise.
    fn counter_clockwise(self) -> Shape {
        match self {
            Shape::Poly(mut p) => {
                let twice_area: f64 = (0..p.len())
                    .map(|i| p[i].perp_dot(p[(i + 1) % p.len()]))
                    .sum();
                if twice_area < 0.0 {
                    p.reverse();
                }
                Shape::Poly(p)
            }
            c => c,
        }
    }

    /// The shape scaled about the origin.
    pub fn scaled(&self, s: f64) -> Shape {
        match self {
            Shape::Poly(p) => Shape::Poly(p.iter().map(|q| *q * s).collect()),
            Shape::Circle(c, r) => Shape::Circle(*c * s, r * s),
        }
    }
}

/// One boundary loop of the union of shapes.
pub struct UnionLoop {
    /// Signed area (positive: outer boundary of a union component).
    pub area: f64,
    /// Fine polygonal approximation.
    pub poly: Vec<DVec2>,
}

/// Shapes in general position and the loops of their union's boundary.
pub struct Shapes {
    /// The shapes.
    pub shapes: Vec<Shape>,
    /// Union boundary loops.
    pub loops: Vec<UnionLoop>,
    /// Loops enclosed by no other loop: they bound the filled union.
    pub top: Vec<usize>,
}

impl Shapes {
    /// Computes the union boundary of `shapes`.
    pub fn new(shapes: Vec<Shape>) -> Shapes {
        // The union is traced along counter-clockwise boundaries.
        let shapes: Vec<Shape> = shapes.into_iter().map(Shape::counter_clockwise).collect();
        struct Piece {
            e: Elem,
            t0: f64,
            t1: f64,
            a: DVec2,
            b: DVec2,
        }
        let mut pieces: Vec<Piece> = Vec::new();
        for (i, s) in shapes.iter().enumerate() {
            for e in s.elems() {
                let mut ts: Vec<f64> = Vec::new();
                for (j, o) in shapes.iter().enumerate() {
                    if j != i {
                        for f in o.elems() {
                            ts.extend(e.hits(&f));
                        }
                    }
                }
                ts.sort_by(f64::total_cmp);
                match e {
                    Elem::Seg(..) => {
                        ts.insert(0, 0.0);
                        ts.push(1.0);
                    }
                    Elem::Circ(..) => {
                        if ts.is_empty() {
                            ts = vec![0.0, TAU];
                        } else {
                            ts.push(ts[0] + TAU);
                        }
                    }
                }
                for w in ts.windows(2) {
                    if w[1] <= w[0] {
                        continue;
                    }
                    let m = e.point(0.5 * (w[0] + w[1]));
                    if shapes
                        .iter()
                        .enumerate()
                        .any(|(j, o)| j != i && o.contains(m))
                    {
                        continue;
                    }
                    pieces.push(Piece {
                        e,
                        t0: w[0],
                        t1: w[1],
                        a: e.point(w[0]),
                        b: e.point(w[1]),
                    });
                }
            }
        }
        let mut used = vec![false; pieces.len()];
        let mut loops = Vec::new();
        for s in 0..pieces.len() {
            if used[s] {
                continue;
            }
            let mut area = 0.0;
            let mut poly = Vec::new();
            let mut k = s;
            loop {
                used[k] = true;
                let p = &pieces[k];
                area += p.e.green(p.t0, p.t1);
                let n = match p.e {
                    Elem::Seg(..) => 1,
                    Elem::Circ(..) => (((p.t1 - p.t0) / TAU * 4096.0).ceil() as usize).max(8),
                };
                for j in 0..n {
                    poly.push(p.e.point(p.t0 + (p.t1 - p.t0) * j as f64 / n as f64));
                }
                // Generic position: exactly one union piece starts here.
                let end = p.b;
                let next = (0..pieces.len())
                    .filter(|&j| !used[j] || j == s)
                    .min_by(|&x, &y| {
                        (pieces[x].a - end)
                            .length()
                            .total_cmp(&(pieces[y].a - end).length())
                    })
                    .expect("union boundary is closed");
                assert!(
                    (pieces[next].a - end).length() < 1e-9,
                    "union boundary does not chain"
                );
                if next == s {
                    break;
                }
                k = next;
            }
            loops.push(UnionLoop { area, poly });
        }
        let top = (0..loops.len())
            .filter(|&i| {
                let p = loops[i].poly[0];
                !(0..loops.len()).any(|j| j != i && in_polygon(p, &loops[j].poly))
            })
            .collect::<Vec<_>>();
        for &i in &top {
            assert!(loops[i].area > 0.0, "an outermost union loop is clockwise");
        }
        Shapes { shapes, loops, top }
    }

    /// Area enclosed by the curves: the union plus the holes it encloses.
    pub fn filled_area(&self) -> f64 {
        self.top.iter().map(|&i| self.loops[i].area).sum()
    }

    /// True if `p` is enclosed by the curves (inside an outermost loop).
    pub fn filled_contains(&self, p: DVec2) -> bool {
        self.top.iter().any(|&i| in_polygon(p, &self.loops[i].poly))
    }

    /// Distance from `p` to the nearest shape boundary (approximate).
    pub fn boundary_distance(&self, p: DVec2) -> f64 {
        self.shapes
            .iter()
            .map(|s| match s {
                Shape::Poly(poly) => polyline_distance(p, poly),
                Shape::Circle(c, r) => ((p - *c).length() - r).abs(),
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// All boundary curves.
    pub fn curves(&self) -> Vec<ProfileCurve> {
        self.shapes.iter().flat_map(|s| s.curves()).collect()
    }
}

/// Axis-aligned segments with integer end points, and the faces of their
/// arrangement as sets of unit cells.
pub struct CellGrid {
    /// Grid extent: coordinates in `[0, n]`.
    pub n: i64,
    /// Segments `(x0, y0, x1, y1)`, horizontal or vertical.
    pub segs: Vec<(i64, i64, i64, i64)>,
}

impl CellGrid {
    /// Bounded faces as lists of cells `(i, j)` = `[i, i+1] × [j, j+1]`.
    pub fn faces(&self) -> Vec<Vec<(i64, i64)>> {
        let mut vblock: HashSet<(i64, i64)> = HashSet::new(); // x = i, y ∈ [j, j+1]
        let mut hblock: HashSet<(i64, i64)> = HashSet::new(); // y = j, x ∈ [i, i+1]
        for &(x0, y0, x1, y1) in &self.segs {
            if x0 == x1 {
                for j in y0.min(y1)..y0.max(y1) {
                    vblock.insert((x0, j));
                }
            } else {
                for i in x0.min(x1)..x0.max(x1) {
                    hblock.insert((i, y0));
                }
            }
        }
        let lo = -1;
        let hi = self.n;
        let mut comp: std::collections::HashMap<(i64, i64), usize> = Default::default();
        let mut faces = Vec::new();
        let mut outside = usize::MAX;
        for i in lo..=hi {
            for j in lo..=hi {
                if comp.contains_key(&(i, j)) {
                    continue;
                }
                let id = faces.len();
                let mut cells = Vec::new();
                let mut q = VecDeque::from([(i, j)]);
                comp.insert((i, j), id);
                let mut touches = false;
                while let Some((a, b)) = q.pop_front() {
                    cells.push((a, b));
                    if a == lo || b == lo || a == hi || b == hi {
                        touches = true;
                    }
                    let nbrs = [
                        ((a + 1, b), !vblock.contains(&(a + 1, b))),
                        ((a - 1, b), !vblock.contains(&(a, b))),
                        ((a, b + 1), !hblock.contains(&(a, b + 1))),
                        ((a, b - 1), !hblock.contains(&(a, b))),
                    ];
                    for ((x, y), open) in nbrs {
                        if open
                            && (lo..=hi).contains(&x)
                            && (lo..=hi).contains(&y)
                            && !comp.contains_key(&(x, y))
                        {
                            comp.insert((x, y), id);
                            q.push_back((x, y));
                        }
                    }
                }
                if touches {
                    outside = id;
                }
                faces.push(cells);
            }
        }
        faces
            .into_iter()
            .enumerate()
            .filter(|&(id, _)| id != outside)
            .map(|(_, c)| c)
            .collect()
    }
}
