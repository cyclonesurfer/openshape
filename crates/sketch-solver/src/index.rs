//! Validation of a sketch and its translation to dense internal indices
//! (design note 0005 §5 step 1). Every problem is reported as
//! [`SolveError::Invalid`] with a message naming the entity or constraint;
//! nothing here panics on bad input.

use std::collections::BTreeMap;
use std::fmt::Debug;

use glam::DVec2;

use crate::model::{Constraint, DistanceKind, Geometry, Sketch, SketchKey, TangentKind};
use crate::spline::{SplineError, SplineMap};
use crate::SolveError;

/// Largest accepted magnitude of a coordinate, radius or length value
/// (metres). Beyond it the model tolerance (1e-7 m) is below what `f64`
/// can resolve after a solve to 1e-12 of the sketch size.
pub(crate) const MAX_MAGNITUDE: f64 = 1e6;

/// A line by the point indices of its ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct LineRef {
    /// Start point index.
    pub s: usize,
    /// End point index.
    pub e: usize,
}

/// A circle or arc.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Circ {
    /// A circle: centre point index and radius index.
    Circle {
        /// Centre.
        c: usize,
        /// Radius index (into the radius unknowns).
        r: usize,
    },
    /// An arc: centre, start and end point indices.
    Arc {
        /// Centre.
        c: usize,
        /// Start.
        s: usize,
        /// End.
        e: usize,
    },
}

impl Circ {
    /// Centre point index.
    pub fn center(&self) -> usize {
        match *self {
            Circ::Circle { c, .. } | Circ::Arc { c, .. } => c,
        }
    }
}

/// One side of a tangency at a shared point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TanSide {
    /// A line: its direction is the tangent.
    Line(LineRef),
    /// An arc: `point − centre` is the normal.
    Arc {
        /// Centre point index.
        c: usize,
    },
    /// A spline at its first (`end == false`) or last control point.
    SplineEnd {
        /// Spline index.
        spline: usize,
        /// True for the end at parameter 1.
        end: bool,
    },
}

/// A validated constraint over internal indices.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Con {
    Coincident(usize, usize),
    OnLine(usize, LineRef),
    OnCirc(usize, Circ),
    OnSpline(usize, usize),
    Horizontal(LineRef),
    Vertical(LineRef),
    HPoints(usize, usize),
    VPoints(usize, usize),
    Parallel(LineRef, LineRef),
    Perpendicular(LineRef, LineRef),
    EqualLines(LineRef, LineRef),
    EqualCirc(Circ, Circ),
    Midpoint(usize, LineRef),
    Concentric(Circ, Circ),
    Symmetric(usize, usize, LineRef),
    TangentAt(TanSide, TanSide, usize),
    TangentLC(LineRef, Circ, bool),
    TangentCC(Circ, Circ, bool),
    Normal(LineRef, Circ),
    Fix(usize, DVec2),
    DistPP(usize, usize, DistanceKind, f64),
    DistPL(usize, LineRef, f64),
    DistLL(LineRef, LineRef, f64),
    Length(LineRef, f64),
    Angle(LineRef, LineRef, f64),
    Radius(Circ, f64),
    Diameter(Circ, f64),
}

/// What an entity is, over internal indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum EntKind {
    /// Point index.
    Point(usize),
    /// A line.
    Line(LineRef),
    /// A circle or an arc.
    Circular(Circ),
    /// Spline index.
    Spline(usize),
}

/// A validated spline entity.
#[derive(Clone, Debug)]
pub(crate) struct SplineEnt {
    /// Linear map from points to curve.
    pub map: SplineMap,
    /// Referenced point indices, in order.
    pub points: Vec<usize>,
}

/// A sketch's entities over dense indices.
#[derive(Clone, Debug)]
pub(crate) struct Index<K> {
    /// Entity key → entity position.
    pub ent: BTreeMap<K, usize>,
    /// Point index → entity position.
    pub points: Vec<usize>,
    /// Position of each point.
    pub positions: Vec<DVec2>,
    /// Circle (radius) index → entity position.
    pub circles: Vec<usize>,
    /// Radius of each circle.
    pub radii: Vec<f64>,
    /// Spline entities.
    pub splines: Vec<SplineEnt>,
    /// Entity position → kind.
    pub kinds: Vec<EntKind>,
}

fn invalid(msg: String) -> SolveError {
    SolveError::Invalid(msg)
}

fn finite_vec(v: DVec2) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.x.abs() <= MAX_MAGNITUDE && v.y.abs() <= MAX_MAGNITUDE
}

impl<K: SketchKey> Index<K> {
    /// Validates the entities of `sketch`.
    pub fn new(sketch: &Sketch<K>) -> Result<Self, SolveError> {
        let mut ent = BTreeMap::new();
        for (i, e) in sketch.entities.iter().enumerate() {
            if ent.insert(e.id, i).is_some() {
                return Err(invalid(format!("duplicate entity key {:?}", e.id)));
            }
        }
        // Points first, so curves can refer to them.
        let n = sketch.entities.len();
        let mut point_of = vec![None; n];
        let mut points = Vec::new();
        let mut positions = Vec::new();
        for (i, e) in sketch.entities.iter().enumerate() {
            if let Geometry::Point { at } = e.geometry {
                if !finite_vec(at) {
                    return Err(invalid(format!(
                        "entity {:?}: point coordinates must be finite and within ±{MAX_MAGNITUDE} m",
                        e.id
                    )));
                }
                point_of[i] = Some(points.len());
                points.push(i);
                positions.push(at);
            }
        }
        let point = |owner: &K, k: &K, role: &str| -> Result<usize, SolveError> {
            match ent.get(k) {
                None => Err(invalid(format!(
                    "entity {owner:?}: {role} {k:?} does not exist"
                ))),
                Some(&i) => point_of[i].ok_or_else(|| {
                    invalid(format!("entity {owner:?}: {role} {k:?} is not a point"))
                }),
            }
        };
        let mut kinds = Vec::with_capacity(n);
        let mut circles = Vec::new();
        let mut radii = Vec::new();
        let mut splines = Vec::new();
        for (i, e) in sketch.entities.iter().enumerate() {
            let id = &e.id;
            let kind = match &e.geometry {
                Geometry::Point { .. } => EntKind::Point(point_of[i].expect("point")),
                Geometry::Line { start, end } => {
                    let s = point(id, start, "start")?;
                    let t = point(id, end, "end")?;
                    if s == t {
                        return Err(invalid(format!(
                            "entity {id:?}: a line cannot start and end at the same point"
                        )));
                    }
                    EntKind::Line(LineRef { s, e: t })
                }
                Geometry::Circle { center, radius } => {
                    let c = point(id, center, "centre")?;
                    if !radius.is_finite() || *radius <= 0.0 || *radius > MAX_MAGNITUDE {
                        return Err(invalid(format!(
                            "entity {id:?}: circle radius must be positive, finite and at most {MAX_MAGNITUDE} m (got {radius})"
                        )));
                    }
                    let r = circles.len();
                    circles.push(i);
                    radii.push(*radius);
                    EntKind::Circular(Circ::Circle { c, r })
                }
                Geometry::Arc { center, start, end } => {
                    let c = point(id, center, "centre")?;
                    let s = point(id, start, "start")?;
                    let t = point(id, end, "end")?;
                    if c == s || c == t || s == t {
                        return Err(invalid(format!(
                            "entity {id:?}: an arc needs distinct centre, start and end points"
                        )));
                    }
                    EntKind::Circular(Circ::Arc { c, s, e: t })
                }
                Geometry::ControlSpline { control, degree } => {
                    let pts = control
                        .iter()
                        .map(|k| point(id, k, "control point"))
                        .collect::<Result<Vec<_>, _>>()?;
                    let deg = usize::try_from(*degree).unwrap_or(usize::MAX);
                    let map =
                        SplineMap::control(pts.len(), deg).map_err(|err| spline_err(id, err))?;
                    splines.push(SplineEnt { map, points: pts });
                    EntKind::Spline(splines.len() - 1)
                }
                Geometry::InterpolatedSpline { through, params } => {
                    let pts = through
                        .iter()
                        .map(|k| point(id, k, "through point"))
                        .collect::<Result<Vec<_>, _>>()?;
                    if pts.len() < 2 {
                        return Err(invalid(format!(
                            "entity {id:?}: an interpolated spline needs at least 2 points"
                        )));
                    }
                    if params.len() != pts.len() {
                        return Err(invalid(format!(
                            "entity {id:?}: {} parameters for {} through points",
                            params.len(),
                            pts.len()
                        )));
                    }
                    let map = SplineMap::interpolated(params).map_err(|err| spline_err(id, err))?;
                    splines.push(SplineEnt { map, points: pts });
                    EntKind::Spline(splines.len() - 1)
                }
            };
            kinds.push(kind);
        }
        Ok(Index {
            ent,
            points,
            positions,
            circles,
            radii,
            splines,
            kinds,
        })
    }

    fn kind_of(&self, owner: &dyn Debug, k: &K) -> Result<EntKind, SolveError> {
        self.ent
            .get(k)
            .map(|&i| self.kinds[i])
            .ok_or_else(|| invalid(format!("constraint {owner:?}: entity {k:?} does not exist")))
    }

    fn point(&self, owner: &dyn Debug, k: &K) -> Result<usize, SolveError> {
        match self.kind_of(owner, k)? {
            EntKind::Point(p) => Ok(p),
            _ => Err(invalid(format!(
                "constraint {owner:?}: entity {k:?} must be a point"
            ))),
        }
    }

    fn line(&self, owner: &dyn Debug, k: &K) -> Result<LineRef, SolveError> {
        match self.kind_of(owner, k)? {
            EntKind::Line(l) => Ok(l),
            _ => Err(invalid(format!(
                "constraint {owner:?}: entity {k:?} must be a line"
            ))),
        }
    }

    fn circular(&self, owner: &dyn Debug, k: &K) -> Result<Circ, SolveError> {
        match self.kind_of(owner, k)? {
            EntKind::Circular(c) => Ok(c),
            _ => Err(invalid(format!(
                "constraint {owner:?}: entity {k:?} must be a circle or an arc"
            ))),
        }
    }

    /// The radius of circle index `r`.
    pub fn radius(&self, r: usize) -> f64 {
        self.radii[r]
    }

    /// Validates one constraint (`id` labels it in messages). `coincident` gives, for each point index,
    /// a representative of its class under all coincident constraints.
    pub fn constraint_dbg(
        &self,
        id: &dyn Debug,
        c: &Constraint<K>,
        coincident: &[usize],
    ) -> Result<Con, SolveError> {
        let distinct = |a: &K, b: &K| -> Result<(), SolveError> {
            if a == b {
                Err(invalid(format!(
                    "constraint {id:?}: both references are the same entity {a:?}"
                )))
            } else {
                Ok(())
            }
        };
        let length = |v: f64, what: &str| -> Result<f64, SolveError> {
            if !v.is_finite() {
                Err(invalid(format!("constraint {id:?}: {what} must be finite")))
            } else if v == 0.0 {
                Err(invalid(format!(
                    "constraint {id:?}: a zero {what} is not allowed; use a coincident constraint instead"
                )))
            } else if v < 0.0 {
                Err(invalid(format!(
                    "constraint {id:?}: {what} must be positive (got {v})"
                )))
            } else if v > MAX_MAGNITUDE {
                Err(invalid(format!(
                    "constraint {id:?}: {what} exceeds {MAX_MAGNITUDE} m"
                )))
            } else {
                Ok(v)
            }
        };
        Ok(match c {
            Constraint::Coincident { a, b } => {
                distinct(a, b)?;
                Con::Coincident(self.point(id, a)?, self.point(id, b)?)
            }
            Constraint::PointOnCurve { point, curve } => {
                let p = self.point(id, point)?;
                match self.kind_of(id, curve)? {
                    EntKind::Line(l) => Con::OnLine(p, l),
                    EntKind::Circular(cc) => {
                        if cc.center() == p {
                            return Err(invalid(format!(
                                "constraint {id:?}: the centre of {curve:?} cannot lie on it"
                            )));
                        }
                        Con::OnCirc(p, cc)
                    }
                    EntKind::Spline(s) => Con::OnSpline(p, s),
                    EntKind::Point(_) => {
                        return Err(invalid(format!(
                            "constraint {id:?}: entity {curve:?} must be a curve"
                        )))
                    }
                }
            }
            Constraint::Horizontal { line } => Con::Horizontal(self.line(id, line)?),
            Constraint::Vertical { line } => Con::Vertical(self.line(id, line)?),
            Constraint::HorizontalPoints { a, b } => {
                distinct(a, b)?;
                Con::HPoints(self.point(id, a)?, self.point(id, b)?)
            }
            Constraint::VerticalPoints { a, b } => {
                distinct(a, b)?;
                Con::VPoints(self.point(id, a)?, self.point(id, b)?)
            }
            Constraint::Parallel { a, b } => {
                distinct(a, b)?;
                Con::Parallel(self.line(id, a)?, self.line(id, b)?)
            }
            Constraint::Perpendicular { a, b } => {
                distinct(a, b)?;
                Con::Perpendicular(self.line(id, a)?, self.line(id, b)?)
            }
            Constraint::Equal { a, b } => {
                distinct(a, b)?;
                match (self.kind_of(id, a)?, self.kind_of(id, b)?) {
                    (EntKind::Line(x), EntKind::Line(y)) => Con::EqualLines(x, y),
                    (EntKind::Circular(x), EntKind::Circular(y)) => Con::EqualCirc(x, y),
                    _ => {
                        return Err(invalid(format!(
                            "constraint {id:?}: equal needs two lines or two circles/arcs"
                        )))
                    }
                }
            }
            Constraint::Midpoint { point, line } => {
                let p = self.point(id, point)?;
                let l = self.line(id, line)?;
                if p == l.s || p == l.e {
                    return Err(invalid(format!(
                        "constraint {id:?}: an end point of {line:?} cannot be its midpoint"
                    )));
                }
                Con::Midpoint(p, l)
            }
            Constraint::Concentric { a, b } => {
                distinct(a, b)?;
                Con::Concentric(self.circular(id, a)?, self.circular(id, b)?)
            }
            Constraint::Symmetric { a, b, axis } => {
                distinct(a, b)?;
                Con::Symmetric(self.point(id, a)?, self.point(id, b)?, self.line(id, axis)?)
            }
            Constraint::Tangent { a, b, kind } => {
                distinct(a, b)?;
                let ka = self.kind_of(id, a)?;
                let kb = self.kind_of(id, b)?;
                match kind {
                    TangentKind::AtPoint { point } => {
                        let p = self.point(id, point)?;
                        let sa = self.tan_side(id, a, ka, p, coincident)?;
                        let sb = self.tan_side(id, b, kb, p, coincident)?;
                        Con::TangentAt(sa, sb, p)
                    }
                    TangentKind::LineCircle { center_left } => match (ka, kb) {
                        (EntKind::Line(l), EntKind::Circular(cc))
                        | (EntKind::Circular(cc), EntKind::Line(l)) => {
                            Con::TangentLC(l, cc, *center_left)
                        }
                        _ => {
                            return Err(invalid(format!(
                                "constraint {id:?}: line-circle tangency needs a line and a circle or arc"
                            )))
                        }
                    },
                    TangentKind::CircleCircle { internal } => match (ka, kb) {
                        (EntKind::Circular(x), EntKind::Circular(y)) => {
                            Con::TangentCC(x, y, *internal)
                        }
                        _ => {
                            return Err(invalid(format!(
                                "constraint {id:?}: circle-circle tangency needs two circles or arcs"
                            )))
                        }
                    },
                }
            }
            Constraint::Normal { line, curve } => {
                Con::Normal(self.line(id, line)?, self.circular(id, curve)?)
            }
            Constraint::Fix { point, at } => {
                if !finite_vec(*at) {
                    return Err(invalid(format!(
                        "constraint {id:?}: the fixed position must be finite and within ±{MAX_MAGNITUDE} m"
                    )));
                }
                Con::Fix(self.point(id, point)?, *at)
            }
            Constraint::Distance { a, b, kind, value } => {
                distinct(a, b)?;
                let v = length(*value, "distance")?;
                match (self.kind_of(id, a)?, self.kind_of(id, b)?, kind) {
                    (EntKind::Point(p), EntKind::Point(q), _) => Con::DistPP(p, q, *kind, v),
                    (EntKind::Point(p), EntKind::Line(l), DistanceKind::Aligned)
                    | (EntKind::Line(l), EntKind::Point(p), DistanceKind::Aligned) => {
                        Con::DistPL(p, l, v)
                    }
                    (EntKind::Line(x), EntKind::Line(y), DistanceKind::Aligned) => {
                        Con::DistLL(x, y, v)
                    }
                    (_, _, DistanceKind::Horizontal | DistanceKind::Vertical) => {
                        return Err(invalid(format!(
                            "constraint {id:?}: horizontal and vertical distances need two points"
                        )))
                    }
                    _ => {
                        return Err(invalid(format!(
                            "constraint {id:?}: a distance needs two points, a point and a line, or two lines"
                        )))
                    }
                }
            }
            Constraint::Length { line, value } => {
                Con::Length(self.line(id, line)?, length(*value, "length")?)
            }
            Constraint::Angle { a, b, value } => {
                distinct(a, b)?;
                let pi = std::f64::consts::PI;
                if !value.is_finite() || *value <= -pi || *value > pi {
                    return Err(invalid(format!(
                        "constraint {id:?}: an angle must lie in (−π, π] (got {value})"
                    )));
                }
                Con::Angle(self.line(id, a)?, self.line(id, b)?, *value)
            }
            Constraint::Radius { curve, value } => {
                Con::Radius(self.circular(id, curve)?, length(*value, "radius")?)
            }
            Constraint::Diameter { curve, value } => {
                Con::Diameter(self.circular(id, curve)?, length(*value, "diameter")?)
            }
        })
    }

    fn tan_side(
        &self,
        id: &dyn Debug,
        curve: &K,
        kind: EntKind,
        p: usize,
        coincident: &[usize],
    ) -> Result<TanSide, SolveError> {
        let at = |q: usize| q == p || coincident[q] == coincident[p];
        let not_end = || {
            invalid(format!(
                "constraint {id:?}: the tangency point is not an end point of {curve:?} (nor coincident with one)"
            ))
        };
        match kind {
            EntKind::Line(l) => {
                if at(l.s) || at(l.e) {
                    Ok(TanSide::Line(l))
                } else {
                    Err(not_end())
                }
            }
            EntKind::Circular(Circ::Arc { c, s, e }) => {
                if at(s) || at(e) {
                    Ok(TanSide::Arc { c })
                } else {
                    Err(not_end())
                }
            }
            EntKind::Circular(Circ::Circle { .. }) => Err(invalid(format!(
                "constraint {id:?}: {curve:?} is a circle and has no end points; use line-circle or circle-circle tangency"
            ))),
            EntKind::Spline(s) => {
                let pts = &self.splines[s].points;
                let first = at(pts[0]);
                let last = at(pts[pts.len() - 1]);
                match (first, last) {
                    (true, false) => Ok(TanSide::SplineEnd {
                        spline: s,
                        end: false,
                    }),
                    (false, true) => Ok(TanSide::SplineEnd {
                        spline: s,
                        end: true,
                    }),
                    (true, true) => Err(invalid(format!(
                        "constraint {id:?}: both ends of {curve:?} are at the tangency point; the end is ambiguous"
                    ))),
                    (false, false) => Err(not_end()),
                }
            }
            EntKind::Point(_) => Err(invalid(format!(
                "constraint {id:?}: entity {curve:?} must be a curve"
            ))),
        }
    }

    /// Validates every constraint, in order. Returns `(position, Con)`.
    pub fn constraints(&self, sketch: &Sketch<K>) -> Result<Vec<Con>, SolveError> {
        let mut seen = BTreeMap::new();
        for c in &sketch.constraints {
            if seen.insert(c.id, ()).is_some() {
                return Err(invalid(format!("duplicate constraint key {:?}", c.id)));
            }
        }
        // Coincidence classes for the tangency end-point rule.
        let mut uf: Vec<usize> = (0..self.points.len()).collect();
        fn find(uf: &mut [usize], mut x: usize) -> usize {
            while uf[x] != x {
                uf[x] = uf[uf[x]];
                x = uf[x];
            }
            x
        }
        for c in &sketch.constraints {
            if let Constraint::Coincident { a, b } = &c.constraint {
                if let (Ok(p), Ok(q)) = (self.point(&c.id, a), self.point(&c.id, b)) {
                    let (rp, rq) = (find(&mut uf, p), find(&mut uf, q));
                    if rp != rq {
                        let (lo, hi) = (rp.min(rq), rp.max(rq));
                        uf[hi] = lo;
                    }
                }
            }
        }
        let classes: Vec<usize> = (0..self.points.len()).map(|p| find(&mut uf, p)).collect();
        sketch
            .constraints
            .iter()
            .map(|c| self.constraint_dbg(&c.id, &c.constraint, &classes))
            .collect()
    }

    /// Validates drag targets: `(point index, target)`.
    pub fn drag(&self, drag: &[(K, DVec2)]) -> Result<Vec<(usize, DVec2)>, SolveError> {
        let mut seen = BTreeMap::new();
        let mut out = Vec::with_capacity(drag.len());
        for (k, at) in drag {
            if seen.insert(*k, ()).is_some() {
                return Err(invalid(format!("drag: point {k:?} is listed twice")));
            }
            let p = match self.ent.get(k).map(|&i| self.kinds[i]) {
                Some(EntKind::Point(p)) => p,
                Some(_) => return Err(invalid(format!("drag: entity {k:?} is not a point"))),
                None => return Err(invalid(format!("drag: entity {k:?} does not exist"))),
            };
            if !finite_vec(*at) {
                return Err(invalid(format!(
                    "drag: the target of {k:?} must be finite and within ±{MAX_MAGNITUDE} m"
                )));
            }
            out.push((p, *at));
        }
        Ok(out)
    }
}

fn spline_err<K: SketchKey>(id: &K, e: SplineError) -> SolveError {
    match e {
        SplineError::Shape(m) | SplineError::Params(m) => invalid(format!("entity {id:?}: {m}")),
    }
}
