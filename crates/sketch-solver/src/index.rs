//! Validation of a sketch and its translation to dense internal indices
//! (design note 0005 §5 step 1). Every problem is reported as
//! [`SolveError::Invalid`] with a message naming the entity or constraint;
//! nothing here panics on bad input.

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;
use std::fmt::Debug;

use glam::DVec2;

use crate::model::{
    check_text, Constraint, DistanceKind, Geometry, Sketch, SketchKey, TangentKind, RHO_MIN,
};
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

/// An ellipse or elliptical arc over point indices (design note 0008 §3):
/// `E(θ) = c + cos θ·(m − c) + sin θ·(n − c)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Ell {
    /// Centre.
    pub c: usize,
    /// End of the first (major) semi-axis.
    pub m: usize,
    /// End of the second (minor) semi-axis.
    pub n: usize,
    /// Start and end of an elliptical arc (`None` for a full ellipse).
    pub arc: Option<(usize, usize)>,
}

/// A conic arc over point indices (design note 0009 §3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Cnc {
    /// Start point.
    pub p0: usize,
    /// Where the end tangents meet.
    pub p1: usize,
    /// End point.
    pub p2: usize,
    /// Conic index (into the shape unknowns).
    pub k: usize,
}

/// One side of a tangency at a shared point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TanSide {
    /// A line, or a conic's leg at its end: its direction is the
    /// tangent.
    Line(LineRef),
    /// An arc: `point − centre` is the normal.
    Arc {
        /// Centre point index.
        c: usize,
    },
    /// An elliptical arc: its tangent at the shared point.
    Ellipse(Ell),
    /// A spline at its first (`end == false`) or last control point.
    SplineEnd {
        /// Spline index.
        spline: usize,
        /// True for the end at parameter 1.
        end: bool,
    },
}

/// One side of a curvature constraint (design note 0010 §4): the curve
/// and which of its ends is at the joint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CurvSide {
    /// A line; `start` if the joint is at `l.s`.
    Line {
        /// The line.
        l: LineRef,
        /// The joint is at its start.
        start: bool,
    },
    /// An arc on centre `c`.
    Arc {
        /// Centre point index.
        c: usize,
        /// The joint is at its start.
        start: bool,
    },
    /// An elliptical arc.
    Ellipse {
        /// The ellipse.
        el: Ell,
        /// The joint is at its start.
        start: bool,
    },
    /// A conic arc.
    Conic {
        /// The conic.
        cn: Cnc,
        /// The joint is at its start `p0`.
        start: bool,
    },
    /// A spline of degree ≥ 2.
    Spline {
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
    OnEllipse(usize, Ell),
    OnConic(usize, Cnc),
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
    /// A smooth joint with equal curvature at the shared point.
    CurvatureAt([CurvSide; 2], usize),
    TangentLC(LineRef, Circ, bool),
    TangentLE(LineRef, Ell, bool),
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
    Rho(Cnc, f64),
    Mirror(LineRef, Derived),
    /// A linear pattern: the direction line and the copies.
    LinearPattern(LineRef, Pattern),
    /// A circular pattern: the centre point, the angle and the copies.
    CircularPattern(usize, f64, Pattern),
    Offset(OffsetCon),
}

/// A curve of an offset, over internal indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OCurve {
    /// A line.
    Line(LineRef),
    /// An arc: centre, start and end point indices.
    Arc {
        /// Centre.
        c: usize,
        /// Start.
        s: usize,
        /// End.
        e: usize,
    },
    /// A circle: centre point and radius indices.
    Circle {
        /// Centre.
        c: usize,
        /// Radius index.
        r: usize,
    },
}

impl OCurve {
    /// Start and end point indices (none for a circle).
    pub fn ends(&self) -> Option<(usize, usize)> {
        match *self {
            OCurve::Line(l) => Some((l.s, l.e)),
            OCurve::Arc { s, e, .. } => Some((s, e)),
            OCurve::Circle { .. } => None,
        }
    }
}

/// One equation placing an offset's image point (design note 0007 §10.4).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OEq {
    /// On piece `k`'s offset curve.
    OnOffset(usize),
    /// On the normal to piece `k`'s source through its end point `v`.
    Normal(usize, usize),
}

/// A validated offset constraint (design note 0007 §10).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct OffsetCon {
    /// Distance (metres).
    pub distance: f64,
    /// `(source, image, left)` pieces.
    pub pieces: Vec<(OCurve, OCurve, bool)>,
    /// Each image end point with its two equations.
    pub points: Vec<(usize, OEq, OEq)>,
    /// `(source, image)` entity positions.
    pub entities: Vec<(usize, usize)>,
    /// Entity positions of the image arcs, whose implicit equation the
    /// offset's equations replace.
    pub arcs: Vec<usize>,
}

/// The point-level content of a derived constraint (design note 0007 §9):
/// which image unknowns follow which source unknowns.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Derived {
    /// `(source, image)` point pairs, each once, image ≠ source.
    pub points: Vec<(usize, usize)>,
    /// Source points shared with their image that must be held where the
    /// map fixes them (for a mirror: on the axis).
    pub fixed: Vec<usize>,
    /// `(source, image)` circle radius indices.
    pub radii: Vec<(usize, usize)>,
    /// `(source, image)` conic indices, whose shapes are equal.
    pub shapes: Vec<(usize, usize)>,
    /// `(source, image)` entity positions.
    pub entities: Vec<(usize, usize)>,
    /// Entity positions of the image arcs, whose implicit equation the
    /// derived equations imply.
    pub arcs: Vec<usize>,
    /// Entity positions of the image ellipses and elliptical arcs, whose
    /// implicit equations the derived equations imply.
    pub ellipses: Vec<usize>,
    /// Entity positions of the image texts, whose box equations the
    /// derived equations imply.
    pub texts: Vec<usize>,
}

/// A copy turned by less than this from a whole number of turns (radians,
/// the model's angular tolerance) would lie on its source.
const WHOLE_TURN: f64 = 1e-9;

/// A validated pattern's copies (design note 0007 §11.1), over internal
/// indices.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Pattern {
    /// `(source, image, k)` point triples, each image once.
    pub points: Vec<(usize, usize, u32)>,
    /// `(source, image)` circle radius indices.
    pub radii: Vec<(usize, usize)>,
    /// `(source, image)` conic indices, whose shapes are equal.
    pub shapes: Vec<(usize, usize)>,
    /// `(source, image)` entity positions.
    pub entities: Vec<(usize, usize)>,
    /// Entity positions of the image arcs, whose implicit equation the
    /// pattern's equations imply.
    pub arcs: Vec<usize>,
    /// Entity positions of the image ellipses and elliptical arcs, whose
    /// implicit equations the pattern's equations imply.
    pub ellipses: Vec<usize>,
    /// Entity positions of the image texts, whose box equations the
    /// pattern's equations imply.
    pub texts: Vec<usize>,
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
    /// An ellipse or an elliptical arc.
    Ellipse(Ell),
    /// A conic arc.
    Conic(Cnc),
    /// Spline index.
    Spline(usize),
    /// Text index.
    Text(usize),
}

/// A validated text entity (ADR 0015).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TextEnt {
    /// Corner point indices: bottom-left, bottom-right, top-right,
    /// top-left in the text's own frame.
    pub c: [usize; 4],
    /// Box width over height.
    pub aspect: f64,
    /// True for a mirror image.
    pub mirrored: bool,
}

/// A validated spline entity.
#[derive(Clone, Debug)]
pub(crate) struct SplineEnt {
    /// Linear map from points to curve.
    pub map: SplineMap,
    /// Referenced point indices, in order.
    pub points: Vec<usize>,
    /// The through-point parameters of an interpolated spline (`None` for
    /// a control spline).
    pub params: Option<Vec<f64>>,
}

/// A sketch's entities over dense indices.
#[derive(Clone, Debug)]
pub(crate) struct Index<K> {
    /// Entity key → entity position.
    pub ent: BTreeMap<K, usize>,
    /// Entity position → key.
    pub keys: Vec<K>,
    /// Point index → entity position.
    pub points: Vec<usize>,
    /// Position of each point.
    pub positions: Vec<DVec2>,
    /// Circle (radius) index → entity position.
    pub circles: Vec<usize>,
    /// Radius of each circle.
    pub radii: Vec<f64>,
    /// Conic (shape) index → entity position.
    pub conics: Vec<usize>,
    /// Shape value ρ of each conic.
    pub rhos: Vec<f64>,
    /// Spline entities.
    pub splines: Vec<SplineEnt>,
    /// Text entities.
    pub texts: Vec<TextEnt>,
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
        let mut conics = Vec::new();
        let mut rhos = Vec::new();
        let mut splines = Vec::new();
        let mut texts = Vec::new();
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
                Geometry::Ellipse {
                    center,
                    major,
                    minor,
                } => {
                    let el = ellipse_points(id, &point, center, major, minor)?;
                    EntKind::Ellipse(el)
                }
                Geometry::EllipticalArc {
                    center,
                    major,
                    minor,
                    start,
                    end,
                } => {
                    let mut el = ellipse_points(id, &point, center, major, minor)?;
                    let s = point(id, start, "start")?;
                    let t = point(id, end, "end")?;
                    if s == el.c || t == el.c || s == t {
                        return Err(invalid(format!(
                            "entity {id:?}: an elliptical arc needs distinct start and end points, neither its centre"
                        )));
                    }
                    el.arc = Some((s, t));
                    EntKind::Ellipse(el)
                }
                Geometry::Conic {
                    start,
                    control,
                    end,
                    rho,
                } => {
                    let p0 = point(id, start, "start")?;
                    let p1 = point(id, control, "control point")?;
                    let p2 = point(id, end, "end")?;
                    if p0 == p1 || p1 == p2 || p0 == p2 {
                        return Err(invalid(format!(
                            "entity {id:?}: a conic needs distinct start, control and end points"
                        )));
                    }
                    if !(*rho >= RHO_MIN && *rho <= 1.0 - RHO_MIN) {
                        return Err(invalid(format!(
                            "entity {id:?}: a conic's ρ must lie in [{RHO_MIN}, 1 − {RHO_MIN}] (got {rho})"
                        )));
                    }
                    let k = conics.len();
                    conics.push(i);
                    rhos.push(*rho);
                    EntKind::Conic(Cnc { p0, p1, p2, k })
                }
                Geometry::ControlSpline {
                    control,
                    degree,
                    periodic,
                } => {
                    let pts = control
                        .iter()
                        .map(|k| point(id, k, "control point"))
                        .collect::<Result<Vec<_>, _>>()?;
                    let deg = usize::try_from(*degree).unwrap_or(usize::MAX);
                    let map = if *periodic {
                        closes_on_first(id, &pts)?;
                        SplineMap::control_periodic(pts.len(), deg)
                    } else {
                        SplineMap::control(pts.len(), deg)
                    }
                    .map_err(|err| spline_err(id, err))?;
                    splines.push(SplineEnt {
                        map,
                        points: pts,
                        params: None,
                    });
                    EntKind::Spline(splines.len() - 1)
                }
                Geometry::InterpolatedSpline {
                    through,
                    params,
                    periodic,
                } => {
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
                    let map = if *periodic {
                        closes_on_first(id, &pts)?;
                        SplineMap::interpolated_periodic(params)
                    } else {
                        SplineMap::interpolated(params)
                    }
                    .map_err(|err| spline_err(id, err))?;
                    splines.push(SplineEnt {
                        map,
                        points: pts,
                        params: Some(params.clone()),
                    });
                    EntKind::Spline(splines.len() - 1)
                }
                Geometry::Text {
                    corners,
                    aspect,
                    mirrored,
                    outline,
                    ..
                } => {
                    check_text(corners, *aspect, outline)
                        .map_err(|why| invalid(format!("entity {id:?}: {why}")))?;
                    let mut c = [0; 4];
                    for (slot, (k, role)) in c.iter_mut().zip(corners.iter().zip([
                        "bottom-left corner",
                        "bottom-right corner",
                        "top-right corner",
                        "top-left corner",
                    ])) {
                        *slot = point(id, k, role)?;
                    }
                    texts.push(TextEnt {
                        c,
                        aspect: *aspect,
                        mirrored: *mirrored,
                    });
                    EntKind::Text(texts.len() - 1)
                }
            };
            kinds.push(kind);
        }
        Ok(Index {
            ent,
            keys: sketch.entities.iter().map(|e| e.id).collect(),
            points,
            positions,
            circles,
            radii,
            conics,
            rhos,
            splines,
            texts,
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
                    EntKind::Ellipse(el) => {
                        if el.c == p {
                            return Err(invalid(format!(
                                "constraint {id:?}: the centre of {curve:?} cannot lie on it"
                            )));
                        }
                        Con::OnEllipse(p, el)
                    }
                    EntKind::Conic(cn) => {
                        if cn.p1 == p {
                            return Err(invalid(format!(
                                "constraint {id:?}: the control point of {curve:?} cannot lie on it"
                            )));
                        }
                        Con::OnConic(p, cn)
                    }
                    EntKind::Point(_) => {
                        return Err(invalid(format!(
                            "constraint {id:?}: entity {curve:?} must be a curve"
                        )))
                    }
                    EntKind::Text(_) => return Err(text_refused(id, curve)),
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
                        (EntKind::Line(l), EntKind::Ellipse(el))
                        | (EntKind::Ellipse(el), EntKind::Line(l)) => {
                            Con::TangentLE(l, el, *center_left)
                        }
                        (EntKind::Conic(_), _) | (_, EntKind::Conic(_)) => {
                            return Err(conic_tangency(id))
                        }
                        _ => {
                            return Err(invalid(format!(
                                "constraint {id:?}: line-circle tangency needs a line and a circle, arc, ellipse or elliptical arc"
                            )))
                        }
                    },
                    TangentKind::CircleCircle { internal } => match (ka, kb) {
                        (EntKind::Circular(x), EntKind::Circular(y)) => {
                            Con::TangentCC(x, y, *internal)
                        }
                        (EntKind::Conic(_), _) | (_, EntKind::Conic(_)) => {
                            return Err(conic_tangency(id))
                        }
                        (EntKind::Ellipse(_), _) | (_, EntKind::Ellipse(_)) => {
                            return Err(invalid(format!(
                                "constraint {id:?}: tangency of an ellipse with a circle or another ellipse is not supported yet; make them tangent at a shared end point instead"
                            )))
                        }
                        _ => {
                            return Err(invalid(format!(
                                "constraint {id:?}: circle-circle tangency needs two circles or arcs"
                            )))
                        }
                    },
                }
            }
            Constraint::Curvature { a, b, point } => {
                distinct(a, b)?;
                let ka = self.kind_of(id, a)?;
                let kb = self.kind_of(id, b)?;
                let p = self.point(id, point)?;
                let ca = self.curv_side(id, a, ka, p, coincident)?;
                let cb = self.curv_side(id, b, kb, p, coincident)?;
                if !matches!(ca, CurvSide::Spline { .. }) && !matches!(cb, CurvSide::Spline { .. })
                {
                    return Err(invalid(format!(
                        "constraint {id:?}: curvature continuity needs a spline on at least one side"
                    )));
                }
                Con::CurvatureAt([ca, cb], p)
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
            Constraint::Rho { curve, value } => {
                let EntKind::Conic(cn) = self.kind_of(id, curve)? else {
                    return Err(invalid(format!(
                        "constraint {id:?}: entity {curve:?} must be a conic"
                    )));
                };
                if !(*value >= RHO_MIN && *value <= 1.0 - RHO_MIN) {
                    return Err(invalid(format!(
                        "constraint {id:?}: a conic's ρ must lie in [{RHO_MIN}, 1 − {RHO_MIN}] (got {value})"
                    )));
                }
                Con::Rho(cn, *value)
            }
            Constraint::Mirror { axis, pairs } => {
                let l = self.line(id, axis)?;
                if pairs.iter().any(|(s, i)| s == axis || i == axis) {
                    return Err(invalid(format!(
                        "constraint {id:?}: the mirror line {axis:?} cannot be mirrored or be an image"
                    )));
                }
                let mut d = self.derived(id, pairs, true)?;
                if let Some(&(_, q)) = d.points.iter().find(|&&(_, q)| q == l.s || q == l.e) {
                    return Err(invalid(format!(
                        "constraint {id:?}: point {:?} is an end of the mirror line and cannot be an image",
                        self.point_key(q)
                    )));
                }
                // A shared point that is (coincident with) an end of the
                // axis is on it already; holding it there again would be
                // an identically satisfied, redundant equation. Likewise
                // one shared point per coincidence class is held: its
                // coincident partners follow.
                let class = |p: usize| coincident.get(p).copied().unwrap_or(p);
                let mut held = BTreeSet::new();
                d.fixed.retain(|&p| {
                    class(p) != class(l.s) && class(p) != class(l.e) && held.insert(class(p))
                });
                Con::Mirror(l, d)
            }
            Constraint::Offset {
                distance,
                pieces,
                smooth,
            } => Con::Offset(self.offset(
                id,
                length(*distance, "offset distance")?,
                pieces,
                smooth,
                coincident,
            )?),
            Constraint::LinearPattern { direction, copies } => {
                let l = self.line(id, direction)?;
                if copies
                    .iter()
                    .any(|(s, i, _)| s == direction || i == direction)
                {
                    return Err(invalid(format!(
                        "constraint {id:?}: the pattern's line {direction:?} cannot be copied or be a copy"
                    )));
                }
                let pt = self.pattern(id, copies, None, coincident)?;
                if let Some(&(_, q, _)) = pt.points.iter().find(|&&(_, q, _)| q == l.s || q == l.e)
                {
                    return Err(invalid(format!(
                        "constraint {id:?}: point {:?} is an end of the pattern's line and cannot be a copy",
                        self.point_key(q)
                    )));
                }
                Con::LinearPattern(l, pt)
            }
            Constraint::CircularPattern {
                centre,
                angle,
                copies,
            } => {
                let c = self.point(id, centre)?;
                if !(angle.is_finite() && *angle != 0.0 && angle.abs() < TAU) {
                    return Err(invalid(format!(
                        "constraint {id:?}: a pattern angle must be finite, non-zero and less than a turn (got {angle})"
                    )));
                }
                if let Some(&(_, _, k)) = copies.iter().find(|&&(_, _, k)| {
                    let t = (f64::from(k) * angle).rem_euclid(TAU);
                    t < WHOLE_TURN || TAU - t < WHOLE_TURN
                }) {
                    return Err(invalid(format!(
                        "constraint {id:?}: copy {k} would be turned a whole number of turns, onto its source"
                    )));
                }
                if copies.iter().any(|(_, i, _)| i == centre) {
                    return Err(invalid(format!(
                        "constraint {id:?}: the pattern's centre {centre:?} cannot be a copy"
                    )));
                }
                let pt = self.pattern(id, copies, Some(c), coincident)?;
                if pt.points.iter().any(|&(_, q, _)| q == c) {
                    return Err(invalid(format!(
                        "constraint {id:?}: the pattern's centre {centre:?} cannot be a copy"
                    )));
                }
                Con::CircularPattern(c, *angle, pt)
            }
        })
    }

    /// Validates a pattern's copies (design note 0007 §11.1) and maps them
    /// to point triples. Each copy number is a map of its own (a
    /// translation or a rotation, both keeping orientation); a source point
    /// may be its own image only at the `centre` of a circular pattern (or
    /// coincident with it), where it gets no equation.
    fn pattern(
        &self,
        id: &dyn Debug,
        copies: &[(K, K, u32)],
        centre: Option<usize>,
        coincident: &[usize],
    ) -> Result<Pattern, SolveError> {
        if copies.is_empty() {
            return Err(invalid(format!("constraint {id:?}: there are no copies")));
        }
        let mut by_k: BTreeMap<u32, Vec<(K, K)>> = BTreeMap::new();
        for (src, img, k) in copies {
            if *k == 0 {
                return Err(invalid(format!(
                    "constraint {id:?}: the copy {img:?} of {src:?} must be copy 1 or later"
                )));
            }
            by_k.entry(*k).or_default().push((*src, *img));
        }
        let class = |p: usize| coincident.get(p).copied().unwrap_or(p);
        let mut out = Pattern::default();
        // Image point → source point, over every copy.
        let mut source_of: BTreeMap<usize, usize> = BTreeMap::new();
        for (k, pairs) in by_k {
            let d = self.derived(id, &pairs, false)?;
            if let Some(&p) = d
                .fixed
                .iter()
                .find(|&&p| centre.is_none_or(|c| class(p) != class(c)))
            {
                return Err(invalid(format!(
                    "constraint {id:?}: point {:?} cannot be its own copy (only the centre of a circular pattern can)",
                    self.point_key(p)
                )));
            }
            for (p, q) in d.points {
                if source_of.insert(q, p).is_some() {
                    return Err(invalid(format!(
                        "constraint {id:?}: point {:?} is in more than one copy",
                        self.point_key(q)
                    )));
                }
                out.points.push((p, q, k));
            }
            out.radii.extend(d.radii);
            out.shapes.extend(d.shapes);
            out.entities.extend(d.entities);
            out.arcs.extend(d.arcs);
            out.ellipses.extend(d.ellipses);
            out.texts.extend(d.texts);
        }
        let mut images = BTreeSet::new();
        if let Some(&(_, i)) = out.entities.iter().find(|&&(_, i)| !images.insert(i)) {
            return Err(invalid(format!(
                "constraint {id:?}: entity {:?} is in more than one copy",
                self.keys[i]
            )));
        }
        if let Some(&(_, q, _)) = out
            .points
            .iter()
            .find(|&&(_, q, _)| source_of.values().any(|&p| p == q))
        {
            return Err(invalid(format!(
                "constraint {id:?}: point {:?} is both a copy and a source",
                self.point_key(q)
            )));
        }
        Ok(out)
    }

    /// Validates an offset's pieces and gives each image end point its two
    /// equations (design note 0007 §10.4): on its piece's offset curve, and
    /// on the next piece's (a corner), or on the normal at its source's end
    /// (a free end, or a smooth joint).
    fn offset(
        &self,
        id: &dyn Debug,
        distance: f64,
        pieces: &[(K, K, bool)],
        smooth: &[K],
        coincident: &[usize],
    ) -> Result<OffsetCon, SolveError> {
        if pieces.is_empty() {
            return Err(invalid(format!(
                "constraint {id:?}: there are no curves to offset"
            )));
        }
        let mut curves = Vec::with_capacity(pieces.len());
        let mut entities: Vec<(usize, usize)> = Vec::with_capacity(pieces.len());
        let mut arcs = Vec::new();
        let mut source_points = BTreeSet::new();
        // Image end point → (piece, true for the piece's end).
        let mut roles: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
        for (k, (src, img, left)) in pieces.iter().enumerate() {
            if src == img {
                return Err(invalid(format!(
                    "constraint {id:?}: entity {src:?} cannot be its own offset"
                )));
            }
            let (ks, ki) = (self.kind_of(id, src)?, self.kind_of(id, img)?);
            let (es, ei) = (self.ent[src], self.ent[img]);
            if entities.iter().any(|&(_, i)| i == ei) {
                return Err(invalid(format!(
                    "constraint {id:?}: entity {img:?} is the offset of more than one curve"
                )));
            }
            let same_centre = |c: usize, c2: usize| {
                if c == c2 {
                    Ok(())
                } else {
                    Err(invalid(format!(
                        "constraint {id:?}: the offset {img:?} of {src:?} must share its centre point"
                    )))
                }
            };
            let (a, b) = match (ks, ki) {
                (EntKind::Line(a), EntKind::Line(b)) => (OCurve::Line(a), OCurve::Line(b)),
                (
                    EntKind::Circular(Circ::Arc { c, s, e }),
                    EntKind::Circular(Circ::Arc {
                        c: c2,
                        s: s2,
                        e: e2,
                    }),
                ) => {
                    same_centre(c, c2)?;
                    arcs.push(ei);
                    (
                        OCurve::Arc { c, s, e },
                        OCurve::Arc {
                            c: c2,
                            s: s2,
                            e: e2,
                        },
                    )
                }
                (
                    EntKind::Circular(Circ::Circle { c, r }),
                    EntKind::Circular(Circ::Circle { c: c2, r: r2 }),
                ) => {
                    same_centre(c, c2)?;
                    (OCurve::Circle { c, r }, OCurve::Circle { c: c2, r: r2 })
                }
                _ => {
                    return Err(invalid(format!(
                        "constraint {id:?}: {src:?} and its offset {img:?} are not two lines, two arcs or two circles"
                    )))
                }
            };
            match a {
                OCurve::Line(l) => source_points.extend([l.s, l.e]),
                OCurve::Arc { c, s, e } => source_points.extend([c, s, e]),
                OCurve::Circle { c, .. } => {
                    source_points.insert(c);
                }
            }
            if let Some((s, e)) = b.ends() {
                roles.entry(s).or_default().push((k, false));
                roles.entry(e).or_default().push((k, true));
            }
            entities.push((es, ei));
            curves.push((a, b, *left));
        }
        let smooth: BTreeSet<usize> = smooth
            .iter()
            .map(|k| self.point(id, k))
            .collect::<Result<_, _>>()?;
        if let Some(&p) = smooth.iter().find(|p| !roles.contains_key(p)) {
            return Err(invalid(format!(
                "constraint {id:?}: point {:?} is not an end of an offset curve",
                self.point_key(p)
            )));
        }
        let class = |p: usize| coincident.get(p).copied().unwrap_or(p);
        // The source point behind an image end.
        let behind = |k: usize, at_end: bool| -> usize {
            let (s, e) = curves[k]
                .0
                .ends()
                .expect("an image with ends has a source with ends");
            if at_end {
                e
            } else {
                s
            }
        };
        let mut points = Vec::with_capacity(roles.len());
        for (&p, r) in &roles {
            if source_points.contains(&p) {
                return Err(invalid(format!(
                    "constraint {id:?}: point {:?} is both an image and a source of the offset",
                    self.point_key(p)
                )));
            }
            let (k, at_end) = r[0];
            let v = behind(k, at_end);
            let (first, second) = match r.as_slice() {
                [_] => (OEq::OnOffset(k), OEq::Normal(k, v)),
                [_, (k2, at_end2)] => {
                    if class(v) != class(behind(*k2, *at_end2)) {
                        return Err(invalid(format!(
                            "constraint {id:?}: the offsets {:?} and {:?} meet where their sources do not",
                            pieces[k].1, pieces[*k2].1
                        )));
                    }
                    if smooth.contains(&p) {
                        (OEq::OnOffset(k), OEq::Normal(k, v))
                    } else {
                        (OEq::OnOffset(k), OEq::OnOffset(*k2))
                    }
                }
                _ => {
                    return Err(invalid(format!(
                        "constraint {id:?}: point {:?} ends more than two offset curves",
                        self.point_key(p)
                    )))
                }
            };
            points.push((p, first, second));
        }
        Ok(OffsetCon {
            distance,
            pieces: curves,
            points,
            entities,
            arcs,
        })
    }

    /// The key of point index `p`.
    fn point_key(&self, p: usize) -> K {
        self.keys[self.points[p]]
    }

    /// Validates the `(source, image)` entity pairs of a derived constraint
    /// (design note 0007 §9) and maps them to point pairs. `reverses`: the
    /// map reverses the sense of rotation (a mirror), so an image arc
    /// starts at the image of its source's end; a mirror is also its own
    /// inverse, so a point pair met both ways round is one relation.
    fn derived(
        &self,
        id: &dyn Debug,
        pairs: &[(K, K)],
        reverses: bool,
    ) -> Result<Derived, SolveError> {
        if pairs.is_empty() {
            return Err(invalid(format!(
                "constraint {id:?}: there are no entities to derive"
            )));
        }
        let mut d = Derived::default();
        let mut point_pairs: Vec<(usize, usize)> = Vec::new();
        for (src, img) in pairs {
            if src == img {
                return Err(invalid(format!(
                    "constraint {id:?}: entity {src:?} cannot be its own image"
                )));
            }
            let (ks, ki) = (self.kind_of(id, src)?, self.kind_of(id, img)?);
            let (es, ei) = (self.ent[src], self.ent[img]);
            if d.entities.iter().any(|&(_, i)| i == ei) {
                return Err(invalid(format!(
                    "constraint {id:?}: entity {img:?} is the image of more than one entity"
                )));
            }
            match (ks, ki) {
                (EntKind::Point(p), EntKind::Point(q)) => point_pairs.push((p, q)),
                (EntKind::Line(a), EntKind::Line(b)) => {
                    point_pairs.extend([(a.s, b.s), (a.e, b.e)]);
                }
                (
                    EntKind::Circular(Circ::Circle { c, r }),
                    EntKind::Circular(Circ::Circle { c: c2, r: r2 }),
                ) => {
                    point_pairs.push((c, c2));
                    d.radii.push((r, r2));
                }
                (
                    EntKind::Circular(Circ::Arc { c, s, e }),
                    EntKind::Circular(Circ::Arc {
                        c: c2,
                        s: s2,
                        e: e2,
                    }),
                ) => {
                    let (s, e) = if reverses { (e, s) } else { (s, e) };
                    point_pairs.extend([(c, c2), (s, s2), (e, e2)]);
                    d.arcs.push(ei);
                }
                (EntKind::Ellipse(a), EntKind::Ellipse(b)) if a.arc.is_some() == b.arc.is_some() => {
                    point_pairs.extend([(a.c, b.c), (a.m, b.m), (a.n, b.n)]);
                    if let (Some((s, e)), Some((s2, e2))) = (a.arc, b.arc) {
                        let (s, e) = if reverses { (e, s) } else { (s, e) };
                        point_pairs.extend([(s, s2), (e, e2)]);
                    }
                    d.ellipses.push(ei);
                }
                // A conic has no orientation: start maps to start even
                // under a reflection (design note 0009 §4.5).
                (EntKind::Conic(a), EntKind::Conic(b)) => {
                    point_pairs.extend([(a.p0, b.p0), (a.p1, b.p1), (a.p2, b.p2)]);
                    d.shapes.push((a.k, b.k));
                }
                (EntKind::Spline(a), EntKind::Spline(b)) => {
                    let (a, b) = (&self.splines[a], &self.splines[b]);
                    if a.map.degree != b.map.degree
                        || a.map.periodic != b.map.periodic
                        || a.points.len() != b.points.len()
                        || a.params != b.params
                    {
                        return Err(invalid(format!(
                            "constraint {id:?}: {src:?} and its image {img:?} are not splines of the same form"
                        )));
                    }
                    point_pairs.extend(a.points.iter().copied().zip(b.points.iter().copied()));
                }
                // Corners keep their names in the text's own frame; a
                // mirror image of text is mirrored text (ADR 0015).
                (EntKind::Text(a), EntKind::Text(b)) => {
                    let (a, b) = (self.texts[a], self.texts[b]);
                    if a.aspect != b.aspect || b.mirrored != (a.mirrored != reverses) {
                        return Err(invalid(format!(
                            "constraint {id:?}: {img:?} is not {} copy of the text {src:?}",
                            if reverses { "a mirrored" } else { "an unmirrored" }
                        )));
                    }
                    point_pairs.extend(a.c.into_iter().zip(b.c));
                    d.texts.push(ei);
                }
                _ => {
                    return Err(invalid(format!(
                        "constraint {id:?}: {src:?} and its image {img:?} are not the same kind of entity"
                    )))
                }
            }
            d.entities.push((es, ei));
        }
        // Each image point follows exactly one source point; a pair met
        // twice (an end shared by two curves) counts once, and so does a
        // pair met both ways round when the map is its own inverse (a
        // curve mirrored from a point of an earlier image back onto that
        // image's source).
        let mut source_of: BTreeMap<usize, usize> = BTreeMap::new();
        for (p, q) in point_pairs {
            if reverses && source_of.get(&p) == Some(&q) {
                continue;
            }
            match source_of.insert(q, p) {
                Some(p0) if p0 != p => {
                    return Err(invalid(format!(
                        "constraint {id:?}: point {:?} would be the image of both {:?} and {:?}",
                        self.point_key(q),
                        self.point_key(p0),
                        self.point_key(p)
                    )))
                }
                Some(_) => continue,
                None => {}
            }
            if p == q {
                d.fixed.push(p);
            } else {
                d.points.push((p, q));
            }
        }
        let sources: BTreeSet<usize> = source_of.values().copied().collect();
        if let Some(&(_, q)) = d.points.iter().find(|(_, q)| sources.contains(q)) {
            return Err(invalid(format!(
                "constraint {id:?}: point {:?} is both an image and a source",
                self.point_key(q)
            )));
        }
        Ok(d)
    }

    /// One side of a curvature constraint at point `p` (design note 0010
    /// §3): which end of `curve` is at `p`.
    fn curv_side(
        &self,
        id: &dyn Debug,
        curve: &K,
        kind: EntKind,
        p: usize,
        coincident: &[usize],
    ) -> Result<CurvSide, SolveError> {
        let at = |q: usize| q == p || coincident[q] == coincident[p];
        let end = |s: usize, e: usize| {
            match (at(s), at(e)) {
            (true, false) => Ok(true),
            (false, true) => Ok(false),
            (true, true) => Err(invalid(format!(
                "constraint {id:?}: both ends of {curve:?} are at the joint; the end is ambiguous"
            ))),
            (false, false) => Err(invalid(format!(
                "constraint {id:?}: the joint is not an end point of {curve:?} (nor coincident with one)"
            ))),
        }
        };
        Ok(match kind {
            EntKind::Line(l) => CurvSide::Line {
                l,
                start: end(l.s, l.e)?,
            },
            EntKind::Circular(Circ::Arc { c, s, e }) => CurvSide::Arc {
                c,
                start: end(s, e)?,
            },
            EntKind::Ellipse(el) => match el.arc {
                Some((s, e)) => CurvSide::Ellipse {
                    el,
                    start: end(s, e)?,
                },
                None => {
                    return Err(invalid(format!(
                        "constraint {id:?}: {curve:?} is an ellipse and has no end points"
                    )))
                }
            },
            EntKind::Conic(cn) => CurvSide::Conic {
                cn,
                start: end(cn.p0, cn.p2)?,
            },
            EntKind::Spline(s) if self.splines[s].map.periodic => {
                return Err(closed_refused(id, curve))
            }
            EntKind::Spline(s) => {
                let def = &self.splines[s];
                if def.map.degree < 2 {
                    return Err(invalid(format!(
                        "constraint {id:?}: {curve:?} has degree 1 and is straight at its ends; curvature continuity needs degree 2 or more"
                    )));
                }
                let pts = &def.points;
                CurvSide::Spline {
                    spline: s,
                    end: !end(pts[0], pts[pts.len() - 1])?,
                }
            }
            EntKind::Circular(Circ::Circle { .. }) => {
                return Err(invalid(format!(
                    "constraint {id:?}: {curve:?} is a circle and has no end points"
                )))
            }
            EntKind::Point(_) => {
                return Err(invalid(format!(
                    "constraint {id:?}: entity {curve:?} must be a curve"
                )))
            }
            EntKind::Text(_) => return Err(text_refused(id, curve)),
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
            EntKind::Ellipse(el) => match el.arc {
                Some((s, e)) if at(s) || at(e) => Ok(TanSide::Ellipse(el)),
                Some(_) => Err(not_end()),
                None => Err(invalid(format!(
                    "constraint {id:?}: {curve:?} is an ellipse and has no end points; use line-ellipse tangency"
                ))),
            },
            EntKind::Conic(cn) => match (at(cn.p0), at(cn.p2)) {
                // The tangent at an end is along the leg to the control
                // point (design note 0009 §4.3).
                (true, false) => Ok(TanSide::Line(LineRef { s: cn.p0, e: cn.p1 })),
                (false, true) => Ok(TanSide::Line(LineRef { s: cn.p1, e: cn.p2 })),
                (true, true) => Err(invalid(format!(
                    "constraint {id:?}: both ends of {curve:?} are at the tangency point; the end is ambiguous"
                ))),
                (false, false) => Err(not_end()),
            },
            EntKind::Spline(s) if self.splines[s].map.periodic => Err(closed_refused(id, curve)),
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
            EntKind::Text(_) => Err(text_refused(id, curve)),
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
        let cons = sketch
            .constraints
            .iter()
            .map(|c| self.constraint_dbg(&c.id, &c.constraint, &classes))
            .collect::<Result<Vec<Con>, _>>()?;
        self.check_derived(sketch, &cons)?;
        Ok(cons)
    }

    /// Checks the derived constraints together: an entity is the image of
    /// at most one, and no entity is derived, through any chain, from its
    /// own image (an image arc's implicit equation is left out because its
    /// source's implies it, which needs every chain to end at a source
    /// that is no image).
    fn check_derived(&self, sketch: &Sketch<K>, cons: &[Con]) -> Result<(), SolveError> {
        // Image entity → (source entity, constraint position).
        let mut source_of: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
        for (j, con) in cons.iter().enumerate() {
            let entities = match con {
                Con::Mirror(_, d) => &d.entities,
                Con::Offset(o) => &o.entities,
                Con::LinearPattern(_, pt) | Con::CircularPattern(_, _, pt) => &pt.entities,
                _ => continue,
            };
            for &(s, i) in entities {
                if let Some((_, j0)) = source_of.insert(i, (s, j)) {
                    return Err(invalid(format!(
                        "constraint {:?}: entity {:?} is already an image of constraint {:?}",
                        sketch.constraints[j].id, self.keys[i], sketch.constraints[j0].id
                    )));
                }
            }
        }
        for (&image, &(mut x, j)) in &source_of {
            for _ in 0..source_of.len() {
                if x == image {
                    return Err(invalid(format!(
                        "constraint {:?}: entity {:?} is derived from its own image",
                        sketch.constraints[j].id, self.keys[image]
                    )));
                }
                match source_of.get(&x) {
                    Some(&(next, _)) => x = next,
                    None => break,
                }
            }
        }
        Ok(())
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

/// The centre and axis ends of an ellipse entity, distinct.
fn ellipse_points<K: SketchKey>(
    id: &K,
    point: &dyn Fn(&K, &K, &str) -> Result<usize, SolveError>,
    center: &K,
    major: &K,
    minor: &K,
) -> Result<Ell, SolveError> {
    let c = point(id, center, "centre")?;
    let m = point(id, major, "major axis end")?;
    let n = point(id, minor, "minor axis end")?;
    if c == m || c == n || m == n {
        return Err(invalid(format!(
            "entity {id:?}: an ellipse needs distinct centre and axis end points"
        )));
    }
    Ok(Ell { c, m, n, arc: None })
}

/// The refusal of tangency to a conic away from its ends (design note
/// 0009 §4.7).
fn conic_tangency(id: &dyn Debug) -> SolveError {
    invalid(format!(
        "constraint {id:?}: tangency of a conic away from its ends is not supported yet; make them tangent at a shared end point instead"
    ))
}

/// Text takes part in constraints only through its box's corners and lines
/// (ADR 0015).
fn text_refused(id: &dyn Debug, curve: &dyn Debug) -> SolveError {
    invalid(format!(
        "constraint {id:?}: {curve:?} is text; constrain its box's corners or lines instead"
    ))
}

/// Checks that a smooth closed spline's point list ends on its first
/// point (ADR 0016).
fn closes_on_first<K: SketchKey>(id: &K, pts: &[usize]) -> Result<(), SolveError> {
    if pts.len() < 2 || pts[0] != pts[pts.len() - 1] {
        return Err(invalid(format!(
            "entity {id:?}: a smooth closed spline must end on its first point"
        )));
    }
    Ok(())
}

/// A tangency or curvature joint at a smooth closed spline, which has no
/// end points (ADR 0016).
fn closed_refused(id: &dyn Debug, curve: &dyn Debug) -> SolveError {
    invalid(format!(
        "constraint {id:?}: {curve:?} is a smooth closed spline and has no end points"
    ))
}

fn spline_err<K: SketchKey>(id: &K, e: SplineError) -> SolveError {
    match e {
        SplineError::Shape(m) | SplineError::Params(m) => invalid(format!("entity {id:?}: {m}")),
    }
}
