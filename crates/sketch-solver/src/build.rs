//! Builds the residual system of a validated sketch (design note 0005 §3–4):
//! unknown slots, scaling, and one or more equations per constraint, in
//! priority order (implicit arc equations first, then constraints in order).

use glam::DVec2;

use crate::conic;
use crate::eq::{CurvTerm, EqKind, PairMode, RTerm, Rhs, SplineDef, TanDir, VecExpr};
use crate::index::{
    Circ, Cnc, Con, CurvSide, Ell, Index, LineRef, OCurve, OEq, OffsetCon, TanSide,
};
use crate::model::{DistanceKind, SketchKey};

/// Lengths are divided by `l` after subtracting the centre `(cx, cy)` from
/// coordinates (radii are only divided).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Scale {
    /// Length unit (metres).
    pub l: f64,
    /// Centre x (metres).
    pub cx: f64,
    /// Centre y (metres).
    pub cy: f64,
}

impl Scale {
    /// No scaling: equations in metres and radians.
    pub fn identity() -> Self {
        Scale {
            l: 1.0,
            cx: 0.0,
            cy: 0.0,
        }
    }

    /// A scaled length.
    pub fn len(&self, v: f64) -> f64 {
        v / self.l
    }

    /// A scaled position.
    pub fn pos(&self, p: DVec2) -> (f64, f64) {
        ((p.x - self.cx) / self.l, (p.y - self.cy) / self.l)
    }
}

/// An equation that can be eliminated by substitution (a pure equality of
/// two unknowns, or an unknown equal to a constant).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Sub {
    /// `slot a = slot b`.
    Equal(usize, usize),
    /// `slot = value` (scaled), with the unscaled value for write-back.
    Const {
        /// Slot.
        slot: usize,
        /// Scaled value.
        value: f64,
        /// Value in metres (for a conic's shape, which is not scaled, the
        /// same as `value`).
        abs: f64,
    },
}

/// One residual equation.
#[derive(Clone, Debug)]
pub(crate) struct Equation {
    /// The residual.
    pub kind: EqKind,
    /// Constraint position, or `None` for an implicit equation of an arc,
    /// ellipse or elliptical arc.
    pub owner: Option<usize>,
    /// Entity position of the arc or ellipse, for implicit equations.
    pub arc: Option<usize>,
    /// True if the residual is an angle (sine/cosine), false for lengths.
    pub angular: bool,
    /// Elimination form, if any.
    pub sub: Option<Sub>,
}

/// The curve a hidden parameter runs along.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ParamCurve {
    /// A spline (by index): the parameter lies in `[0, 1]`.
    Spline(usize),
    /// An ellipse: the parameter is its eccentric anomaly, periodic and
    /// unbounded (design note 0008 §4.2).
    Ellipse(Ell),
    /// A conic arc: the parameter lies in `[0, 1]` (design note 0009 §4.2).
    Conic(Cnc),
}

/// A hidden curve parameter placing a point on a curve: a point-on-spline,
/// point-on-ellipse or point-on-conic constraint, or an elliptical arc's
/// own end.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Param {
    /// The curve.
    pub curve: ParamCurve,
    /// The point on it.
    pub point: usize,
}

/// The residual system.
#[derive(Clone, Debug)]
pub(crate) struct Problem {
    /// Number of point entities.
    pub n_points: usize,
    /// Number of circle radii.
    pub n_circles: usize,
    /// Number of conic shapes `u = ln w` (design note 0009 §4.1).
    pub n_shapes: usize,
    /// Hidden curve parameters, in slot order.
    pub params: Vec<Param>,
    /// Equations in priority order.
    pub eqs: Vec<Equation>,
    /// Splines.
    pub splines: Vec<SplineDef>,
    /// Directions shorter than this (in scaled units) are degenerate.
    pub degen: f64,
}

/// The slot layout: two per point, then circle radii, conic shapes and
/// hidden curve parameters.
impl Problem {
    /// Total number of slots.
    pub fn n_slots(&self) -> usize {
        self.param_base() + self.params.len()
    }

    /// Slot of circle radius `r`.
    pub fn radius_slot(&self, r: usize) -> usize {
        2 * self.n_points + r
    }

    /// Slot of conic shape `k`.
    pub fn shape_slot(&self, k: usize) -> usize {
        2 * self.n_points + self.n_circles + k
    }

    /// Slot of hidden parameter `k`.
    pub fn param_slot(&self, k: usize) -> usize {
        self.param_base() + k
    }

    fn param_base(&self) -> usize {
        2 * self.n_points + self.n_circles + self.n_shapes
    }

    /// True for hidden-parameter slots.
    pub fn is_param(&self, slot: usize) -> bool {
        slot >= self.param_base()
    }

    /// True for hidden parameters confined to `[0, 1]` (on open splines and
    /// conics; the parameters of an ellipse and of a smooth closed spline
    /// are periodic and unbounded).
    pub fn is_bounded_param(&self, slot: usize) -> bool {
        self.is_param(slot)
            && match self.params[slot - self.param_base()].curve {
                ParamCurve::Ellipse(_) => false,
                ParamCurve::Spline(s) => !self.splines[s].map.periodic,
                ParamCurve::Conic(_) => true,
            }
    }

    /// True for circle-radius slots.
    pub fn is_radius(&self, slot: usize) -> bool {
        slot >= 2 * self.n_points && slot < 2 * self.n_points + self.n_circles
    }

    /// True for conic-shape slots.
    pub fn is_shape(&self, slot: usize) -> bool {
        slot >= 2 * self.n_points + self.n_circles && slot < self.param_base()
    }
}

/// Relative threshold below which a direction is degenerate.
pub(crate) const DEGENERATE: f64 = 1e-10;

/// The scaling of design note 0005 §3: `L` is the larger of the bounding-box
/// diagonal (points, circle extents, fixed positions and drag targets) and
/// the largest length dimension, or 1e-6 m if both are zero; the centre is
/// the bounding-box centre.
pub(crate) fn scale_of<K: SketchKey>(
    index: &Index<K>,
    cons: &[Con],
    drag: &[(usize, DVec2)],
) -> Scale {
    let mut lo = DVec2::splat(f64::INFINITY);
    let mut hi = DVec2::splat(f64::NEG_INFINITY);
    let mut add = |p: DVec2| {
        lo = lo.min(p);
        hi = hi.max(p);
    };
    for p in &index.positions {
        add(*p);
    }
    for k in &index.kinds {
        match *k {
            crate::index::EntKind::Circular(Circ::Circle { c, r }) => {
                let center = index.positions[c];
                let rad = index.radius(r);
                add(center - DVec2::splat(rad));
                add(center + DVec2::splat(rad));
            }
            crate::index::EntKind::Ellipse(el) => {
                // The vertices opposite the axis ends.
                let p = &index.positions;
                add(p[el.c] * 2.0 - p[el.m]);
                add(p[el.c] * 2.0 - p[el.n]);
            }
            _ => {}
        }
    }
    let mut max_dim = 0.0f64;
    for c in cons {
        match c {
            Con::Fix(_, at) => add(*at),
            Con::DistPP(_, _, _, v)
            | Con::DistPL(_, _, v)
            | Con::DistLL(_, _, v)
            | Con::Length(_, v)
            | Con::Radius(_, v)
            | Con::Diameter(_, v) => max_dim = max_dim.max(*v),
            Con::Offset(o) => max_dim = max_dim.max(o.distance),
            _ => {}
        }
    }
    for (_, t) in drag {
        add(*t);
    }
    let (diag, center) = if lo.x <= hi.x {
        ((hi - lo).length(), (lo + hi) * 0.5)
    } else {
        (0.0, DVec2::ZERO)
    };
    let mut l = diag.max(max_dim);
    if !(l.is_finite() && l > 0.0) {
        l = 1e-6;
    }
    Scale {
        l,
        cx: center.x,
        cy: center.y,
    }
}

struct Builder<'a> {
    scale: Scale,
    radius_base: usize,
    shape_base: usize,
    param_base: usize,
    splines: &'a [SplineDef],
    eqs: Vec<Equation>,
    params: Vec<Param>,
}

impl Builder<'_> {
    fn x(&self, p: usize) -> usize {
        2 * p
    }
    fn y(&self, p: usize) -> usize {
        2 * p + 1
    }

    fn push(&mut self, owner: usize, kind: EqKind, angular: bool, sub: Option<Sub>) {
        self.eqs.push(Equation {
            kind,
            owner: Some(owner),
            arc: None,
            angular,
            sub,
        });
    }

    fn equal_slots(&mut self, owner: usize, a: usize, b: usize) {
        self.push(
            owner,
            EqKind::Lin {
                terms: vec![(a, 1.0), (b, -1.0)],
                rhs: 0.0,
            },
            false,
            Some(Sub::Equal(a, b)),
        );
    }

    fn dir(l: LineRef) -> VecExpr {
        VecExpr::diff(l.e, l.s)
    }

    fn radius(&self, c: Circ) -> RTerm {
        match c {
            Circ::Circle { r, .. } => RTerm::Slot(self.radius_base + r),
            Circ::Arc { c, s, .. } => RTerm::Norm(VecExpr::diff(s, c)),
        }
    }

    /// A new hidden parameter placing `p` on `curve`; returns its slot.
    fn param(&mut self, curve: ParamCurve, p: usize) -> usize {
        self.params.push(Param { curve, point: p });
        self.param_base + self.params.len() - 1
    }

    /// A hidden parameter placing `p` on ellipse `el`, with its two
    /// equations (design note 0008 §4.2), owned by constraint `owner` or,
    /// for an elliptical arc's own ends, implicit to entity `ent`.
    fn on_ellipse(&mut self, owner: Option<usize>, ent: Option<usize>, p: usize, el: Ell) {
        let t = self.param_base + self.params.len();
        self.params.push(Param {
            curve: ParamCurve::Ellipse(el),
            point: p,
        });
        for axis in 0..2 {
            self.eqs.push(Equation {
                kind: EqKind::EllipsePoint {
                    c: el.c,
                    m: el.m,
                    n: el.n,
                    t,
                    p,
                    axis,
                },
                owner,
                arc: ent,
                angular: false,
                sub: None,
            });
        }
    }

    /// The direction of one side of a tangency at point `p`, and whether it
    /// is a tangent (true) or a normal.
    fn tan_dir(&self, side: TanSide, p: usize) -> (TanDir, bool) {
        match side {
            TanSide::Ellipse(el) => (
                TanDir::Ellipse {
                    c: el.c,
                    m: el.m,
                    n: el.n,
                    p,
                },
                true,
            ),
            TanSide::Line(l) => (TanDir::Lin(Self::dir(l)), true),
            TanSide::Arc { c } => (TanDir::Lin(VecExpr::diff(p, c)), false),
            TanSide::SplineEnd { spline, end } => {
                let def = &self.splines[spline];
                let n = def.map.count;
                let (j0, j1) = if end { (n - 2, n - 1) } else { (0, 1) };
                let w0 = def.map.control_weights(j0);
                let w1 = def.map.control_weights(j1);
                let terms = def
                    .points
                    .iter()
                    .enumerate()
                    .filter_map(|(k, &q)| {
                        let c = w1[k] - w0[k];
                        (c != 0.0).then_some((q, c))
                    })
                    .collect();
                (TanDir::Lin(VecExpr(terms)), true)
            }
        }
    }

    /// The tangency row at shared point `p` (design note 0005 §4).
    fn tangent_at(&mut self, owner: usize, a: TanSide, b: TanSide, p: usize) {
        let (da, ta) = self.tan_dir(a, p);
        let (db, tb) = self.tan_dir(b, p);
        // Tangent vs tangent and normal vs normal: parallel (cross);
        // tangent vs normal: perpendicular (dot).
        let mode = if ta == tb {
            PairMode::Cross
        } else {
            PairMode::Dot
        };
        let kind = match (da, db) {
            (TanDir::Lin(a), TanDir::Lin(b)) => EqKind::DirPair { a, b, mode },
            (a, b) => EqKind::TanDirs { a, b, mode },
        };
        self.push(owner, kind, true, None);
    }

    /// One side of the curvature row at shared point `p` (design note 0010
    /// §4).
    fn curv_term(&self, side: CurvSide, p: usize) -> CurvTerm {
        let sign = |away: bool| if away { 1.0 } else { -1.0 };
        match side {
            CurvSide::Line { l, start } => CurvTerm {
                d1: TanDir::Lin(Self::dir(l)),
                sign: sign(start),
                d2: VecExpr(Vec::new()),
                shape: None,
            },
            CurvSide::Arc { c, start } => CurvTerm {
                d1: TanDir::Perp(VecExpr::diff(p, c)),
                sign: sign(start),
                d2: VecExpr::diff(c, p),
                shape: None,
            },
            CurvSide::Ellipse { el, start } => CurvTerm {
                d1: TanDir::Ellipse {
                    c: el.c,
                    m: el.m,
                    n: el.n,
                    p,
                },
                sign: sign(start),
                d2: VecExpr::diff(el.c, p),
                shape: None,
            },
            CurvSide::Conic { cn, start } => {
                let (near, far) = if start {
                    (cn.p0, cn.p2)
                } else {
                    (cn.p2, cn.p0)
                };
                CurvTerm {
                    d1: TanDir::Lin(VecExpr::diff(cn.p1, near)),
                    sign: 1.0,
                    d2: VecExpr::diff(far, cn.p1),
                    shape: Some(self.shape_base + cn.k),
                }
            }
            CurvSide::Spline { spline, end } => CurvTerm::spline(&self.splines[spline], end),
        }
    }

    /// `q` is the mirror image of `p` about line `l`: their midpoint is on
    /// the line and `q − p` is normal to it.
    fn symmetric(&mut self, owner: usize, p: usize, q: usize, l: LineRef) {
        self.push(
            owner,
            EqKind::SignedDist {
                dir: Self::dir(l),
                v: VecExpr(vec![(p, 0.5), (q, 0.5), (l.s, -1.0)]),
                abs: false,
                rhs: Rhs::zero(),
            },
            false,
            None,
        );
        self.push(
            owner,
            EqKind::Proj {
                dir: Self::dir(l),
                v: VecExpr::diff(q, p),
            },
            false,
            None,
        );
    }

    /// The equation of design note 0007 §10.4 placing image point `x` of
    /// offset `o`.
    fn offset_eq(&mut self, owner: usize, o: &OffsetCon, x: usize, eq: OEq) {
        match eq {
            OEq::OnOffset(k) => {
                let (src, _, left) = o.pieces[k];
                let side = if left { 1.0 } else { -1.0 };
                let d = self.scale.len(o.distance);
                let kind = match src {
                    // Left of a line is positive signed distance.
                    OCurve::Line(l) => EqKind::SignedDist {
                        dir: Self::dir(l),
                        v: VecExpr::diff(x, l.s),
                        abs: false,
                        rhs: Rhs::constant(side * d),
                    },
                    // Left of a counter-clockwise arc is towards its centre.
                    OCurve::Arc { c, s, .. } => EqKind::Norm {
                        k: 1.0,
                        v: VecExpr::diff(x, c),
                        rhs: Rhs {
                            terms: vec![
                                (1.0, RTerm::Norm(VecExpr::diff(s, c))),
                                (-side, RTerm::Const(d)),
                            ],
                            abs: false,
                        },
                    },
                    OCurve::Circle { .. } => unreachable!("a circle's offset has no end points"),
                };
                self.push(owner, kind, false, None);
            }
            OEq::Normal(k, v) => {
                let kind = match o.pieces[k].0 {
                    OCurve::Line(l) => EqKind::Proj {
                        dir: Self::dir(l),
                        v: VecExpr::diff(x, v),
                    },
                    OCurve::Arc { c, .. } => EqKind::SignedDist {
                        dir: VecExpr::diff(v, c),
                        v: VecExpr::diff(x, c),
                        abs: false,
                        rhs: Rhs::zero(),
                    },
                    OCurve::Circle { .. } => unreachable!("a circle's offset has no end points"),
                };
                self.push(owner, kind, false, None);
            }
        }
    }

    fn constraint(&mut self, owner: usize, con: &Con) {
        let sc = self.scale;
        match *con {
            Con::Coincident(p, q) => {
                self.equal_slots(owner, self.x(p), self.x(q));
                self.equal_slots(owner, self.y(p), self.y(q));
            }
            Con::OnLine(p, l) => self.push(
                owner,
                EqKind::SignedDist {
                    dir: Self::dir(l),
                    v: VecExpr::diff(p, l.s),
                    abs: false,
                    rhs: Rhs::zero(),
                },
                false,
                None,
            ),
            Con::OnCirc(p, c) => {
                let r = self.radius(c);
                self.push(
                    owner,
                    EqKind::Norm {
                        k: 1.0,
                        v: VecExpr::diff(p, c.center()),
                        rhs: Rhs::term(1.0, r),
                    },
                    false,
                    None,
                )
            }
            Con::OnEllipse(p, el) => self.on_ellipse(Some(owner), None, p, el),
            Con::OnConic(p, cn) => {
                let t = self.param(ParamCurve::Conic(cn), p);
                let u = self.shape_base + cn.k;
                for axis in 0..2 {
                    self.push(
                        owner,
                        EqKind::ConicPoint {
                            p0: cn.p0,
                            p1: cn.p1,
                            p2: cn.p2,
                            u,
                            t,
                            p,
                            axis,
                        },
                        false,
                        None,
                    );
                }
            }
            Con::OnSpline(p, s) => {
                let t = self.param(ParamCurve::Spline(s), p);
                for axis in 0..2 {
                    self.push(
                        owner,
                        EqKind::SplinePoint {
                            spline: s,
                            t,
                            p,
                            axis,
                        },
                        false,
                        None,
                    );
                }
            }
            Con::Horizontal(l) => self.equal_slots(owner, self.y(l.s), self.y(l.e)),
            Con::Vertical(l) => self.equal_slots(owner, self.x(l.s), self.x(l.e)),
            Con::HPoints(a, b) => self.equal_slots(owner, self.y(a), self.y(b)),
            Con::VPoints(a, b) => self.equal_slots(owner, self.x(a), self.x(b)),
            Con::Parallel(a, b) => self.push(
                owner,
                EqKind::DirPair {
                    a: Self::dir(a),
                    b: Self::dir(b),
                    mode: PairMode::Cross,
                },
                true,
                None,
            ),
            Con::Perpendicular(a, b) => self.push(
                owner,
                EqKind::DirPair {
                    a: Self::dir(a),
                    b: Self::dir(b),
                    mode: PairMode::Dot,
                },
                true,
                None,
            ),
            Con::EqualLines(a, b) => self.push(
                owner,
                EqKind::Norm {
                    k: 1.0,
                    v: Self::dir(a),
                    rhs: Rhs::term(1.0, RTerm::Norm(Self::dir(b))),
                },
                false,
                None,
            ),
            Con::EqualCirc(a, b) => match (a, b) {
                (Circ::Circle { r: ra, .. }, Circ::Circle { r: rb, .. }) => {
                    self.equal_slots(owner, self.radius_base + ra, self.radius_base + rb)
                }
                (Circ::Arc { c, s, .. }, other) => {
                    // ‖s − c‖ − R_b
                    let rb = self.radius(other);
                    self.push(
                        owner,
                        EqKind::Norm {
                            k: 1.0,
                            v: VecExpr::diff(s, c),
                            rhs: Rhs::term(1.0, rb),
                        },
                        false,
                        None,
                    )
                }
                (circle @ Circ::Circle { .. }, Circ::Arc { c, s, .. }) => {
                    // r_a − ‖s_b − c_b‖ = −‖s_b − c_b‖ − (−r_a)
                    let ra = self.radius(circle);
                    self.push(
                        owner,
                        EqKind::Norm {
                            k: -1.0,
                            v: VecExpr::diff(s, c),
                            rhs: Rhs::term(-1.0, ra),
                        },
                        false,
                        None,
                    )
                }
            },
            Con::Midpoint(p, l) => {
                for axis in 0..2 {
                    self.push(
                        owner,
                        EqKind::Lin {
                            terms: vec![
                                (2 * p + axis, 1.0),
                                (2 * l.s + axis, -0.5),
                                (2 * l.e + axis, -0.5),
                            ],
                            rhs: 0.0,
                        },
                        false,
                        None,
                    );
                }
            }
            Con::Concentric(a, b) => {
                let (ca, cb) = (a.center(), b.center());
                self.equal_slots(owner, self.x(ca), self.x(cb));
                self.equal_slots(owner, self.y(ca), self.y(cb));
            }
            Con::Symmetric(p, q, l) => self.symmetric(owner, p, q, l),
            Con::TangentAt(a, b, p) => self.tangent_at(owner, a, b, p),
            Con::CurvatureAt([ca, cb], p) => {
                let (a, b) = (self.curv_term(ca, p), self.curv_term(cb, p));
                self.push(
                    owner,
                    EqKind::AwayAngle {
                        a: a.clone(),
                        b: b.clone(),
                        value: std::f64::consts::PI,
                    },
                    true,
                    None,
                );
                self.push(owner, EqKind::Curvature { a, b }, false, None);
            }
            Con::TangentLC(l, c, left) => {
                let r = self.radius(c);
                let sign = if left { 1.0 } else { -1.0 };
                self.push(
                    owner,
                    EqKind::SignedDist {
                        dir: Self::dir(l),
                        v: VecExpr::diff(c.center(), l.s),
                        abs: false,
                        rhs: Rhs::term(sign, r),
                    },
                    false,
                    None,
                );
            }
            Con::TangentLE(l, el, left) => self.push(
                owner,
                EqKind::EllipseSupport {
                    dir: Self::dir(l),
                    v: VecExpr::diff(el.c, l.s),
                    a: VecExpr::diff(el.m, el.c),
                    b: VecExpr::diff(el.n, el.c),
                    sign: if left { 1.0 } else { -1.0 },
                },
                false,
                None,
            ),
            Con::TangentCC(a, b, internal) => {
                let (ra, rb) = (self.radius(a), self.radius(b));
                let rhs = if internal {
                    Rhs {
                        terms: vec![(1.0, ra), (-1.0, rb)],
                        abs: true,
                    }
                } else {
                    Rhs {
                        terms: vec![(1.0, ra), (1.0, rb)],
                        abs: false,
                    }
                };
                self.push(
                    owner,
                    EqKind::Norm {
                        k: 1.0,
                        v: VecExpr::diff(a.center(), b.center()),
                        rhs,
                    },
                    false,
                    None,
                );
            }
            Con::Normal(l, c) => self.push(
                owner,
                EqKind::SignedDist {
                    dir: Self::dir(l),
                    v: VecExpr::diff(c.center(), l.s),
                    abs: false,
                    rhs: Rhs::zero(),
                },
                false,
                None,
            ),
            Con::Fix(p, at) => {
                let (sx, sy) = sc.pos(at);
                for (slot, value, abs) in [(self.x(p), sx, at.x), (self.y(p), sy, at.y)] {
                    self.push(
                        owner,
                        EqKind::Lin {
                            terms: vec![(slot, 1.0)],
                            rhs: value,
                        },
                        false,
                        Some(Sub::Const { slot, value, abs }),
                    );
                }
            }
            Con::DistPP(p, q, kind, v) => {
                let v = sc.len(v);
                match kind {
                    DistanceKind::Aligned => self.push(
                        owner,
                        EqKind::Norm {
                            k: 1.0,
                            v: VecExpr::diff(q, p),
                            rhs: Rhs::constant(v),
                        },
                        false,
                        None,
                    ),
                    DistanceKind::Horizontal | DistanceKind::Vertical => {
                        let axis = usize::from(kind == DistanceKind::Vertical);
                        self.push(
                            owner,
                            EqKind::AbsLin {
                                terms: vec![(2 * q + axis, 1.0), (2 * p + axis, -1.0)],
                                rhs: v,
                            },
                            false,
                            None,
                        )
                    }
                }
            }
            Con::DistPL(p, l, v) => self.push(
                owner,
                EqKind::SignedDist {
                    dir: Self::dir(l),
                    v: VecExpr::diff(p, l.s),
                    abs: true,
                    rhs: Rhs::constant(sc.len(v)),
                },
                false,
                None,
            ),
            Con::DistLL(a, b, v) => self.push(
                owner,
                EqKind::SignedDist {
                    dir: Self::dir(b),
                    v: VecExpr::diff(a.s, b.s),
                    abs: true,
                    rhs: Rhs::constant(sc.len(v)),
                },
                false,
                None,
            ),
            Con::Length(l, v) => self.push(
                owner,
                EqKind::Norm {
                    k: 1.0,
                    v: Self::dir(l),
                    rhs: Rhs::constant(sc.len(v)),
                },
                false,
                None,
            ),
            Con::Angle(a, b, v) => self.push(
                owner,
                EqKind::DirPair {
                    a: Self::dir(a),
                    b: Self::dir(b),
                    mode: PairMode::AngleSin { value: v },
                },
                true,
                None,
            ),
            Con::LinearPattern(l, ref pt) => {
                // Design note 0007 §11.1: q − p − k·(e − s) = 0 per axis,
                // the identity in the image point q.
                for &(p, q, k) in &pt.points {
                    let k = f64::from(k);
                    for (q, p, e, s) in [
                        (self.x(q), self.x(p), self.x(l.e), self.x(l.s)),
                        (self.y(q), self.y(p), self.y(l.e), self.y(l.s)),
                    ] {
                        self.push(
                            owner,
                            EqKind::Lin {
                                terms: vec![(q, 1.0), (p, -1.0), (e, -k), (s, k)],
                                rhs: 0.0,
                            },
                            false,
                            None,
                        );
                    }
                }
                for &(r, r2) in &pt.radii {
                    self.equal_slots(owner, self.radius_base + r, self.radius_base + r2);
                }
                for &(k, k2) in &pt.shapes {
                    self.equal_slots(owner, self.shape_base + k, self.shape_base + k2);
                }
            }
            Con::CircularPattern(c, angle, ref pt) => {
                // Design note 0007 §11.1: q − c − R(k·θ)(p − c) = 0, linear
                // for the given angle and a rotation in q.
                for &(p, q, k) in &pt.points {
                    let (sn, cs) = (f64::from(k) * angle).sin_cos();
                    let (xp, yp, xq, yq, xc, yc) = (
                        self.x(p),
                        self.y(p),
                        self.x(q),
                        self.y(q),
                        self.x(c),
                        self.y(c),
                    );
                    for terms in [
                        vec![(xq, 1.0), (xc, cs - 1.0), (xp, -cs), (yp, sn), (yc, -sn)],
                        vec![(yq, 1.0), (yc, cs - 1.0), (xp, -sn), (xc, sn), (yp, -cs)],
                    ] {
                        self.push(owner, EqKind::Lin { terms, rhs: 0.0 }, false, None);
                    }
                }
                for &(r, r2) in &pt.radii {
                    self.equal_slots(owner, self.radius_base + r, self.radius_base + r2);
                }
                for &(k, k2) in &pt.shapes {
                    self.equal_slots(owner, self.shape_base + k, self.shape_base + k2);
                }
            }
            Con::Mirror(l, ref d) => {
                // Design note 0007 §9: two equations per image point, whose
                // Jacobian block in that point (half the axis normal, the
                // axis direction) is regular whatever else holds it.
                for &(p, q) in &d.points {
                    self.symmetric(owner, p, q, l);
                }
                for &p in &d.fixed {
                    self.push(
                        owner,
                        EqKind::SignedDist {
                            dir: Self::dir(l),
                            v: VecExpr::diff(p, l.s),
                            abs: false,
                            rhs: Rhs::zero(),
                        },
                        false,
                        None,
                    );
                }
                for &(r, r2) in &d.radii {
                    self.equal_slots(owner, self.radius_base + r, self.radius_base + r2);
                }
                for &(k, k2) in &d.shapes {
                    self.equal_slots(owner, self.shape_base + k, self.shape_base + k2);
                }
            }
            Con::Offset(ref o) => {
                // Design note 0007 §10.4: two equations per image end
                // point, regular in it (two offset curves crossing, or an
                // offset curve and a normal to it); an image circle's
                // radius is its source's less or more the distance.
                for &(x, a, b) in &o.points {
                    self.offset_eq(owner, o, x, a);
                    self.offset_eq(owner, o, x, b);
                }
                for &(src, img, left) in &o.pieces {
                    if let (OCurve::Circle { r, .. }, OCurve::Circle { r: r2, .. }) = (src, img) {
                        let side = if left { 1.0 } else { -1.0 };
                        self.push(
                            owner,
                            EqKind::Lin {
                                terms: vec![
                                    (self.radius_base + r2, 1.0),
                                    (self.radius_base + r, -1.0),
                                ],
                                rhs: -side * sc.len(o.distance),
                            },
                            false,
                            None,
                        );
                    }
                }
            }
            Con::Radius(c, v) | Con::Diameter(c, v) => {
                let diameter = matches!(con, Con::Diameter(..));
                let k = if diameter { 2.0 } else { 1.0 };
                let sv = sc.len(v);
                match c {
                    Circ::Circle { r, .. } => {
                        let slot = self.radius_base + r;
                        self.push(
                            owner,
                            EqKind::Lin {
                                terms: vec![(slot, k)],
                                rhs: sv,
                            },
                            false,
                            Some(Sub::Const {
                                slot,
                                value: sv / k,
                                abs: v / k,
                            }),
                        )
                    }
                    Circ::Arc { c, s, .. } => self.push(
                        owner,
                        EqKind::Norm {
                            k,
                            v: VecExpr::diff(s, c),
                            rhs: Rhs::constant(sv),
                        },
                        false,
                        None,
                    ),
                }
            }
            Con::Rho(cn, v) => {
                // u = logit(ρ), linear in the unscaled shape unknown, so it
                // is eliminated like a radius dimension (design note 0009
                // §4.4).
                let slot = self.shape_base + cn.k;
                let value = conic::logit(v);
                self.push(
                    owner,
                    EqKind::Lin {
                        terms: vec![(slot, 1.0)],
                        rhs: value,
                    },
                    false,
                    Some(Sub::Const {
                        slot,
                        value,
                        abs: value,
                    }),
                )
            }
        }
    }
}

/// Each point's coincidence class: the smallest point index coincident
/// with it through the sketch's coincidences.
fn coincidence_classes(n_points: usize, cons: &[Con]) -> Vec<usize> {
    fn find(uf: &mut [usize], mut x: usize) -> usize {
        while uf[x] != x {
            uf[x] = uf[uf[x]];
            x = uf[x];
        }
        x
    }
    let mut uf: Vec<usize> = (0..n_points).collect();
    for c in cons {
        if let Con::Coincident(p, q) = *c {
            let (rp, rq) = (find(&mut uf, p), find(&mut uf, q));
            uf[rp.max(rq)] = rp.min(rq);
        }
    }
    (0..n_points).map(|p| find(&mut uf, p)).collect()
}

/// Builds the residual system with the given scaling.
pub(crate) fn build<K: SketchKey>(index: &Index<K>, cons: &[Con], scale: Scale) -> Problem {
    let n_points = index.points.len();
    let n_circles = index.circles.len();
    let n_shapes = index.conics.len();
    let splines: Vec<SplineDef> = index
        .splines
        .iter()
        .map(|s| SplineDef {
            map: s.map.clone(),
            points: s.points.clone(),
        })
        .collect();
    let mut b = Builder {
        scale,
        radius_base: 2 * n_points,
        shape_base: 2 * n_points + n_circles,
        param_base: 2 * n_points + n_circles + n_shapes,
        splines: &splines,
        eqs: Vec::new(),
        params: Vec::new(),
    };
    // Implicit equations first: every arc keeps |end − centre| = |start −
    // centre|, except derived images, whose source's equation implies it.
    let derived_arcs: std::collections::BTreeSet<usize> = cons
        .iter()
        .flat_map(|c| match c {
            Con::Mirror(_, d) => d.arcs.clone(),
            Con::Offset(o) => o.arcs.clone(),
            Con::LinearPattern(_, pt) | Con::CircularPattern(_, _, pt) => pt.arcs.clone(),
            _ => Vec::new(),
        })
        .collect();
    for (i, k) in index.kinds.iter().enumerate() {
        if derived_arcs.contains(&i) {
            continue;
        }
        if let crate::index::EntKind::Circular(Circ::Arc { c, s, e }) = *k {
            b.eqs.push(Equation {
                kind: EqKind::Norm {
                    k: 1.0,
                    v: VecExpr::diff(e, c),
                    rhs: Rhs::term(1.0, RTerm::Norm(VecExpr::diff(s, c))),
                },
                owner: None,
                arc: Some(i),
                angular: false,
                sub: None,
            });
        }
    }
    // Every ellipse keeps its axes perpendicular and every elliptical arc
    // its ends on its ellipse (design note 0008 §4.1–4.2), except derived
    // images, whose source's equations imply theirs. An end coincident
    // with a point already on the same ellipse lies there through that
    // point: an axis end (on it by definition), or an end of another arc
    // of it, as where the pieces of a split arc meet (0008 §8). The
    // ellipse is the same when its three points are, up to coincidence.
    let derived_ellipses: std::collections::BTreeSet<usize> = cons
        .iter()
        .flat_map(|c| match c {
            Con::Mirror(_, d) => d.ellipses.clone(),
            Con::LinearPattern(_, pt) | Con::CircularPattern(_, _, pt) => pt.ellipses.clone(),
            _ => Vec::new(),
        })
        .collect();
    let class = coincidence_classes(n_points, cons);
    let mut held = std::collections::BTreeSet::new();
    for (i, k) in index.kinds.iter().enumerate() {
        if let crate::index::EntKind::Ellipse(el) = *k {
            if !derived_ellipses.contains(&i) {
                let key = (class[el.c], class[el.m], class[el.n]);
                held.insert((class[el.m], key));
                held.insert((class[el.n], key));
            }
        }
    }
    for (i, k) in index.kinds.iter().enumerate() {
        let crate::index::EntKind::Ellipse(el) = *k else {
            continue;
        };
        if derived_ellipses.contains(&i) {
            continue;
        }
        b.eqs.push(Equation {
            kind: EqKind::DirPair {
                a: VecExpr::diff(el.m, el.c),
                b: VecExpr::diff(el.n, el.c),
                mode: PairMode::Dot,
            },
            owner: None,
            arc: Some(i),
            angular: true,
            sub: None,
        });
        if let Some((s, e)) = el.arc {
            let key = (class[el.c], class[el.m], class[el.n]);
            for q in [s, e] {
                if held.insert((class[q], key)) {
                    b.on_ellipse(None, Some(i), q, el);
                }
            }
        }
    }
    // Every text keeps its box a rectangle of its aspect ratio (ADR 0015),
    // except derived images: with u = top-left − bottom-left and R the
    // quarter turn clockwise (anticlockwise when mirrored),
    // bottom-right − bottom-left = aspect·R(u) and top-right − top-left =
    // bottom-right − bottom-left. Linear, so exact in one step.
    let derived_texts: std::collections::BTreeSet<usize> = cons
        .iter()
        .flat_map(|c| match c {
            Con::Mirror(_, d) => d.texts.clone(),
            Con::LinearPattern(_, pt) | Con::CircularPattern(_, _, pt) => pt.texts.clone(),
            _ => Vec::new(),
        })
        .collect();
    for (i, k) in index.kinds.iter().enumerate() {
        let crate::index::EntKind::Text(t) = *k else {
            continue;
        };
        if derived_texts.contains(&i) {
            continue;
        }
        let text = index.texts[t];
        let [bl, br, tr, tl] = text.c;
        let ks = if text.mirrored {
            -text.aspect
        } else {
            text.aspect
        };
        let (x, y) = (|p: usize| 2 * p, |p: usize| 2 * p + 1);
        let rows = [
            // br.x − bl.x − ks·(tl.y − bl.y) = 0
            vec![(x(br), 1.0), (x(bl), -1.0), (y(tl), -ks), (y(bl), ks)],
            // br.y − bl.y + ks·(tl.x − bl.x) = 0
            vec![(y(br), 1.0), (y(bl), -1.0), (x(tl), ks), (x(bl), -ks)],
            // tr − tl − br + bl = 0
            vec![(x(tr), 1.0), (x(tl), -1.0), (x(br), -1.0), (x(bl), 1.0)],
            vec![(y(tr), 1.0), (y(tl), -1.0), (y(br), -1.0), (y(bl), 1.0)],
        ];
        for terms in rows {
            b.eqs.push(Equation {
                kind: EqKind::Lin { terms, rhs: 0.0 },
                owner: None,
                arc: Some(i),
                angular: false,
                sub: None,
            });
        }
    }
    for (pos, con) in cons.iter().enumerate() {
        b.constraint(pos, con);
    }
    let Builder { eqs, params, .. } = b;
    Problem {
        n_points,
        n_circles,
        n_shapes,
        params,
        eqs,
        splines,
        degen: DEGENERATE,
    }
}
