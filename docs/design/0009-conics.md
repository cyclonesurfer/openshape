# Design note 0009 — Conic arcs

Status: accepted (2026-10-12). Extends design notes 0005 (solver), 0006
(regions), 0007 (editing tools) and 0008 (ellipses, whose kernel form this
reuses).

## 1. Scope

A conic arc sketch entity defined the way lofting has defined conics since
the 1940s: its two end points, the point where the tangents at the ends
meet, and a shape value ρ (rho). Its data, the solver's equations and
diagnostics, the ρ dimension, how the kernel receives it exactly, the
drawing tool, how the editing tools treat it, and the tests.

## 2. Literature

* R. A. Liming, *Practical Analytic Geometry with Applications to
  Aircraft*, Macmillan, 1944: conics through two points with given end
  tangents and a shoulder point, the "rho conic" of lofting practice.
* I. D. Faux & M. J. Pratt, *Computational Geometry for Design and
  Manufacture*, Ellis Horwood, 1979, ch. 2 (conics in design; the
  shoulder point and the conic's discriminating value).
* E. T. Y. Lee, "The rational Bézier representation for conics", in G.
  Farin (ed.), *Geometric Modeling: Algorithms and New Trends*, SIAM,
  1987: every conic arc inside its tangent triangle is a rational
  quadratic Bézier with end weights 1 and middle weight w > 0; elliptic
  for w < 1, parabolic for w = 1, hyperbolic for w > 1; its implicit
  equation in barycentric coordinates.
* Piegl & Tiller, *The NURBS Book*, 2nd ed., §7.2–7.3 (conics as rational
  quadratic Bézier curves; the shoulder point and the middle weight;
  splitting a conic arc and restoring the standard form), §4.2 (affine
  invariance of rational curves).
* G. Farin, *Curves and Surfaces for CAGD*, 5th ed., ch. 12 (conic
  sections; the shoulder tangent is parallel to the chord).
* Design note 0005 §3–4 (unknowns, hidden curve parameters, scaling,
  tangency rows), 0007 §9–§11 (derived images), 0008 §6 (rational
  quadratics in the kernel).

## 3. Data

```text
Geometry::Conic { start, control, end, rho }
```

`start`, `control` and `end` are point entities P₀, P₁, P₂; P₁ is where the
tangents at the ends meet. ρ ∈ (0, 1) is a number stored on the entity
like a circle's radius: document data that the solver treats as an
unknown and writes back after solving (`SketchEdit::SetRho`, as
`SetRadius`).

The curve is the rational quadratic Bézier with weights 1, w, 1,

  C(t) = ((1 − t)² P₀ + 2t(1 − t) w P₁ + t² P₂) / D(t),
  D(t) = (1 − t)² + 2t(1 − t) w + t²,  t ∈ [0, 1],

with w = ρ/(1 − ρ). Its shoulder point C(½) is

  S = M + ρ (P₁ − M),  M = (P₀ + P₂)/2,

so ρ is how far the curve reaches from the chord's midpoint towards P₁
along that median (Liming; P&T §7.3): ρ < ½ elliptic, ρ = ½ parabolic,
ρ > ½ hyperbolic. The arc runs from P₀ (t = 0) to P₂ (t = 1), tangent to
P₁ − P₀ at its start and to P₂ − P₁ at its end, inside the triangle
P₀P₁P₂ (convex hull property, w > 0). It has no orientation convention:
unlike circular and elliptical arcs, nothing depends on which end is the
start.

Degrees of freedom: 6 point coordinates and ρ = 7.

## 4. Equations

### 4.1 The shape unknown

Options: (a) ρ itself; (b) u = ln w = ln(ρ/(1 − ρ)), the logit of ρ.

With (a) a Newton step can leave (0, 1), where w is negative or infinite
and D(t) can vanish: the equations themselves break, not just the
answer. With (b) every real u is a valid conic (w = eᵘ > 0), the
equations are smooth everywhere, and ρ = 1/(1 + e⁻ᵘ).

Recommendation and choice: (b). u is dimensionless and, like hidden
curve parameters, is not scaled by the sketch size (0005 §3).

### 4.2 Points on the conic

A point p on a conic gets a hidden parameter t ∈ [0, 1] (clamped after
each step, as on a spline) and two residual rows C(t) − p, lengths
scaled as all lengths. With R₀ = (1 − t)²/D, R₁ = 2t(1 − t)w/D,
R₂ = t²/D (so C = Σ Rₖ Pₖ and Σ Rₖ = 1):

* ∂C/∂Pₖ = Rₖ I (constant multiples of the identity, as for splines);
* ∂C/∂t = C′(t), the derivative of the rational curve,
  (N′ − C D′)/D with N the numerator;
* ∂C/∂u = w ∂C/∂w = w · 2t(1 − t)(P₁ − C)/D.

Its start value is the parameter of the closest point of the arc, by
sampling, golden-section search and safeguarded Newton as for splines
(0005), and Gauss–Newton re-seeds it the same way when it stalls (0008
§4.2). The samples are 32 uniform in t, then intervals whose chord is
longer than 1/64 of the polyline are halved (up to 2048 samples): a
conic with ρ near 1 runs most of its length within a few hundredths of t
of its ends. When the closest point is an end of the arc, the search
returns that end exactly: golden section alone stops a rounding error
short of it, and Newton's step from there points out of [0, 1]. Net effect: −1 DOF. `PointOnCurve` on a conic means on the arc
t ∈ [0, 1], as on a spline (the full conic continues past the ends; the
arc is the entity).

### 4.3 Tangency at a shared end

`Tangent { kind: AtPoint }` at the start uses the direction P₁ − P₀, at
the end P₂ − P₁, in the existing normalized-cross residual (0005 §4).
Both are linear in the points and independent of ρ, as for a spline's
end tangent.

### 4.4 The ρ dimension

`Rho { curve, value }`: since u is the logit of ρ, a ρ dimension is
the linear row u − logit(v), and like a radius dimension it is
eliminated before Newton's iterations (u is substituted by the constant
logit(v)), so it costs no iteration and has no conditioning of its own.
v must lie in [10⁻⁶, 1 − 10⁻⁶] (§4.6).

After solving, ρ is written back from u. ρ(logit(v)) can differ from v
in the last bit, so the write-back is exact: a conic whose u equals,
bit for bit, the logit of a ρ dimension's value, or of its own stored ρ,
gets that value back; any other gets ρ(u). Dimensioned and untouched
conics therefore keep their ρ exactly, and so do their mirror and
pattern images (§4.5).

### 4.5 Derived images

Mirror and both patterns map a conic to the conic on the images of its
three points with the same ρ: reflections, rotations and translations
are affine and rational curves are affine invariant (P&T §4.2), so the
image of C(t) is the image conic at the same t. The image's u is held
equal to its source's in the derived constraint, as an image circle's
radius is to its source's (0007 §9). No ends swap: a conic has no
orientation. The class of unknowns held equal starts from the source's
value (the dragged member's, if any), never from a mean over the class,
so an image never pulls its source's ρ towards its own.

### 4.6 Collapse

After solving (0005 §9's final check), three kinds of collapse are
reported, never returned as a success, each with the constraint rows
that caused it for the conflict diagnosis:

* closed: its two ends together (the chord has zero length);
* flat: its control point within the solver's degeneracy distance of
  the line through its ends, which includes the control point on an end:
  the curve would be a segment along its chord;
* shape: ρ outside [10⁻⁶, 1 − 10⁻⁶]: the curve would lie within a
  millionth of the triangle's size of its chord or of its control
  polygon.

### 4.7 Not yet

Tangency to a line or a circle away from the ends needs a hidden contact
parameter (the point on the line or circle and the tangent parallel: two
rows, one hidden unknown). As for splines, it is refused until it is
built, and so are sketch fillets that would need it.

## 5. Constraints

Supported: coincident (its points), point on curve (§4.2), tangent at a
shared end (§4.3), ρ (§4.4), mirror and both patterns (§4.5), and
everything on its three points and on lines between them.

Refused before solving, with an error naming the constraint: radius,
diameter, equal, concentric, normal, length, distance to the curve,
offset, and tangency not at a shared end (§4.7).

## 6. Kernel

The rational B-spline of degree 2 with knots 0, 0, 0, 1, 1, 1, control
points P₀, P₁, P₂ and weights 1, w, 1: exactly the curve of §3, in the
form the kernel already handles in profiles, regions, extrude and revolve
(0008 §6). Nothing is approximated.

## 7. Tool

**Conic**: click the start, the end, then the point where the end
tangents meet, then the shoulder. The shoulder is the conic's highest
point above its chord (its tangent there is parallel to the chord), so
the fourth click sets ρ to its height above the chord over the control
point's, ρ = (P₂ − P₀) × (q − P₀) / (P₂ − P₀) × (P₁ − P₀): the curve
reaches up to the pointer. It is kept within [0.02, 0.98] while drawing;
the rubber band shows the conic following the pointer, with the control
polygon and the median. Before the third click the band shows the
parabola (ρ = ½) on the pointer. The third click is refused while it is
within a thousandth of the chord's length of the line through the ends.
Snapped clicks attach their points as every tool does. Makes the three
points and the conic, nothing else; the ρ is free until dimensioned.

* The control polygon is drawn dashed, as a control-point spline's.
* The dimension tool on a conic makes its ρ dimension, typed as a plain
  number; its label reads "ρ 0.42" beside the shoulder point.
* Dragging the curve moves its shoulder with the pointer: ρ changes by
  the pointer's change of height above the chord over the control
  point's height, from wherever the curve was grabbed (so grabbing it
  does not make it jump), within [0.02, 0.98] or the range that already
  held it, and the sketch re-solves, as dragging a circle sets its
  radius. A dimensioned ρ does not move. Dragging a point moves the
  point; a free ρ is kept.
* Snapping onto a conic finds the closest point of the arc, starting
  from the drawn polyline's nearest point.
* A line drawn from a conic's end along its tangent is inferred tangent,
  and a tangent arc can start there.

## 8. Editing tools

* Mirror and both patterns copy conics (§4.5).
* Trim, extend and split: every piece of a conic arc is a conic arc
  (subdivide the homogeneous control points by de Casteljau, then restore
  end weights 1 with w′ = w_mid/√(w₀ w₂), P&T §7.3); as built in §8.1.
* Fillet: refused (§4.7).
* Offset: refused; the offset of a conic is not a conic (Farouki & Neff,
  as for ellipses, 0008 §2).

### 8.1 As built (trim, extend, split)

* **Pieces.** The arc of the conic between parameters a and b, inside the
  arc or past its ends, comes from the polar form (blossom) of the
  quadratic's homogeneous control points H₀ = (P₀, 1), H₁ = (wP₁, w),
  H₂ = (P₂, 1) (Farin, *Curves and Surfaces for CAGD*, ch. 4 and 13):

    b(u, v) = (1 − u)(1 − v) H₀ + ((1 − u) v + u (1 − v)) H₁ + u v H₂,

  the piece's homogeneous points are b(a, a), b(a, b), b(b, b), and its
  end weights are made 1 again with w′ = |w₁|/√(w₀ w₂) (P&T §7.3). It is
  a conic arc of the sketch only when w₀w₂ > 0 and w₀w₁ > 0 (it does not
  pass through a point at infinity or turn half a turn), its ρ =
  w′/(1 + w′) is in range (§4.6), and its control point is more than the
  model tolerance (0.1 µm) from its chord. The last refuses a few microns
  of a gently curved conic, which to the tolerance is a straight segment
  whose ρ means nothing; the fuzz below found the solver flattening such
  a piece altogether. Each refusal is a message, and nothing changes.
* **Parameters.** A cut's parameter is the point's t on the rational form,
  found exactly from its barycentric coordinates (τ₀, τ₁, τ₂) in the
  control triangle. Lee's relation τ₁/τ₀ = 2wt/(1 − t) (§2) gives
  t = τ₁/(τ₁ + 2wτ₀), or 2wτ₂/(τ₁ + 2wτ₂) when |τ₂| > |τ₀|, the form
  that keeps its digits near that end; it holds past the arc's ends too.
  Pieces are ordered, chosen and previewed by t as on a line.
* **Corner and shoulder.** A kept piece's control point moves to the
  piece's own (it stays on the kept end's tangent line, where the new
  end's tangent meets it) and its ρ is set to the piece's. Constraints on
  the old control point and a ρ dimension are removed and counted: they
  described a corner and a shoulder that are no longer there.
* **Points on the conic.** A point-on-curve constraint on the conic goes
  to the piece that still contains its point (a point within the cut gap
  of a piece's end counts as on it), and is removed with a piece that
  goes.
* **Trim.** An end piece moves that end to the cut, held on the cutter
  as for lines (0007 §7.4a). A middle piece leaves a second conic from
  the far cut to the old end, with its own control point and ρ, which
  takes over the old end's constraints. The two are not tied to one conic
  afterwards: they share no defining point, and saying "the same conic"
  needs a constraint the solver does not have yet.
* **Split.** The two pieces are joined by a coincidence and held tangent
  there (`Tangent` at the joint), as they were; their curvatures can
  differ once either is moved.
* **Extend.** An end runs on along the same conic (from the start, as
  the reversed arc does from its end). How far it may run is found by
  a geometric scan of t − 1 from 10⁻⁴ by factors of 1.25 while the arc
  from the kept end to t is a conic arc of the sketch, its control point
  and end stay within four times the sketch's extent of its start, and
  the added arc from the old end to t is one too; the first curve that
  added arc meets is the new end. An ellipse closing on itself or a
  hyperbola running to its asymptote stops the scan; an end that cannot
  grow at all is refused with a message.
* **Tangent touches.** The fuzz's tangent line into a conic's start
  showed the kernel reporting that one touch as three contacts, 0, 1.9
  and 3.7 µm along the line. Extend then moved the line's end to the
  second one and removed its coincidence and tangency, and a later trim
  of the conic kept the 2 µm between. The kernel now finishes such
  contacts at the point of contact (0007 §6).
* **Tests.** The solver tests check pieces against the conic's implicit
  equation, inside the arc and past its ends, with their end points and
  tangents; parameters round-trip for points on the curve; the whole arc
  is its own piece and pieces join; and the limits of growth (an
  ellipse's t < 1/(1 − w), a hyperbola's root of
  1 + 2t(t − 1)(1 − w)). `crates/ui/tests/sketch_conic.rs` drives each
  edit through the editor for elliptic, parabolic and hyperbolic conics:
  every piece on the old conic, segment areas of the regions left against
  the closed form r(w), the constraints handed on and removed, degrees of
  freedom, no redundancy, previews, and the refusals (nothing to extend
  to, a piece too short to keep its shape). Normalising a piece's weight
  wrongly fails four of them. The trim fuzz (`crates/ui/tests/sketch_edit.rs`)
  now also draws random conics, directly and with the conic tool's
  builder with an optional line into its start held tangent there, and
  checks that every conic an edit makes or changes lies on the one it cut
  (nine points each, to 1e-9 m). In a 5,000-case run per property, 2,477
  conic edits committed (1,067 trims, 1,044 splits, 366 extensions),
  every one solved, without redundancy, on the conic it cut, with nothing
  unplanned moving.

## 9. Tests

* Solver: DOF 7, and 6 with a ρ dimension; a point on a conic removes
  one; points on conics checked against Lee's implicit equation in
  barycentric coordinates (τ₀, τ₁, τ₂) of the control triangle,
  τ₁² = 4w² τ₀ τ₂, with all τ ≥ 0, which is independent of the solver's
  parametric form; tangency at both ends with lines and arcs; ρ
  dimensions across elliptic, parabolic and hyperbolic values; images add
  no freedom and no redundancy and keep ρ; a dimension outside its range
  is refused; a conic forced flat is a conflict, not a result; a fuzz of
  random conics and constraints that either solve and pass the
  independent checks or fail with a diagnosis.
* Kernel: the area between the chord and the arc is
  r(w) times the control triangle's area, with

    r(w) = w (arccos w − w √(1 − w²)) / (1 − w²)^{3/2}   (w < 1),
    r(1) = 2/3 (Archimedes),
    r(w) = w (w √(w² − 1) − arccosh w) / (w² − 1)^{3/2}  (w > 1),

  derived for w = cos α from the circular segment of half-angle α
  (triangle area sin³α/cos α, segment α − sin α cos α), continued
  analytically to w > 1, and affine invariant; checked to 1e-12 relative
  and against an independent quadrature. A quarter ellipse as a conic
  (ρ = 1/(1 + √2)) bounds πab/4 with its two semi-axes; a parabolic arc
  revolved about its axis is a paraboloid of volume πr²h/2.
* Editor: the tool by clicking (ρ from the shoulder, refusals), the ρ
  dimension typed, dragging the curve changes ρ, mirror and patterns,
  tangency and inference at the ends, refusals of the editing tools, one
  undo step per edit.

### 9.1 Measured

* Solver fuzz: the generator's conic motif (a random triangle, ρ in
  [0.15, 0.85], optionally a ρ dimension, and one of a point on the arc,
  a line continuing from its end, or a tangent arc from its start) mixed
  into the existing motifs: 20,000 cases per property (hold, stay put,
  perturbed, large) with no failure. Perturbed re-solves, 100,000 seeds
  each, without and with conics: 0 and 0 failures at 1%, 1 and 2 at 2%,
  19 and 19 at 5%; every failure reported, no wrong answer. Conics add
  no measurable failure rate.
* Kernel: segment areas to 1e-12 relative for weights from 0.05 to 6 and
  sizes from millimetres to metres; the closed form agrees with an
  independent quadrature to 1e-10; 2,000 random extruded segments (u in
  [−6, 6]) valid, watertight, area to 1e-11 and volume to 1e-9; the
  paraboloid's and hyperboloid's volume and area to 1e-9.
