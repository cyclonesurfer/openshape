# ADR 0001 — Geometry kernel: write our own, behind a `Kernel` trait

Status: accepted (Phase 0)

## Context

The kernel is the largest and riskiest part of OpenShape. It must be pure
Rust, GPL-3.0-compatible, f64, robust, and eventually support booleans,
fillets, offsets, shelling, healing and STEP-quality geometry. Candidates:

| | truck | fornjot | curvo |
|---|---|---|---|
| License | Apache-2.0 | 0BSD | MIT |
| Scope | B-rep + NURBS, meshing, booleans (`truck-shapeops`), STEP output | Code-first CAD; B-rep with lines/circles; deliberately narrow, re-scoped repeatedly | NURBS curves/surfaces, interpolation, curve booleans, offsets; no solid B-rep |
| Strengths | Most complete Rust B-rep; good modular crates | Very careful validation culture; simple data model | Clean, well-tested NURBS math (NURBS Book algorithms) |
| Weaknesses for us | Booleans fragile on tangencies/coincident faces; topology generic over geometry makes persistent-naming hooks and tolerance handling awkward; no fillets/offsets; heavy generics in public API | No NURBS; no fillets; API churn; not aimed at feature-based CAD | Not a solid modeller; nalgebra-based types |

All three are license-compatible, so the choice is technical.

## Options

1. **Fork truck** and extend it. Fastest to a first boolean; inherits its
   topology model, tolerance approach and boolean architecture, which we would
   have to rewrite to meet the validation, naming and robustness requirements.
2. **Depend on truck/curvo as libraries.** Leaks their types into our
   features unless wrapped; their release cadence would drive ours.
3. **Write our own kernel from scratch**, using truck, fornjot and curvo (and
   the literature) as references, behind a narrow `Kernel` trait.

## Decision

Option 3. `crates/kernel` is a from-scratch B-rep modeller; free-form math
lives in the kernel-private `crates/nurbs` (Piegl & Tiller algorithms). The
rest of the application uses only the `Kernel` trait and neutral types. A
different kernel (or a truck adapter for comparison testing) can be plugged in
behind the same trait.

## Consequences

* We own the topology data model, so persistent-naming generation tags
  (ADR 0003), tolerance handling (ADR 0006) and validation are designed in.
* More up-front work, mitigated by a golden corpus and property tests from
  day one, and by consulting the reference implementations.
* Every kernel operation validates its result (closed loops, manifold edge
  use, edges on surfaces, Euler–Poincaré, orientation) and fails instead of
  returning a corrupt body.

References: L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed., Springer 1997;
M. Mäntylä, *An Introduction to Solid Modeling*, 1988; C. Hoffmann,
*Geometric and Solid Modeling*, 1989.
