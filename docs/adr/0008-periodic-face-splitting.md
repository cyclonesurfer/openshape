# ADR 0008 — Periodic surfaces are split so no face spans more than π

Status: accepted (Phase 1)

## Context

Faces on periodic surfaces (cylinders, cones, spheres, tori, revolved
splines) can wrap around, which complicates parameter-space algorithms
(tessellation, projection, trimming, booleans): seam edges are used twice by
one face, and parameter values jump by 2π.

## Options

1. Periodic faces with seam edges (Parasolid/ACIS style): fewer faces; every
   parameter-space algorithm needs unwrapping logic.
2. Split periodic faces so each spans at most π in its periodic directions,
   with each face's surface frame starting at its own boundary.

## Decision

Option 2. Circles and arcs over π in profiles are split into ≤ π edges;
closed single-curve profile loops are split in two; revolves are split into
⌈θ/π⌉ pieces (at least two for full revolutions). Periodic parameters use the
range [−π/2, 3π/2), which strictly contains every face's [0, ≤π] domain, so
projection never wraps.

## Consequences

Simple, robust parameter-space code. A full cylinder has two side faces; the
UI and persistent naming present them as one logical face where needed
(naming design, ADR 0003). Booleans may later merge faces on the same
surface where the split is not needed.
