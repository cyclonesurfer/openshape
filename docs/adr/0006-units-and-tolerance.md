# ADR 0006 — Units, numerics and tolerance

Status: accepted (Phase 0)

* All geometry is `f64`. Internal length unit is the **metre**, angles are
  **radians**. The UI shows millimetres by default; documents will carry a
  display-unit preference.
* `core-model::units` provides `Quantity` with a dimension (length, angle,
  dimensionless) and parsing of values with units (`"25 mm"`, `"1.5in"`,
  `"30 deg"`). Full expressions with units (variables, functions) arrive with
  the Variable Studio in Phase 8 and reuse this type.
* The kernel's `Tolerance` is configurable per kernel instance: default
  linear resolution 1e-7 m (0.1 µm), angular 1e-9 rad. Points closer than
  the linear tolerance are coincident; validation allows 100× tolerance for
  "lies on" checks to absorb projection round-off.
* Topology-deciding 2D predicates use Shewchuk's adaptive robust predicates
  (`robust` crate); triangulation uses `spade` (robust CDT).
