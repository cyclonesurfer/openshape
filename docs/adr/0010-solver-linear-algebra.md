# ADR 0010 — Linear algebra for solvers

Status: accepted (Phase 2)

## Context

The sketch solver (Phase 2), the mate solver (Phase 6) and later the
simulation crates need linear algebra: least-squares and minimum-norm solves,
rank decisions and null spaces. Solver output becomes model geometry, and the
document must regenerate deterministically: the same binary must produce
bit-identical geometry for the same op log on every machine and every run, or
pixel-exact undo/redo, golden tests and persistent naming tie-breaks become
flaky.

## Options

1. **nalgebra** (Apache-2.0). Mature and general. Large dynamic products go
   through `matrixmultiply`, which picks AVX/FMA kernels at run time, so
   results can differ between CPUs running the same binary.
2. **faer** (MIT). Fastest pure-Rust option with dense and sparse
   factorizations. Runtime SIMD dispatch and optional multithreading make
   results depend on the CPU and thread count; the API changes often; it
   pulls in a large dependency tree.
3. **In-house dense kernels** in a small `openshape-linalg` crate:
   Householder QR with column pivoting, complete orthogonal decomposition
   (least-squares and minimum-norm solves), one-sided Jacobi SVD, written in
   plain scalar Rust with a fixed operation order.

## Decision

Option 3 for dense work. The algorithms are textbook (Golub & Van Loan,
*Matrix Computations*, 4th ed., §5.4 rank-revealing QR and complete orthogonal
decomposition, §8.6 Jacobi SVD; Demmel & Veselić 1992 on the accuracy of
one-sided Jacobi), a few hundred lines, and fully under our control. Rust
never fuses `a * b + c` into an FMA unless asked (`mul_add`), so scalar code
gives the same bits on every x86-64 and AArch64 machine. Sketch systems
decompose into small components (design note 0005), so dense O(n³) cost is
acceptable.

Sparse direct solvers (large assemblies, FEA) get their own ADR when needed;
faer is the leading candidate, used with a fixed thread count and SIMD level
where reproducibility matters.

## Consequences

* `crates/linalg` (`openshape-linalg`) with property tests against
  reconstruction identities (A·P = Q·R, A = U·Σ·Vᵀ, orthogonality,
  minimum-norm optimality) on random and rank-deficient matrices.
* Transcendental functions (`sin`, `atan2`, …) still come from the platform's
  libm, which may differ in the last bit between platforms. Derived geometry
  is therefore bit-exact per platform, not across platforms. Nothing may
  depend on cross-platform bit equality of derived geometry (the op log is the
  source of truth; persistent names are topological, not numeric).
