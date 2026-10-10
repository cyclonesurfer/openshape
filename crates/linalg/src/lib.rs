//! Deterministic dense linear algebra for OpenShape's solvers (ADR 0010):
//! Householder QR with column pivoting, complete orthogonal decomposition for
//! least-squares and minimum-norm solves, and one-sided Jacobi SVD, in plain
//! scalar Rust with a fixed operation order.
//!
//! * [`Matrix`]: dense, column-major `f64` matrix.
//! * [`QrPivoted`]: `A·P = Q·R`, rank-revealing (Businger & Golub; Golub &
//!   Van Loan §5.4.1), with numerical rank by a relative threshold.
//! * [`CompleteOrthogonal`]: `A·P = Q·[T 0; 0 0]·Z`, giving minimum-norm
//!   least-squares solutions ([`lstsq_min_norm`], [`weighted_min_norm`]) and
//!   orthonormal null-space bases ([`null_space`]).
//! * [`svd`]: one-sided Jacobi SVD (Hestenes; Demmel & Veselić).
//! * [`select_columns_in_order`]: a maximal independent set of columns,
//!   earlier columns first.
//!
//! **Determinism.** Every routine uses only `+ − × ÷ √` in an order fixed by
//! the code: no fused multiply-add, no threads, no run-time SIMD dispatch.
//! The same inputs give bit-identical outputs on every run and every
//! IEEE-754 machine.
//!
//! Entries are expected to lie well inside the floating-point range (about
//! `1e±150`); the solvers feed scaled, order-one values.

mod cod;
mod householder;
mod matrix;
mod qr;
mod select;
mod svd;
mod vector;

pub use cod::{
    lstsq_min_norm, null_space, projection_norm, residual, weighted_min_norm,
    CompleteOrthogonal, WeightedSolution,
};
pub use matrix::Matrix;
pub use qr::QrPivoted;
pub use select::{select_columns_in_order, ColumnSelection};
pub use svd::{svd, Svd};
pub use vector::{axpy, dot, norm2, norm_inf};
