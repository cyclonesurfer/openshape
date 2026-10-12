//! NURBS curves and surfaces for the OpenDatum geometry kernel.
//!
//! This crate is kernel-private: only `opendatum-kernel` depends on it. It
//! implements the core algorithms of *The NURBS Book* (Piegl & Tiller, 2nd
//! edition) in `f64`:
//!
//! - [`KnotVector`]: knot span search (A2.1), basis functions (A2.2) and
//!   their derivatives (A2.3).
//! - [`NurbsCurve`]: rational evaluation (A4.1) and derivatives (A3.2 with
//!   A4.2), knot insertion (A5.1), splitting, reversal, affine transforms,
//!   exact circular arcs (A7.1), global interpolation (A9.1), point inversion
//!   (Eqs. 6.3-6.5), arc length by adaptive Gauss-Legendre quadrature.
//! - [`NurbsSurface`]: rational evaluation (A4.3) and derivatives (A3.6 with
//!   A4.4), normals, point inversion (Eqs. 6.6-6.8), surfaces of revolution
//!   (A8.1), extrusions, isoparametric curves and affine transforms.
//!
//! Control points are stored in homogeneous form `(w x, w y, w z, w)`. All
//! constructors validate their input and return [`NurbsError`] instead of
//! panicking; evaluation clamps parameters to the domain.

mod curve;
mod error;
mod knots;
mod surface;
mod util;

pub use curve::NurbsCurve;
pub use error::NurbsError;
pub use knots::KnotVector;
pub use surface::NurbsSurface;

/// Re-export of the `glam` vector and matrix types used in this crate's API.
pub use glam::{DMat4, DVec2, DVec3, DVec4};
