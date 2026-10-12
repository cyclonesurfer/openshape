//! Curvature combs on sketch curves: at points evenly spaced along a curve,
//! a tooth along the normal, as long as the curvature there and on the
//! side away from the centre of curvature, with the tips joined. A comb
//! shows at a glance where a curve bends more or less, where it changes
//! side (an inflection, where the teeth cross over), and whether two
//! curves continue each other's curvature at a joint (their teeth there
//! are equally long: G2, design note 0010) or only their direction (the
//! teeth jump).
//!
//! The samples come from the kernel's exact first and second derivatives
//! ([`Kernel::curvature_samples`]); this module only scales them into
//! teeth for drawing.

use openshape_core_model::sketch::Geometry;
use openshape_core_model::{DVec2, SketchId};
use openshape_features::kernel_curve;
use openshape_kernel::{CurvatureSample, Kernel};
use openshape_sketch_solver::ResolvedCurve;

/// Teeth along each combed curve (one more sample than intervals).
pub const TEETH: usize = 96;

/// Screen length of the longest tooth among the combs shown, in points.
pub const LONGEST_TOOTH: f32 = 72.0;

/// A curvature comb in plane coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct Comb {
    /// Each tooth, from the curve to its tip.
    pub teeth: Vec<[DVec2; 2]>,
    /// The line through the tips, in order along the curve.
    pub envelope: Vec<DVec2>,
}

/// True for geometry that has a curvature to show: circles, arcs,
/// ellipses, elliptical arcs, conics and splines (not points, lines or
/// text).
pub fn has_curvature(g: &Geometry<SketchId>) -> bool {
    matches!(
        g,
        Geometry::Circle { .. }
            | Geometry::Arc { .. }
            | Geometry::Ellipse { .. }
            | Geometry::EllipticalArc { .. }
            | Geometry::Conic { .. }
            | Geometry::ControlSpline { .. }
            | Geometry::InterpolatedSpline { .. }
    )
}

/// Curvature samples of a solved sketch curve ([`TEETH`] intervals), or
/// none if the kernel cannot evaluate it (a degenerate curve).
pub fn samples(kernel: &impl Kernel, curve: &ResolvedCurve) -> Vec<CurvatureSample> {
    kernel_curve(curve)
        .and_then(|c| kernel.curvature_samples(&c, TEETH))
        .unwrap_or_default()
}

/// The largest |curvature| among `samples` (0 for none).
pub fn peak<'a>(samples: impl IntoIterator<Item = &'a CurvatureSample>) -> f64 {
    samples
        .into_iter()
        .map(|s| s.curvature.abs())
        .filter(|k| k.is_finite())
        .fold(0.0, f64::max)
}

/// The comb of `samples` with teeth `scale` long per unit of curvature
/// (plane length × metres): each tooth leaves its point along the normal
/// on the side away from the centre of curvature, |κ|·`scale` long, so a
/// counter-clockwise turn (κ > 0, centre on the left) puts it on the
/// right.
pub fn comb(samples: &[CurvatureSample], scale: f64) -> Comb {
    let teeth: Vec<[DVec2; 2]> = samples
        .iter()
        .map(|s| [s.point, s.point - s.tangent.perp() * (s.curvature * scale)])
        .collect();
    let envelope = teeth.iter().map(|t| t[1]).collect();
    Comb { teeth, envelope }
}

#[cfg(test)]
mod tests {
    use super::*;
    use openshape_kernel::{OsKernel, ProfileCurve};

    fn kernel_samples(c: &ProfileCurve) -> Vec<CurvatureSample> {
        OsKernel::default().curvature_samples(c, 8).unwrap()
    }

    #[test]
    fn teeth_stand_outside_a_counter_clockwise_circle() {
        let (centre, r) = (DVec2::new(0.01, 0.02), 0.005);
        let s = kernel_samples(&ProfileCurve::Circle {
            center: centre,
            radius: r,
        });
        let c = comb(&s, 1e-5);
        assert_eq!(c.teeth.len(), 9);
        for t in &c.teeth {
            // Radially outward, |κ|·scale = (1/r)·1e-5 = 2 mm long.
            let out = (t[0] - centre).normalize();
            assert!(((t[1] - t[0]) - out * 0.002).length() < 1e-12, "{t:?}");
        }
        assert_eq!(c.envelope, c.teeth.iter().map(|t| t[1]).collect::<Vec<_>>());
    }

    #[test]
    fn teeth_stand_outside_a_clockwise_arc_too() {
        let s = kernel_samples(&ProfileCurve::Arc {
            center: DVec2::ZERO,
            radius: 0.01,
            start_angle: 0.5,
            sweep_angle: -2.0,
        });
        let c = comb(&s, 1e-4);
        for t in &c.teeth {
            let out = t[0].normalize();
            assert!(((t[1] - t[0]) - out * 0.01).length() < 1e-12, "{t:?}");
        }
    }

    #[test]
    fn a_straight_line_has_no_teeth_to_speak_of() {
        let s = kernel_samples(&ProfileCurve::Line {
            start: DVec2::ZERO,
            end: DVec2::new(0.03, 0.04),
        });
        let c = comb(&s, 1.0);
        assert!(c.teeth.iter().all(|t| t[0] == t[1]));
        assert_eq!(peak(&s), 0.0);
    }

    #[test]
    fn the_peak_is_the_tightest_bend() {
        let s = kernel_samples(&ProfileCurve::Circle {
            center: DVec2::ZERO,
            radius: 0.004,
        });
        assert!((peak(&s) - 250.0).abs() < 1e-9);
        assert_eq!(peak(&[]), 0.0);
    }
}
