//! Feature evaluation: turns a Part Studio's feature list into bodies by
//! calling the geometry kernel through the [`Kernel`] trait.
//!
//! Each feature is a pure function of its definition (and, from Phase 3, of
//! the prior state). A failing feature is reported with its error and later
//! features still regenerate. The incremental, cached, cancellable
//! regeneration engine (`regen`) replaces the sequential loop in Phase 3.

use openshape_core_model::{CurveDef, Feature, FeatureId, PartStudio, PlaneRef, ProfileDef};
use openshape_kernel::{
    DVec3, ExtrudeSpec, Frame, Kernel, KernelError, Profile, ProfileCurve, RevolveSpec,
};

/// The outcome of one feature.
#[derive(Clone, Debug, PartialEq)]
pub enum FeatureStatus {
    /// Regenerated successfully.
    Ok,
    /// Skipped because it is suppressed.
    Suppressed,
    /// Failed; the message is shown in the feature tree.
    Failed(String),
}

/// A body produced by a feature.
#[derive(Clone, Debug)]
pub struct PartBody<B> {
    /// The producing feature.
    pub feature: FeatureId,
    /// The body.
    pub body: B,
}

/// The result of regenerating a Part Studio.
#[derive(Clone, Debug)]
pub struct StudioResult<B> {
    /// Status per feature, in feature order.
    pub statuses: Vec<(FeatureId, FeatureStatus)>,
    /// Bodies (parts) in creation order.
    pub bodies: Vec<PartBody<B>>,
}

impl<B> StudioResult<B> {
    /// Status of a feature.
    pub fn status(&self, id: FeatureId) -> Option<&FeatureStatus> {
        self.statuses.iter().find(|(f, _)| *f == id).map(|(_, s)| s)
    }
}

/// Converts a model plane to a kernel frame.
pub fn plane_frame(p: &PlaneRef) -> Result<Frame, KernelError> {
    Ok(match p {
        PlaneRef::Top => Frame::top(),
        PlaneRef::Front => Frame::front(),
        PlaneRef::Right => Frame::right(),
        PlaneRef::Custom {
            origin,
            normal,
            x_dir,
        } => {
            if !(origin.is_finite() && normal.is_finite() && x_dir.is_finite())
                || normal.length() < 1e-12
            {
                return Err(KernelError::InvalidInput("invalid custom plane".into()));
            }
            Frame::from_normal(*origin, *normal, *x_dir)
        }
    })
}

/// Converts a model profile to a kernel profile.
pub fn profile(def: &ProfileDef) -> Result<Profile, KernelError> {
    let plane = plane_frame(&def.plane)?;
    let mut curves = Vec::with_capacity(def.curves.len());
    for c in &def.curves {
        curves.push(match c {
            CurveDef::Line { start, end } => ProfileCurve::Line {
                start: *start,
                end: *end,
            },
            CurveDef::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => ProfileCurve::Arc {
                center: *center,
                radius: *radius,
                start_angle: *start_angle,
                sweep_angle: *sweep_angle,
            },
            CurveDef::Circle { center, radius } => ProfileCurve::Circle {
                center: *center,
                radius: *radius,
            },
            CurveDef::Spline { through } => {
                if through.len() < 2 {
                    return Err(KernelError::InvalidInput(
                        "a spline needs at least two points".into(),
                    ));
                }
                ProfileCurve::interpolate(through)?
            }
        });
    }
    Ok(Profile { plane, curves })
}

/// Evaluates one feature.
pub fn evaluate_feature<K: Kernel>(kernel: &K, feature: &Feature) -> Result<K::Body, KernelError> {
    match feature {
        Feature::Box { corner, size } => kernel.make_box(&Frame::WORLD.translated(*corner), *size),
        Feature::Extrude {
            profile: p,
            distance,
            symmetric,
        } => kernel.extrude(
            &profile(p)?,
            &ExtrudeSpec {
                distance: *distance,
                symmetric: *symmetric,
            },
        ),
        Feature::Revolve {
            profile: p,
            axis_origin,
            axis_dir,
            angle,
        } => kernel.revolve(
            &profile(p)?,
            &RevolveSpec {
                axis_origin: *axis_origin,
                axis_dir: *axis_dir,
                angle: *angle,
            },
        ),
    }
}

/// Regenerates a Part Studio from scratch.
pub fn regenerate<K: Kernel>(kernel: &K, studio: &PartStudio) -> StudioResult<K::Body> {
    let mut out = StudioResult {
        statuses: Vec::with_capacity(studio.features.len()),
        bodies: Vec::new(),
    };
    for f in &studio.features {
        if f.suppressed {
            out.statuses.push((f.id, FeatureStatus::Suppressed));
            continue;
        }
        match evaluate_feature(kernel, &f.feature) {
            Ok(body) => {
                out.statuses.push((f.id, FeatureStatus::Ok));
                out.bodies.push(PartBody {
                    feature: f.id,
                    body,
                });
            }
            Err(e) => out
                .statuses
                .push((f.id, FeatureStatus::Failed(e.to_string()))),
        }
    }
    out
}

/// Convenience: a default box feature in metres.
pub fn default_box() -> Feature {
    Feature::Box {
        corner: DVec3::ZERO,
        size: DVec3::new(0.04, 0.03, 0.02),
    }
}
