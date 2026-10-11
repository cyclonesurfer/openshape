//! Feature evaluation: turns a Part Studio's feature list into sketches and
//! bodies by calling the sketch solver and the geometry kernel (through the
//! [`Kernel`] trait).
//!
//! A failing feature is reported with its error and later features still
//! regenerate. The incremental, cached, cancellable regeneration engine
//! (`regen`) replaces the sequential loop in Phase 3.

pub mod regions;

use openshape_core_model::{
    CurveDef, Feature, FeatureId, PartStudio, PlaneRef, ProfileDef, ProfileSource, RegionKey,
    RevolveAxis, SketchDef, SketchId,
};
use openshape_kernel::{
    DVec2, DVec3, ExtrudeSpec, Frame, Kernel, KernelError, Profile, ProfileCurve, RevolveSpec,
    SketchArrangement,
};
use openshape_sketch_solver::{resolve_curves, solve, ResolvedCurve, SolveOptions, SolveStatus};

pub use openshape_sketch_solver::{EntityState, Solution};

/// The outcome of one feature.
#[derive(Clone, Debug, PartialEq)]
pub enum FeatureStatus {
    /// Regenerated successfully.
    Ok,
    /// Regenerated, with a problem the user should fix (e.g. redundant
    /// sketch constraints).
    Warning(String),
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

/// A regenerated sketch.
#[derive(Clone, Debug)]
pub struct SketchResult {
    /// The sketch feature.
    pub feature: FeatureId,
    /// The sketch plane.
    pub frame: Frame,
    /// The solver's result (solved geometry and diagnostics).
    pub solution: Solution<SketchId>,
    /// Every curve with its explicit geometry and construction flag.
    pub curves: Vec<(SketchId, bool, ResolvedCurve)>,
    /// Regions of the non-construction curves, or why they could not be
    /// computed.
    pub regions: Result<SketchRegions, String>,
}

/// The regions of a sketch.
#[derive(Clone, Debug)]
pub struct SketchRegions {
    /// The kernel curves, one per non-construction sketch curve.
    pub curves: Vec<ProfileCurve>,
    /// The sketch curve behind each kernel curve.
    pub owners: Vec<SketchId>,
    /// The planar arrangement.
    pub arrangement: SketchArrangement,
    /// The name of each region, in arrangement order.
    pub keys: Vec<RegionKey>,
}

/// The result of regenerating a Part Studio.
#[derive(Clone, Debug)]
pub struct StudioResult<B> {
    /// Status per feature, in feature order.
    pub statuses: Vec<(FeatureId, FeatureStatus)>,
    /// Bodies (parts) in creation order.
    pub bodies: Vec<PartBody<B>>,
    /// Sketches in feature order (including failed solves, which keep their
    /// input geometry).
    pub sketches: Vec<SketchResult>,
}

impl<B> Default for StudioResult<B> {
    fn default() -> Self {
        StudioResult {
            statuses: Vec::new(),
            bodies: Vec::new(),
            sketches: Vec::new(),
        }
    }
}

impl<B> StudioResult<B> {
    /// Status of a feature.
    pub fn status(&self, id: FeatureId) -> Option<&FeatureStatus> {
        self.statuses.iter().find(|(f, _)| *f == id).map(|(_, s)| s)
    }

    /// The regenerated sketch of a sketch feature.
    pub fn sketch(&self, id: FeatureId) -> Option<&SketchResult> {
        self.sketches.iter().find(|s| s.feature == id)
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

/// Converts an inline model profile to a kernel profile.
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

/// Converts a solved sketch curve to a kernel curve.
pub fn kernel_curve(c: &ResolvedCurve) -> ProfileCurve {
    match c {
        ResolvedCurve::Line { start, end } => ProfileCurve::Line {
            start: *start,
            end: *end,
        },
        ResolvedCurve::Circle { center, radius } => ProfileCurve::Circle {
            center: *center,
            radius: *radius,
        },
        ResolvedCurve::Arc {
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
        ResolvedCurve::Spline {
            degree,
            knots,
            points,
        } => ProfileCurve::Spline {
            degree: *degree,
            knots: knots.clone(),
            points: points.clone(),
            weights: None,
        },
    }
}

/// Solves a sketch feature and computes its regions.
pub fn evaluate_sketch<K: Kernel>(
    kernel: &K,
    feature: FeatureId,
    def: &SketchDef,
) -> (SketchResult, FeatureStatus) {
    let solution = solve(&def.sketch, &SolveOptions::default());
    let frame = match plane_frame(&def.plane) {
        Ok(f) => f,
        Err(e) => {
            let result = SketchResult {
                feature,
                frame: Frame::top(),
                solution,
                curves: Vec::new(),
                regions: Err(e.to_string()),
            };
            return (result, FeatureStatus::Failed(e.to_string()));
        }
    };
    let mut status = match &solution.status {
        SolveStatus::Solved => FeatureStatus::Ok,
        SolveStatus::Redundant => FeatureStatus::Warning(format!(
            "{} redundant constraint{}",
            solution.redundant.len(),
            if solution.redundant.len() == 1 {
                ""
            } else {
                "s"
            }
        )),
        SolveStatus::Failed(e) => FeatureStatus::Failed(e.to_string()),
    };
    let curves = match resolve_curves(&solution.sketch) {
        Ok(c) => c,
        Err(e) => {
            if !matches!(status, FeatureStatus::Failed(_)) {
                status = FeatureStatus::Failed(e.to_string());
            }
            Vec::new()
        }
    };
    let regions = if matches!(status, FeatureStatus::Failed(_)) {
        Err("the sketch did not solve".to_string())
    } else {
        sketch_regions(kernel, &curves)
    };
    (
        SketchResult {
            feature,
            frame,
            solution,
            curves,
            regions,
        },
        status,
    )
}

fn sketch_regions<K: Kernel>(
    kernel: &K,
    curves: &[(SketchId, bool, ResolvedCurve)],
) -> Result<SketchRegions, String> {
    let (owners, kcurves): (Vec<SketchId>, Vec<ProfileCurve>) = curves
        .iter()
        .filter(|(_, construction, _)| !construction)
        .map(|(id, _, c)| (*id, kernel_curve(c)))
        .unzip();
    let arrangement = kernel.sketch_regions(&kcurves).map_err(|e| e.to_string())?;
    let keys = regions::region_keys(&arrangement, &owners);
    Ok(SketchRegions {
        curves: kcurves,
        owners,
        arrangement,
        keys,
    })
}

/// Builds the kernel profile of a [`ProfileSource`].
fn source_profile<K: Kernel>(
    kernel: &K,
    source: &ProfileSource,
    sketches: &[SketchResult],
) -> Result<(Profile, Option<usize>), KernelError> {
    match source {
        ProfileSource::Inline(def) => Ok((profile(def)?, None)),
        ProfileSource::Sketch { sketch, regions } => {
            let index = sketches
                .iter()
                .position(|s| s.feature == *sketch)
                .ok_or_else(|| {
                    KernelError::InvalidInput(
                        "the profile's sketch is missing, suppressed or later in the list".into(),
                    )
                })?;
            let s = &sketches[index];
            let sr = s
                .regions
                .as_ref()
                .map_err(|e| KernelError::InvalidInput(format!("sketch regions: {e}")))?;
            if regions.is_empty() {
                return Err(KernelError::InvalidInput(
                    "no profile regions selected".into(),
                ));
            }
            let mut picked = Vec::new();
            for key in regions {
                let i = regions::resolve_region(kernel, sr, key).ok_or_else(|| {
                    KernelError::InvalidInput("a profile region no longer exists".into())
                })?;
                if !picked.contains(&i) {
                    picked.push(i);
                }
            }
            let curves = kernel.region_profile(&sr.arrangement, &sr.curves, &picked)?;
            Ok((
                Profile {
                    plane: s.frame,
                    curves,
                },
                Some(index),
            ))
        }
    }
}

fn revolve_axis(
    axis: &RevolveAxis,
    sketch: Option<&SketchResult>,
) -> Result<(DVec2, DVec2), KernelError> {
    match axis {
        RevolveAxis::Inline { origin, dir } => Ok((*origin, *dir)),
        RevolveAxis::SketchLine { line } => {
            let s = sketch.ok_or_else(|| {
                KernelError::InvalidInput("a sketch-line axis needs a sketch profile".into())
            })?;
            match s.curves.iter().find(|(id, _, _)| id == line) {
                Some((_, _, ResolvedCurve::Line { start, end })) => Ok((*start, *end - *start)),
                Some(_) => Err(KernelError::InvalidInput("the axis is not a line".into())),
                None => Err(KernelError::InvalidInput(
                    "the axis line no longer exists".into(),
                )),
            }
        }
    }
}

/// Evaluates one solid feature. `sketches` are the sketches regenerated so
/// far (earlier in the feature list).
pub fn evaluate_feature<K: Kernel>(
    kernel: &K,
    feature: &Feature,
    sketches: &[SketchResult],
) -> Result<K::Body, KernelError> {
    match feature {
        Feature::Sketch(_) => Err(KernelError::InvalidInput(
            "a sketch is not a solid feature".into(),
        )),
        Feature::Box { corner, size } => kernel.make_box(&Frame::WORLD.translated(*corner), *size),
        Feature::Extrude {
            profile: p,
            distance,
            symmetric,
        } => kernel.extrude(
            &source_profile(kernel, p, sketches)?.0,
            &ExtrudeSpec {
                distance: *distance,
                symmetric: *symmetric,
            },
        ),
        Feature::Revolve {
            profile: p,
            axis,
            angle,
        } => {
            let (prof, sketch) = source_profile(kernel, p, sketches)?;
            let (axis_origin, axis_dir) = revolve_axis(axis, sketch.map(|i| &sketches[i]))?;
            kernel.revolve(
                &prof,
                &RevolveSpec {
                    axis_origin,
                    axis_dir,
                    angle: *angle,
                },
            )
        }
    }
}

/// Regenerates a Part Studio from scratch.
pub fn regenerate<K: Kernel>(kernel: &K, studio: &PartStudio) -> StudioResult<K::Body> {
    let mut out = StudioResult {
        statuses: Vec::with_capacity(studio.features.len()),
        bodies: Vec::new(),
        sketches: Vec::new(),
    };
    for f in &studio.features {
        if f.suppressed {
            out.statuses.push((f.id, FeatureStatus::Suppressed));
            continue;
        }
        if let Feature::Sketch(def) = &f.feature {
            let (result, status) = evaluate_sketch(kernel, f.id, def);
            out.sketches.push(result);
            out.statuses.push((f.id, status));
            continue;
        }
        match evaluate_feature(kernel, &f.feature, &out.sketches) {
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
