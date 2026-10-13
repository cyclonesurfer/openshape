//! Feature evaluation: turns a Part Studio's feature list into sketches and
//! bodies by calling the sketch solver and the geometry kernel (through the
//! [`Kernel`] trait).
//!
//! A failing feature is reported with its error and later features still
//! regenerate. The incremental, cached, cancellable regeneration engine
//! (`regen`) replaces the sequential loop in Phase 3.
//!
//! Every part carries the persistent names of its faces, edges and
//! vertices (design note 0013), and features refer to the parts of earlier
//! features through them: a sketch on a planar face, an extrude up to one
//! (ADR 0019).

pub mod regions;

use std::cell::OnceCell;
use std::f64::consts::TAU;

use opendatum_core_model::naming::{Hint, Name, Reference};
use opendatum_core_model::sketch::{EdgeKey, VertexKey};
use opendatum_core_model::{
    CurveDef, Feature, FeatureId, PartStudio, PlaneRef, ProfileDef, ProfileSource, RegionKey,
    RevolveAxis, SketchDef, SketchId,
};
use opendatum_kernel::{
    BodyMesh, DVec2, DVec3, ExtrudeSpec, FacePlane, Frame, Kernel, KernelError, Origins, Profile,
    ProfileCurve, RevolveSpec, SketchArrangement, TessellationOptions,
};
use opendatum_naming::{BodyNames, BodySource, Made, NamingError, Pick, PickedRegion, Resolved};
use opendatum_sketch_solver::{resolve_curves, solve, ResolvedCurve, SolveOptions, SolveStatus};

pub use opendatum_sketch_solver::{EntityState, Solution};

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
    /// Not regenerated because it is after the rollback bar.
    RolledBack,
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
    /// The names of its faces, edges and vertices (design note 0013), or
    /// why they could not be made; references to the body then fail with
    /// that reason (ADR 0019 §5).
    pub names: Result<BodyNames, NamingError>,
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

fn invalid(m: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(m.into())
}

/// Converts a standard or custom model plane to a kernel frame. A plane on
/// a face needs the parts it refers to: see [`sketch_frame`].
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
        PlaneRef::Face(_) => {
            return Err(invalid("a plane on a face needs the parts made before it"))
        }
    })
}

/// The frame of a sketch or profile plane; a plane on a face is resolved
/// against `bodies`, the parts made before the referring feature, and gets
/// the frame [`face_frame`] gives it.
pub fn sketch_frame<K: Kernel>(
    kernel: &K,
    p: &PlaneRef,
    bodies: &[PartBody<K::Body>],
) -> Result<Frame, KernelError> {
    match p {
        PlaneRef::Face(r) => Ok(face_frame(&reference_plane(
            kernel,
            r,
            bodies,
            "sketch face",
        )?)),
        _ => plane_frame(p),
    }
}

/// The frame of a sketch on a planar face (ADR 0019 §2): its normal is the
/// face's outward normal; its origin is the world origin projected onto the
/// plane; its x axis is world X projected onto it for a face that faces
/// mostly up or down (|n·Z| at least |n·X| and |n·Y|), else Z × n, which
/// keeps the sketch upright. Faces parallel to the standard planes get the
/// Top, Front and Right frames (or their opposites), offset.
pub fn face_frame(plane: &FacePlane) -> Frame {
    let n = plane.normal;
    let origin = n * plane.point.dot(n);
    let x_hint = if n.z.abs() >= n.x.abs().max(n.y.abs()) {
        DVec3::X
    } else {
        DVec3::Z.cross(n)
    };
    Frame::from_normal(origin, n, x_hint)
}

/// The feature whose part a name is of (until booleans, the faces around
/// an edge or vertex are all of one feature).
fn name_feature(name: &Name) -> Option<FeatureId> {
    match name {
        Name::Face(f) => Some(f.feature),
        Name::Edge(e) => e.faces.first().map(|f| f.feature),
        Name::Vertex(v) => v.faces.first().map(|f| f.feature),
    }
}

/// The mesh stored hints are checked against: as fine as the display mesh
/// the application picks on.
fn hint_mesh() -> TessellationOptions {
    TessellationOptions {
        chord_tolerance: 2e-5,
        max_angle: 12f64.to_radians(),
    }
}

/// How far a stored hint may lie from an entity's mesh and still be on it
/// (0.1 mm): the chord deviation of the mesh it was picked on and of
/// [`hint_mesh`], with room to spare. Hints only break ties (design note
/// 0013 §6): a hint this close to two candidates is ambiguous, never
/// resolved to either.
pub const HINT_TOLERANCE: f64 = 1e-4;

/// Distance from `p` to the triangle `abc` (Ericson, *Real-Time Collision
/// Detection*, §5.1.5).
fn triangle_distance(p: DVec3, a: DVec3, b: DVec3, c: DVec3) -> f64 {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.length();
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return bp.length();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return (p - (a + ab * (d1 / (d1 - d3)))).length();
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return cp.length();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return (p - (a + ac * (d2 / (d2 - d6)))).length();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).length();
    }
    let denom = va + vb + vc;
    if denom.abs() < 1e-300 {
        // A degenerate triangle: the nearest of its sides.
        return [(a, b), (b, c), (c, a)]
            .iter()
            .map(|&(s, t)| segment_distance(p, s, t))
            .fold(f64::INFINITY, f64::min);
    }
    let (v, w) = (vb / denom, vc / denom);
    (p - (a + ab * v + ac * w)).length()
}

fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared().max(1e-300)).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

/// Distance from `p` to a kernel entity's mesh (infinite for an index the
/// mesh does not have).
pub fn mesh_distance(mesh: &BodyMesh, pick: Pick, p: DVec3) -> f64 {
    match pick {
        Pick::Face(f) => mesh.faces.get(f.0 as usize).map_or(f64::INFINITY, |fm| {
            fm.triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| fm.positions[i as usize]);
                    triangle_distance(p, a, b, c)
                })
                .fold(f64::INFINITY, f64::min)
        }),
        Pick::Edge(e) => mesh.edges.get(e.0 as usize).map_or(f64::INFINITY, |pl| {
            if let [q] = pl.as_slice() {
                return (p - *q).length();
            }
            pl.windows(2)
                .map(|w| segment_distance(p, w[0], w[1]))
                .fold(f64::INFINITY, f64::min)
        }),
        Pick::Vertex(v) => mesh
            .vertices
            .get(v.0 as usize)
            .map_or(f64::INFINITY, |q| (p - *q).length()),
    }
}

/// Resolves a reference against `bodies`, the parts made before the
/// referring feature (design note 0013 §6): the index of the part it is
/// to, and what it is now. The error completes a sentence about the
/// reference ("the sketch face …").
pub fn resolve_reference<K: Kernel>(
    kernel: &K,
    reference: &Reference,
    bodies: &[PartBody<K::Body>],
) -> Result<(usize, Resolved), String> {
    let feature = name_feature(&reference.name).ok_or("names no feature")?;
    let b = bodies
        .iter()
        .position(|b| b.feature == feature)
        .ok_or("is on a part that is missing, suppressed, failed or later in the list")?;
    let names = bodies[b]
        .names
        .as_ref()
        .map_err(|e| format!("cannot be found: its part's faces cannot be referenced ({e})"))?;
    let mesh: OnceCell<Option<BodyMesh>> = OnceCell::new();
    let lies_on = |pick: Pick, hint: &Hint| {
        mesh.get_or_init(|| kernel.tessellate(&bodies[b].body, &hint_mesh()).ok())
            .as_ref()
            .is_some_and(|m| mesh_distance(m, pick, hint.point) <= HINT_TOLERANCE)
    };
    let resolved = names
        .resolve(reference, &lies_on)
        .map_err(|e| e.to_string())?;
    Ok((b, resolved))
}

/// The plane of the planar face a carrier reference resolves to: all its
/// pieces (a face split since it was picked) must lie in one plane. `what`
/// names the reference in errors.
fn reference_plane<K: Kernel>(
    kernel: &K,
    reference: &Reference,
    bodies: &[PartBody<K::Body>],
    what: &str,
) -> Result<FacePlane, KernelError> {
    let (b, resolved) = resolve_reference(kernel, reference, bodies)
        .map_err(|e| invalid(format!("the {what} {e}")))?;
    let tol = kernel.tolerance();
    let mut plane: Option<FacePlane> = None;
    for pick in resolved.pieces {
        let Pick::Face(f) = pick else {
            return Err(invalid(format!("the {what} is not a face")));
        };
        let fp = kernel
            .face_plane(&bodies[b].body, f)
            .ok_or_else(|| invalid(format!("the {what} is not planar")))?;
        match plane {
            None => plane = Some(fp),
            Some(q) => {
                if q.normal.cross(fp.normal).length() > tol.angular
                    || q.normal.dot(fp.normal) < 0.0
                    || (fp.point - q.point).dot(q.normal).abs() > tol.linear
                {
                    return Err(invalid(format!(
                        "the {what} is now several faces in different planes"
                    )));
                }
            }
        }
    }
    plane.ok_or_else(|| invalid(format!("the {what} no longer exists")))
}

/// The signed distance along `frame`'s normal from its plane to the plane
/// of the face `reference` resolves to, which must be parallel to it
/// (ADR 0019 §3).
fn distance_to_face<K: Kernel>(
    kernel: &K,
    frame: &Frame,
    reference: &Reference,
    bodies: &[PartBody<K::Body>],
) -> Result<f64, KernelError> {
    let what = "face to extrude up to";
    let plane = reference_plane(kernel, reference, bodies, what)?;
    let tol = kernel.tolerance();
    if plane.normal.cross(frame.z).length() > tol.angular {
        return Err(invalid(format!(
            "the {what} is not parallel to the sketch plane (other faces need surface \
             intersection, not yet available)"
        )));
    }
    let d = (plane.point - frame.origin).dot(frame.z);
    if d.abs() <= tol.linear {
        return Err(invalid(format!("the {what} is in the sketch plane")));
    }
    Ok(d)
}

/// Converts an inline model profile on a standard or custom plane to a
/// kernel profile (see [`profile_in`] for one on a face).
pub fn profile(def: &ProfileDef) -> Result<Profile, KernelError> {
    curves_on(def, plane_frame(&def.plane)?)
}

/// Converts an inline model profile to a kernel profile, resolving a plane
/// on a face against `bodies` (the parts made before the feature).
pub fn profile_in<K: Kernel>(
    kernel: &K,
    def: &ProfileDef,
    bodies: &[PartBody<K::Body>],
) -> Result<Profile, KernelError> {
    curves_on(def, sketch_frame(kernel, &def.plane, bodies)?)
}

fn curves_on(def: &ProfileDef, plane: Frame) -> Result<Profile, KernelError> {
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

/// Converts a solved sketch curve to a kernel curve. Ellipses and
/// elliptical arcs become exact rational quadratic splines
/// ([`ProfileCurve::elliptical_arc`]), a full ellipse running
/// counter-clockwise from its major end; so do conics
/// ([`ProfileCurve::conic`]). Fails only for a degenerate ellipse or conic,
/// which [`resolve_curves`] never returns.
pub fn kernel_curve(c: &ResolvedCurve) -> Result<ProfileCurve, KernelError> {
    Ok(match c {
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
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        } => {
            // θ increases counter-clockwise when the minor axis is left of
            // the major one.
            let sweep = if major.perp_dot(*minor) < 0.0 {
                -TAU
            } else {
                TAU
            };
            ProfileCurve::elliptical_arc(*center, *major, *minor, 0.0, sweep)?
        }
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            start_param,
            sweep_param,
        } => ProfileCurve::elliptical_arc(*center, *major, *minor, *start_param, *sweep_param)?,
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => ProfileCurve::conic(*start, *control, *end, *weight)?,
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
    })
}

/// Solves a sketch feature and computes its regions. `bodies` are the
/// parts made before it (for a sketch on a face).
pub fn evaluate_sketch<K: Kernel>(
    kernel: &K,
    feature: FeatureId,
    def: &SketchDef,
    bodies: &[PartBody<K::Body>],
) -> (SketchResult, FeatureStatus) {
    let solution = solve(&def.sketch, &SolveOptions::default());
    let frame = match sketch_frame(kernel, &def.plane, bodies) {
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
    let mut owners = Vec::new();
    let mut kcurves = Vec::new();
    for (id, _, c) in curves.iter().filter(|(_, construction, _)| !construction) {
        owners.push(*id);
        kcurves.push(kernel_curve(c).map_err(|e| e.to_string())?);
    }
    let arrangement = kernel.sketch_regions(&kcurves).map_err(|e| e.to_string())?;
    let keys = regions::region_keys(&arrangement, &owners);
    Ok(SketchRegions {
        curves: kcurves,
        owners,
        arrangement,
        keys,
    })
}

/// A sweep's profile and what its curves were made from.
struct Swept<'a> {
    profile: Profile,
    /// The index of the profile's sketch among the regenerated sketches.
    sketch: Option<usize>,
    /// The sketch edge of each profile curve; made-up ones for an inline
    /// profile (see [`inline_keys`]).
    keys: Vec<EdgeKey>,
    /// Each selected region as the feature stores it and as it is now
    /// (none for an inline profile).
    regions: Vec<(&'a RegionKey, &'a RegionKey)>,
}

impl Swept<'_> {
    /// Names the body swept from the profile.
    fn names<K: Kernel>(
        &self,
        kernel: &K,
        feature: FeatureId,
        body: &K::Body,
    ) -> Result<BodyNames, NamingError> {
        let origins = kernel.origins(body);
        let inline;
        let regions: Vec<PickedRegion> = if self.sketch.is_some() {
            self.regions
                .iter()
                .map(|&(stored, current)| PickedRegion { stored, current })
                .collect()
        } else {
            inline = inline_regions(origins, &self.keys)?;
            inline
                .iter()
                .map(|k| PickedRegion {
                    stored: k,
                    current: k,
                })
                .collect()
        };
        let source = BodySource {
            feature,
            made: Made::Sweep {
                curves: &self.keys,
                regions: &regions,
            },
        };
        BodyNames::new(origins, &kernel.adjacency(body), &source)
    }
}

/// Made-up sketch edges for the curves of an inline profile of `feature`:
/// curve `i` is an edge of its own entity, numbered `i` in the feature's
/// op (ADR 0019 §4), so the names stay as long as the curve order does.
fn inline_keys(feature: FeatureId, n: usize) -> Vec<EdgeKey> {
    (0..n)
        .map(|i| {
            let entity = SketchId {
                op: feature.0,
                index: i as u32,
            };
            let end = |ordinal| VertexKey {
                entities: vec![entity],
                ordinal,
                parts: Vec::new(),
            };
            EdgeKey {
                entity,
                forward: true,
                from: end(0),
                to: end(1),
                part: 0,
            }
        })
        .collect()
}

/// The regions of an inline profile, as the kernel found them, keyed by
/// the made-up edges of their loops (each loop starting at its smallest).
fn inline_regions(origins: &Origins, keys: &[EdgeKey]) -> Result<Vec<RegionKey>, NamingError> {
    origins
        .regions
        .iter()
        .map(|loops| {
            let loops = loops
                .iter()
                .map(|pieces| {
                    let mut l: Vec<EdgeKey> = Vec::new();
                    for p in pieces {
                        let k = keys
                            .get(p.curve as usize)
                            .ok_or(NamingError::MissingCurve(p.curve))?;
                        if l.last() != Some(k) {
                            l.push(k.clone());
                        }
                    }
                    while l.len() > 1 && l.first() == l.last() {
                        l.pop();
                    }
                    if let Some(i) = (0..l.len()).min_by(|&a, &b| l[a].cmp(&l[b])) {
                        l.rotate_left(i);
                    }
                    Ok(l)
                })
                .collect::<Result<Vec<_>, NamingError>>()?;
            Ok(RegionKey {
                loops,
                interior: DVec2::ZERO,
            })
        })
        .collect()
}

/// Builds the kernel profile of a [`ProfileSource`] of `feature`.
fn source_profile<'a, K: Kernel>(
    kernel: &K,
    feature: FeatureId,
    source: &'a ProfileSource,
    sketches: &'a [SketchResult],
    bodies: &[PartBody<K::Body>],
) -> Result<Swept<'a>, KernelError> {
    match source {
        ProfileSource::Inline(def) => {
            let profile = profile_in(kernel, def, bodies)?;
            let keys = inline_keys(feature, profile.curves.len());
            Ok(Swept {
                profile,
                sketch: None,
                keys,
                regions: Vec::new(),
            })
        }
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
            let mut selected = Vec::new();
            for key in regions {
                let i = regions::resolve_region(kernel, sr, key).ok_or_else(|| {
                    KernelError::InvalidInput("a profile region no longer exists".into())
                })?;
                if !picked.contains(&i) {
                    picked.push(i);
                    selected.push((key, &sr.keys[i]));
                }
            }
            let tagged = kernel.region_profile_edges(&sr.arrangement, &sr.curves, &picked)?;
            let edges: Vec<_> = tagged.iter().map(|(_, e)| *e).collect();
            Ok(Swept {
                profile: Profile {
                    plane: s.frame,
                    curves: tagged.into_iter().map(|(c, _)| c).collect(),
                },
                sketch: Some(index),
                keys: regions::edge_keys(&sr.arrangement, &sr.owners, &edges),
                regions: selected,
            })
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

/// Evaluates the solid feature `id`, naming its part. `sketches` are the
/// sketches and `bodies` the parts regenerated so far (earlier in the
/// feature list).
pub fn evaluate_feature<K: Kernel>(
    kernel: &K,
    id: FeatureId,
    feature: &Feature,
    sketches: &[SketchResult],
    bodies: &[PartBody<K::Body>],
) -> Result<PartBody<K::Body>, KernelError> {
    let part = |body: K::Body, names: Result<BodyNames, NamingError>| PartBody {
        feature: id,
        body,
        names,
    };
    match feature {
        Feature::Sketch(_) => Err(KernelError::InvalidInput(
            "a sketch is not a solid feature".into(),
        )),
        Feature::Box { corner, size } => {
            let body = kernel.make_box(&Frame::WORLD.translated(*corner), *size)?;
            let source = BodySource {
                feature: id,
                made: Made::Box {
                    z_positive: size.z > 0.0,
                },
            };
            let names = BodyNames::new(kernel.origins(&body), &kernel.adjacency(&body), &source);
            Ok(part(body, names))
        }
        Feature::Extrude {
            profile: p,
            distance,
            symmetric,
            up_to,
        } => {
            let swept = source_profile(kernel, id, p, sketches, bodies)?;
            let (distance, symmetric) = match up_to {
                None => (*distance, *symmetric),
                Some(_) if *symmetric => {
                    return Err(invalid("a symmetric extrude cannot end at a face"))
                }
                Some(r) => (
                    distance_to_face(kernel, &swept.profile.plane, r, bodies)?,
                    false,
                ),
            };
            let body = kernel.extrude(
                &swept.profile,
                &ExtrudeSpec {
                    distance,
                    symmetric,
                },
            )?;
            let names = swept.names(kernel, id, &body);
            Ok(part(body, names))
        }
        Feature::Revolve {
            profile: p,
            axis,
            angle,
        } => {
            let swept = source_profile(kernel, id, p, sketches, bodies)?;
            let (axis_origin, axis_dir) = revolve_axis(axis, swept.sketch.map(|i| &sketches[i]))?;
            let body = kernel.revolve(
                &swept.profile,
                &RevolveSpec {
                    axis_origin,
                    axis_dir,
                    angle: *angle,
                },
            )?;
            let names = swept.names(kernel, id, &body);
            Ok(part(body, names))
        }
    }
}

/// Regenerates a Part Studio from scratch, up to its rollback bar; the
/// features after the bar are reported as [`FeatureStatus::RolledBack`].
pub fn regenerate<K: Kernel>(kernel: &K, studio: &PartStudio) -> StudioResult<K::Body> {
    let mut out = StudioResult {
        statuses: Vec::with_capacity(studio.features.len()),
        bodies: Vec::new(),
        sketches: Vec::new(),
    };
    let bar = studio.rollback_index();
    for (i, f) in studio.features.iter().enumerate() {
        if i >= bar {
            out.statuses.push((f.id, FeatureStatus::RolledBack));
            continue;
        }
        if f.suppressed {
            out.statuses.push((f.id, FeatureStatus::Suppressed));
            continue;
        }
        if let Feature::Sketch(def) = &f.feature {
            let (result, status) = evaluate_sketch(kernel, f.id, def, &out.bodies);
            out.sketches.push(result);
            out.statuses.push((f.id, status));
            continue;
        }
        match evaluate_feature(kernel, f.id, &f.feature, &out.sketches, &out.bodies) {
            Ok(part) => {
                let status = match &part.names {
                    Ok(_) => FeatureStatus::Ok,
                    Err(e) => {
                        FeatureStatus::Warning(format!("its faces cannot be referenced: {e}"))
                    }
                };
                out.statuses.push((f.id, status));
                out.bodies.push(part);
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
