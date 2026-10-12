//! Sketch features through regeneration: solving, regions, and extrudes and
//! revolves whose profiles are sketch regions named by `RegionKey`s
//! (ADR 0011, design notes 0005 and 0006).

use std::f64::consts::{PI, TAU};

use openshape_core_model::sketch::{Constraint, ConstraintEntry, Entity, Geometry};
use openshape_core_model::{
    DVec2, Feature, FeatureId, InsertAt, ModelOp, OpId, PlaneRef, ProfileSource, RegionKey,
    ReplicaId, RevolveAxis, Session, SketchData, SketchDef, SketchEdit, SketchId,
};
use openshape_features::{regenerate, regions::resolve_region, FeatureStatus, StudioResult};
use openshape_kernel::{Body, Kernel, OsKernel};
use openshape_sketch_solver::SolveStatus;

const MM: f64 = 1e-3;

fn insert(s: &mut Session, name: &str, f: Feature) -> FeatureId {
    let element = s.first_studio().unwrap();
    FeatureId(s.apply(ModelOp::InsertFeature {
        element,
        at: InsertAt::End,
        name: name.into(),
        feature: f,
    }))
}

fn new_sketch(s: &mut Session, plane: PlaneRef) -> FeatureId {
    insert(
        s,
        "Sketch",
        Feature::Sketch(SketchDef::new(plane, SketchData::default())),
    )
}

/// Builds one edit batch with batch-local ids.
#[derive(Default)]
struct Batch {
    edits: Vec<SketchEdit>,
    next: u32,
}

impl Batch {
    fn id(&mut self) -> SketchId {
        self.next += 1;
        SketchId::new_in_batch(self.next - 1)
    }
    fn entity(&mut self, geometry: Geometry<SketchId>, construction: bool) -> SketchId {
        let id = self.id();
        self.edits.push(SketchEdit::AddEntity(Entity {
            id,
            geometry,
            construction,
        }));
        id
    }
    fn point(&mut self, x: f64, y: f64) -> SketchId {
        self.entity(
            Geometry::Point {
                at: DVec2::new(x, y) * MM,
            },
            false,
        )
    }
    fn line(&mut self, a: SketchId, b: SketchId) -> SketchId {
        self.entity(Geometry::Line { start: a, end: b }, false)
    }
    fn c(&mut self, constraint: Constraint<SketchId>) -> SketchId {
        let id = self.id();
        self.edits.push(SketchEdit::AddConstraint(ConstraintEntry {
            id,
            constraint,
        }));
        id
    }
    /// Applies the batch; returns a function mapping batch ids to real ids.
    fn apply(self, s: &mut Session, sketch: FeatureId) -> impl Fn(SketchId) -> SketchId {
        let element = s.first_studio().unwrap();
        let op = s.apply(ModelOp::EditSketch {
            element,
            feature: sketch,
            edits: self.edits,
        });
        assert!(
            s.document().conflicts.iter().all(|c| c.op != op),
            "{:?}",
            s.document().conflicts
        );
        move |id: SketchId| SketchId {
            op,
            index: id.index,
        }
    }
}

/// Ids of a rectangle's parts.
struct Rect {
    lines: [SketchId; 4],
    width: SketchId,
    height: SketchId,
}

/// A fully defined rectangle: corner fixed at `(x, y)` mm, `w` × `h` mm.
fn rectangle(s: &mut Session, sketch: FeatureId, x: f64, y: f64, w: f64, h: f64) -> Rect {
    let mut b = Batch::default();
    let c = [(x, y), (x + w, y), (x + w, y + h), (x, y + h)];
    let mut ends = Vec::new();
    for i in 0..4 {
        let p = b.point(c[i].0, c[i].1);
        let q = b.point(c[(i + 1) % 4].0, c[(i + 1) % 4].1);
        ends.push((p, q));
    }
    let lines: Vec<SketchId> = ends.iter().map(|(p, q)| b.line(*p, *q)).collect();
    for i in 0..4 {
        b.c(Constraint::Coincident {
            a: ends[i].1,
            b: ends[(i + 1) % 4].0,
        });
    }
    b.c(Constraint::Horizontal { line: lines[0] });
    b.c(Constraint::Vertical { line: lines[1] });
    b.c(Constraint::Horizontal { line: lines[2] });
    b.c(Constraint::Vertical { line: lines[3] });
    b.c(Constraint::Fix {
        point: ends[0].0,
        at: DVec2::new(x, y) * MM,
    });
    let width = b.c(Constraint::Length {
        line: lines[0],
        value: w * MM,
    });
    let height = b.c(Constraint::Length {
        line: lines[1],
        value: h * MM,
    });
    let r = b.apply(s, sketch);
    Rect {
        lines: [r(lines[0]), r(lines[1]), r(lines[2]), r(lines[3])],
        width: r(width),
        height: r(height),
    }
}

fn regen(s: &Session) -> StudioResult<Body> {
    let el = s.document().element(s.first_studio().unwrap()).unwrap();
    regenerate(&OsKernel::default(), &el.studio)
}

fn keys(s: &Session, sketch: FeatureId) -> Vec<RegionKey> {
    regen(s)
        .sketch(sketch)
        .unwrap()
        .regions
        .as_ref()
        .expect("regions")
        .keys
        .clone()
}

fn extrude(s: &mut Session, sketch: FeatureId, regions: Vec<RegionKey>, d: f64) -> FeatureId {
    insert(
        s,
        "Extrude",
        Feature::Extrude {
            profile: ProfileSource::Sketch { sketch, regions },
            distance: d * MM,
            symmetric: false,
        },
    )
}

/// Extrudes every region of the sketch.
fn extrude_all(s: &mut Session, sketch: FeatureId, d: f64) -> FeatureId {
    let k = keys(s, sketch);
    extrude(s, sketch, k, d)
}

fn volume(r: &StudioResult<Body>, f: FeatureId) -> f64 {
    let b = r
        .bodies
        .iter()
        .find(|b| b.feature == f)
        .unwrap_or_else(|| panic!("no body; status {:?}", r.status(f)));
    OsKernel::default()
        .mass_properties(&b.body, 1e-10)
        .unwrap()
        .volume
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs(),
        "{actual} vs {expected}"
    );
}

fn set_value(s: &mut Session, sketch: FeatureId, constraint: SketchId, mm: f64) -> OpId {
    let element = s.first_studio().unwrap();
    s.apply(ModelOp::EditSketch {
        element,
        feature: sketch,
        edits: vec![SketchEdit::SetValue {
            constraint,
            value: mm * MM,
        }],
    })
}

#[test]
fn a_fully_defined_rectangle_extrudes() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    rectangle(&mut s, sk, -20.0, -10.0, 40.0, 20.0);
    let r = regen(&s);
    let sr = r.sketch(sk).unwrap();
    assert_eq!(r.status(sk), Some(&FeatureStatus::Ok));
    assert_eq!(sr.solution.status, SolveStatus::Solved);
    assert_eq!(sr.solution.dof, 0);
    let regions = sr.regions.as_ref().unwrap();
    assert_eq!(regions.keys.len(), 1);
    assert_close(regions.arrangement.regions[0].area, 800.0 * MM * MM);
    let e = extrude_all(&mut s, sk, 10.0);
    let r = regen(&s);
    assert_eq!(r.status(e), Some(&FeatureStatus::Ok));
    assert_close(volume(&r, e), 8000.0 * MM.powi(3));
    let body = &r.bodies.iter().find(|b| b.feature == e).unwrap().body;
    let t = OsKernel::default().topology(body);
    assert_eq!((t.faces, t.edges, t.vertices), (6, 12, 8));
}

#[test]
fn the_extrude_follows_dimension_changes() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Front);
    let rect = rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let key = keys(&s, sk);
    let e = extrude(&mut s, sk, key.clone(), 5.0);
    set_value(&mut s, sk, rect.width, 25.0);
    set_value(&mut s, sk, rect.height, 30.0);
    let r = regen(&s);
    assert_eq!(r.status(e), Some(&FeatureStatus::Ok));
    assert_close(volume(&r, e), 25.0 * 30.0 * 5.0 * MM.powi(3));
    // The name is topological: it still matches exactly.
    assert_eq!(keys(&s, sk)[0].loops, key[0].loops);
}

#[test]
fn unrelated_geometry_does_not_disturb_the_picked_region() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let e = extrude_all(&mut s, sk, 10.0);
    let mut b = Batch::default();
    let c = b.point(100.0, 0.0);
    b.entity(
        Geometry::Circle {
            center: c,
            radius: 5.0 * MM,
        },
        false,
    );
    let _ = b.apply(&mut s, sk);
    let r = regen(&s);
    assert_eq!(
        r.sketch(sk).unwrap().regions.as_ref().unwrap().keys.len(),
        2
    );
    assert_close(volume(&r, e), 8000.0 * MM.powi(3));
}

#[test]
fn a_split_region_resolves_to_the_part_holding_its_interior_point() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    let rect = rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let key = keys(&s, sk);
    let e = extrude(&mut s, sk, key.clone(), 10.0);
    // A vertical line from the bottom edge to the top edge at x = 10 mm.
    let mut b = Batch::default();
    let p = b.point(10.0, 0.0);
    let q = b.point(10.0, 20.0);
    let l = b.line(p, q);
    b.c(Constraint::PointOnCurve {
        point: p,
        curve: rect.lines[0],
    });
    b.c(Constraint::PointOnCurve {
        point: q,
        curve: rect.lines[2],
    });
    b.c(Constraint::Vertical { line: l });
    let _ = b.apply(&mut s, sk);
    let r = regen(&s);
    let regions = r.sketch(sk).unwrap().regions.as_ref().unwrap();
    assert_eq!(regions.keys.len(), 2);
    let i = resolve_region(&OsKernel::default(), regions, &key[0]).expect("resolves");
    let picked = &regions.arrangement.regions[i];
    // The region chosen is the one containing the stored interior point.
    let interior = key[0].interior;
    let expected = if interior.x < 10.0 * MM { 200.0 } else { 600.0 };
    assert_close(picked.area, expected * MM * MM);
    assert_close(volume(&r, e), expected * 10.0 * MM.powi(3));
}

#[test]
fn deleted_regions_fail_clearly() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    let rect = rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let e = extrude_all(&mut s, sk, 10.0);
    let element = s.first_studio().unwrap();
    s.apply(ModelOp::EditSketch {
        element,
        feature: sk,
        edits: vec![SketchEdit::Delete {
            items: vec![rect.lines[2]],
        }],
    });
    let r = regen(&s);
    match r.status(e) {
        Some(FeatureStatus::Failed(m)) => assert!(m.contains("no longer exists"), "{m}"),
        other => panic!("{other:?}"),
    }
    assert!(r.bodies.iter().all(|b| b.feature != e));
}

#[test]
fn a_sketch_that_does_not_solve_fails_its_extrude() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    let rect = rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let e = extrude_all(&mut s, sk, 10.0);
    // A vertical constraint on a horizontal line.
    let mut b = Batch::default();
    b.c(Constraint::Vertical {
        line: rect.lines[0],
    });
    let _ = b.apply(&mut s, sk);
    let r = regen(&s);
    assert!(matches!(r.status(sk), Some(FeatureStatus::Failed(_))));
    match r.status(e) {
        Some(FeatureStatus::Failed(m)) => assert!(m.contains("did not solve"), "{m}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn redundant_constraints_warn_but_still_build() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    let rect = rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let e = extrude_all(&mut s, sk, 10.0);
    // Opposite sides of a rectangle are already equal.
    let mut b = Batch::default();
    b.c(Constraint::Equal {
        a: rect.lines[0],
        b: rect.lines[2],
    });
    let _ = b.apply(&mut s, sk);
    let r = regen(&s);
    match r.status(sk) {
        Some(FeatureStatus::Warning(m)) => assert!(m.contains("redundant"), "{m}"),
        other => panic!("{other:?}"),
    }
    assert_eq!(r.status(e), Some(&FeatureStatus::Ok));
}

#[test]
fn a_suppressed_sketch_fails_its_extrude() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Top);
    rectangle(&mut s, sk, 0.0, 0.0, 40.0, 20.0);
    let e = extrude_all(&mut s, sk, 10.0);
    let element = s.first_studio().unwrap();
    s.apply(ModelOp::SetSuppressed {
        element,
        feature: sk,
        suppressed: true,
    });
    let r = regen(&s);
    assert_eq!(r.status(sk), Some(&FeatureStatus::Suppressed));
    match r.status(e) {
        Some(FeatureStatus::Failed(m)) => assert!(m.contains("suppressed"), "{m}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn an_annulus_between_concentric_circles() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Right);
    let mut b = Batch::default();
    let c1 = b.point(0.0, 0.0);
    let c2 = b.point(0.0, 0.0);
    b.entity(
        Geometry::Circle {
            center: c1,
            radius: 20.0 * MM,
        },
        false,
    );
    b.entity(
        Geometry::Circle {
            center: c2,
            radius: 8.0 * MM,
        },
        false,
    );
    b.c(Constraint::Coincident { a: c1, b: c2 });
    let _ = b.apply(&mut s, sk);
    let r = regen(&s);
    let regions = r.sketch(sk).unwrap().regions.as_ref().unwrap();
    assert_eq!(regions.keys.len(), 2);
    let ring = regions
        .arrangement
        .regions
        .iter()
        .position(|g| !g.holes.is_empty())
        .expect("the annulus has a hole");
    let e = extrude(&mut s, sk, vec![regions.keys[ring].clone()], 7.0);
    let r = regen(&s);
    assert_close(volume(&r, e), PI * (400.0 - 64.0) * 7.0 * MM.powi(3));
    // Picking both regions fills the disc.
    let both = keys(&s, sk);
    let e2 = extrude(&mut s, sk, both, 7.0);
    assert_close(volume(&regen(&s), e2), PI * 400.0 * 7.0 * MM.powi(3));
}

#[test]
fn revolve_about_a_sketch_line() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = new_sketch(&mut s, PlaneRef::Front);
    rectangle(&mut s, sk, 10.0, 0.0, 10.0, 30.0);
    // A construction line along the sketch's y axis.
    let mut b = Batch::default();
    let p = b.point(0.0, -10.0);
    let q = b.point(0.0, 40.0);
    let axis = b.entity(Geometry::Line { start: p, end: q }, true);
    let circle_center = b.point(100.0, 100.0);
    let circle = b.entity(
        Geometry::Circle {
            center: circle_center,
            radius: 1.0 * MM,
        },
        true,
    );
    let r = b.apply(&mut s, sk);
    let (axis, circle) = (r(axis), r(circle));
    let region = keys(&s, sk);
    assert_eq!(region.len(), 1, "construction geometry bounds no region");
    let rev = insert(
        &mut s,
        "Revolve",
        Feature::Revolve {
            profile: ProfileSource::Sketch {
                sketch: sk,
                regions: region.clone(),
            },
            axis: RevolveAxis::SketchLine { line: axis },
            angle: TAU,
        },
    );
    let res = regen(&s);
    assert_eq!(res.status(rev), Some(&FeatureStatus::Ok));
    // Pappus: 2π × centroid radius (15 mm) × area (300 mm²).
    assert_close(volume(&res, rev), TAU * 15.0 * 300.0 * MM.powi(3));
    // An axis that is not a line is refused.
    let bad = insert(
        &mut s,
        "Revolve 2",
        Feature::Revolve {
            profile: ProfileSource::Sketch {
                sketch: sk,
                regions: region,
            },
            axis: RevolveAxis::SketchLine { line: circle },
            angle: TAU,
        },
    );
    match regen(&s).status(bad) {
        Some(FeatureStatus::Failed(m)) => assert!(m.contains("not a line"), "{m}"),
        other => panic!("{other:?}"),
    }
}
